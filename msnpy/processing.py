#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import warnings
import operator
import collections
from typing import Sequence
import networkx as nx
import numpy as np
import h5py
from dimspy.portals import mzml_portal
from dimspy.portals import hdf5_portal
from dimspy.portals import thermo_raw_portal
from dimspy.process.replicate_processing import average_replicate_scans
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_filters import filter_attr
from dimspy.process.peak_filters import filter_ringing
from dimspy.process.peak_filters import filter_mz_ranges
from .filters import validate_injection_time_ms1, filter_by_replicates, filter_by_isolation


def hdf5_peaklists_to_txt(filename: str, path_out: str, delimiter:str = "\t"):

    """

    :param filename:
    :param path_out:
    :param delimiter:
    """

    if not os.path.isfile(filename):
        raise IOError('HDF5 database [%s] does not exist' % filename)
    if not h5py.is_hdf5(filename):
        raise IOError('input file [%s] is not a valid HDF5 database' % filename)

    if not os.path.isdir(path_out):
        raise IOError("File or Directory does not exist:".format(path_out))

    obj = hdf5_portal.load_peaklists_from_hdf5(filename)
    if "#" in obj[0].ID:
        fns = set([pl.ID.split("#")[0] for pl in obj])
        sub_ids = [pl.ID.split("#")[1] for pl in obj]
        for fn in fns:
            with open(os.path.join(path_out, os.path.splitext(fn)[0] + ".txt"), "w") as pk_out:
                for i, pl in enumerate(obj):
                    if fn in pl.ID:
                        pl.add_attribute("event", pl.full_shape[0] * [sub_ids[i]], flagged_only=False, on_index=3)
                        str_out = pl.to_str(delimiter=delimiter)
                        if i > 0:
                            pk_out.write(str_out[str_out.index('\n'):])
                        else:
                            pk_out.write(str_out)
                        pl.drop_attribute("event")
    else:
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(delimiter=delimiter))
    return


def mz_tolerance(mz: float, tol: float, unit: str = "ppm"):

    """

    :param mz: mz value
    :param tol: tolerance
    :param unit: ppm or da
    :return:
    :rtype: float
    """

    if unit.lower() == "ppm":
        return mz * (1 - (float(tol) * 0.000001)), mz * (1 + (float(tol) * 0.000001))
    elif unit.lower() == "da":
        return mz - float(tol), mz + float(tol)
    else:
        raise ValueError("Incorrect unit type (options: ppm or da)")


def create_graphs_from_scan_ids(scan_dependents: list, scan_events: dict, ion_injection_times: dict):

    """
    Create Directed Graph from scan dependent relationships

    :param scan_dependents:
    :param scan_events:
    :param ion_injection_times:
    :return:
    :rtype:
    """

    graphs = []
    G = nx.OrderedDiGraph()
    G.add_edges_from(sorted(list(scan_dependents), key=operator.itemgetter(0, 1)))
    for subgraph in [G.subgraph(c) for c in nx.weakly_connected_components(G)]:

        edges = sorted(list(subgraph.edges()), key=operator.itemgetter(0, 1))
        nodes = sorted(subgraph.nodes())

        replicates_within, its =  collections.OrderedDict(), collections.OrderedDict()
        for n in nodes:
            replicates_within.setdefault(scan_events[n], []).append(n)
            its.setdefault(scan_events[n], []).append(ion_injection_times[n])

        G = nx.OrderedDiGraph()
        for rw in replicates_within:

            scan_info = [(None, None, 0.0)]
            scan_info.extend(re.findall(r'([\w\.-]+)@([a-zA-Z]+)(\d+\.\d+)', rw))

            G.add_node(rw,
                       scanids=replicates_within[rw],
                       mslevel=len(scan_info),
                       coltype=scan_info[-1][1],
                       colenergy=float(scan_info[-1][2]),
                       injectiontimes=its[rw],
                       flag=True)
        G.add_edges_from([(scan_events[e[0]], scan_events[e[1]]) for e in edges])
        graphs.append(G)

    return graphs


def merge_ms1_scans(graphs: list):

    """

    :param graphs:
    :return:
    :rtype:
    """

    scan_ids = collections.OrderedDict()
    for G in graphs:
        root = list(nx.topological_sort(G))[0]
        scan_ids.setdefault(root, []).extend(G.node[root]["scanids"])
    for G in graphs:
        root = list(nx.topological_sort(G))[0]
        G.node[root]["scanids"] = scan_ids[root]
    return graphs


def create_templates(graphs: list, nh: int):

    """
    Create a 'master' graph that include all the experimental trees
    Loop through all the subgraphs/graphs

    :param graphs:
    :param nh:
    :return:
    :rtype:
    """

    templates = []
    for G in graphs:
        # Validate if the root node represents a scan event without fragmentation
        # Check if a subgraph, with a user defined number of nodes, exist in the list of templates
        # The nodes (scan events) are matched based on the order they have been collected
        if list(G.edges())[0:nh - 1] not in [list(g.edges())[0:nh - 1] for g in templates]:
            # Create a initial template with a particular number of nodes / edges
            Gt = nx.OrderedDiGraph()
            Gt.add_edges_from(list(G.edges())[0:nh - 1])
            for n in Gt.nodes():
                scan_info = re.findall(r'([\w\.-]+)@([a-zA-Z]+)(\d+\.\d+)', n)
                Gt.node[n]["scanids"] = list()
                Gt.node[n]["mslevel"] = len(scan_info) + 1
                if len(scan_info) == 0:
                    Gt.node[n]["coltype"] = None
                    Gt.node[n]["colenergy"] = None
                else:
                    Gt.node[n]["coltype"] = scan_info[-1][1]
                    Gt.node[n]["colenergy"] = float(scan_info[-1][2])
                Gt.node[n]["template"] = True
                Gt.node[n]["flag"] = True
            templates.append(Gt)
    return templates


def group_by_template(graphs: list, templates: list):

    """

    :param graphs:
    :param templates:
    :return:
    :rtype:
    """

    master_graphs = [G.copy() for G in templates]
    for G in graphs:
        for Gt in templates:
            if G.subgraph(Gt.nodes()).number_of_edges() == Gt.number_of_edges() and \
                          sorted(G.subgraph(Gt.nodes()).nodes()) == sorted(Gt.nodes()):

                i = templates.index(Gt)

                for e in G.edges():
                    for j in range(0, 2):
                        # update master_graphs add nodes/edges or update scanids
                        if e[j] not in master_graphs[i].nodes():
                            master_graphs[i].add_node(e[j],
                                                      scanids=G.node[e[j]]["scanids"],
                                                      mslevel=G.node[e[j]]["mslevel"],
                                                      coltype=G.node[e[j]]["coltype"],
                                                      colenergy=G.node[e[j]]["colenergy"],
                                                      flag=G.node[e[j]]["flag"],
                                                      template=False)
                        else:
                            for scan_id in G.node[e[j]]["scanids"]:
                                if scan_id not in master_graphs[i].node[e[j]]["scanids"]:
                                    master_graphs[i].node[e[j]]["scanids"].append(scan_id)

                    if e not in master_graphs[i].edges():
                        master_graphs[i].add_edge(e[0], e[1])

    return master_graphs


def assign_precursor(peaklist: PeakList, header_frag: str, tolerance: float = 0.5):

    """


    :param peaklist:
    :param header_frag:
    :param tolerance:
    :return:
    :rtype:
    """

    prec_at_energy = re.findall(r'([\w\.-]+)@([\w\.-]+)', header_frag)
    subset = []
    for i, mz in enumerate(peaklist.mz):
        if mz >= float(prec_at_energy[-1][0]) - tolerance and mz <= float(prec_at_energy[-1][0]) + tolerance:
            subset.append((mz, peaklist.intensity[i]))

    if len(subset) > 0:
        s = sorted(subset, key=lambda x: x[1])[-1]
        return s[0], s[1]
    else:
        return None, None


def group_scans(filename: str, nh: int = 2, min_replicates: int = 1, report: str = None,
                max_injection_time: float = None, merge_ms1: bool = False, split: bool = False, remove: bool = True):

    """

    :param filename:
    :param nh:
    :param min_replicates:
    :param report:
    :param max_injection_time:
    :param merge_ms1:
    :param split:
    :param remove:
    :return:
    """

    if filename.lower().endswith(".mzml"):
        d = mzml_portal.Mzml(filename)
    elif filename.lower().endswith(".raw"):
        d = thermo_raw_portal.ThermoRaw(filename)
    else:
        raise IOError("Incorrect file format: {}".format(os.path.basename(filename)))

    si = d.scan_ids()
    sd = d.scan_dependents()
    sit = d.ion_injection_times()

    graphs = create_graphs_from_scan_ids(sd, si, sit)

    for G in list(graphs):
        h = list(nx.topological_sort(G))[0]
        if G.node[h]["mslevel"] > 1:
            warnings.warn("MS1 scan missing. The following scans ids have been removed: {}".format([G.node[n]["scanids"] for n in G.nodes()]))
            graphs.remove(G)

    if max_injection_time:
        for G in list(graphs):
            if not validate_injection_time_ms1(G, max_injection_time):
                scan_id_ms1 = G.nodes[list(nx.topological_sort(G))[0]]["scanids"]
                warnings.warn("Injection time MS1 {} > Maximum injection time for MS1. The following scan ids have been removed: {}".format(scan_id_ms1, [G.node[n]["scanids"] for n in G.nodes()]))
                graphs.remove(G)

    if not split:
        templates = create_templates(graphs, nh)
        groups = group_by_template(graphs, templates)
    else:
        groups = graphs
        for G in groups: nx.set_node_attributes(G, False, 'template')

    for i, G in enumerate(groups):
        G.graph['id'] = i + 1

    if merge_ms1:
        # Merge all MS1 scans across a run/sample
        groups = merge_ms1_scans(groups)

    # flag attribute set to False if not pass filter
    groups = [filter_by_replicates(G, min_replicates) for G in groups]
    groups = [filter_by_isolation(G) for G in groups]

    if report is not None:

        with open(report, "w") as out:
            out.write("tree_id\tevent\ttemplate\tscan_ids\tscans\tflag\n")

            for G in groups:
                if report is not None:
                    for n in G.nodes(data=True):
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(G.graph['id'],
                                                                    n[0],
                                                                    int(n[1]["template"]),
                                                                    n[1]["scanids"],
                                                                    len(n[1]["scanids"]),
                                                                    int(n[1]["flag"])))

    if remove:
        for G in list(groups):
            h = list(nx.topological_sort(G))[0]
            if not G.node[h]["flag"]:
                groups.remove(G)
                continue

            for n in G.nodes(data=True):
                if not n[1]["flag"]:
                    G.remove_node(n[0])
                else:
                    del n[1]['flag']

    if len(groups) == 0:
        warnings.warn("No scan events remaining after filtering. Remove MS data file or alter parameters.")

    return groups


def process_scans(filename: str, groups: list, function_noise: str, snr_thres: float, ppm: float,
                  min_fraction: float = None, rsd_thres: float = None, normalise: bool = False,
                  ringing_thres: float = None, exclusion_list: dict = {}, report: str = None,
                  block_size: int = 5000, ncpus: int = None):

    """

    :param filename:
    :param groups:
    :param function_noise:
    :param snr_thres:
    :param ppm:
    :param min_fraction:
    :param rsd_thres:
    :param normalise:
    :param ringing_thres:
    :param exclusion_list:
    :param report:
    :param block_size: number of peaks in each clustering block.
    :param ncpus: number of CPUs for parallel clustering. Default = None, indicating using as many as possible
    :return: List of (average) PeakList objects (DIMSpy)
    :rtype: Sequence[PeakList]
    """

    print()
    print(os.path.basename(filename))

    if filename.lower().endswith(".mzml"):
        run = mzml_portal.Mzml(filename)
    elif filename.lower().endswith(".raw"):
        run = thermo_raw_portal.ThermoRaw(filename)
    else:
        raise IOError("Incorrect file format: {}".format(os.path.basename(filename)))

    mz_ranges = []
    if exclusion_list is not None and len(exclusion_list) > 0:
        mz_ranges = [mz_tolerance(mz, ppm) for mz in exclusion_list]

    if normalise:
        rsd_on_attr = "intensity_norm"
        rsd_label = "rsd_intensity_norm"
    else:
        rsd_on_attr = "intensity"
        rsd_label = "rsd"

    if report is not None:
        out = open(report, "w")
        out.write("tree_id\tevent\tscans\tpeaks\tmedian_{}\n".format(rsd_label))

    # Check for MS1 scans with the same scan_ids (grouped) to avoid redundant processing
    ms1_headers, temp_scan_ids = collections.OrderedDict(), []
    for G in groups:
        n = list(G.nodes(data=True))[0]
        if temp_scan_ids.count(n[1]["scanids"]) > 1 and n[0] not in ms1_headers:
            ms1_headers[n[0]] = None
        temp_scan_ids.append(n[1]["scanids"])

    pls_avg = []

    for G in groups:
        nodes = G.nodes(data=True)
        print("Processing scans....")
        print("\n".join(map(str, [n[0] for n in nodes])))
        print()
        for n in nodes:
            pls_scans = [run.peaklist(scan_id, function_noise=function_noise) for scan_id in n[1]["scanids"]]
            # Check for MS1 scan available with the same scan_ids (grouped) to avoid redundant processing
            if n[0] in ms1_headers and ms1_headers[n[0]] is not None:
                copy_ms1 = ms1_headers[n[0]].copy()
                # update id
                copy_ms1.ID = "{}#{}:{}".format(os.path.basename(filename), G.graph['id'], n[0])
                pls_avg.append(copy_ms1)
                nscans, n_peaks, median_rsd = len(pls_scans), copy_ms1.shape[0], np.nanmedian(copy_ms1.get_attribute(rsd_label))
            else:
                if ringing_thres is not None and float(ringing_thres) > 0.0:
                    #print "Removing ringing artifacts....."
                    pls_scans = [filter_ringing(pl, threshold=ringing_thres, bin_size=1.0) if len(pl.mz) > 0 else pl for pl in pls_scans]

                pls_scans = [filter_attr(pl, "snr", min_threshold=snr_thres) if len(pl.mz) > 0 else pl for pl in pls_scans]

                if normalise:
                    # print "Normalise by Total Ion Current (TIC)....."
                    pls_scans = [pl.add_attribute("intensity_norm", pl.get_attribute("intensity", False) / pl.metadata["tic"], flagged_only=False, on_index=2) if len(pl.mz) > 0 else pl for pl in pls_scans]

                #print "Aligning, averaging and filtering peaks....."
                nscans, n_peaks, median_rsd = len(pls_scans), 0, "NA"

                if sum(pl.shape[0] for pl in pls_scans) == 0:
                    warnings.warn("No scan data available for {}".format(n[0]))
                else:
                    if len(pls_scans) == 1:
                        pl_avg = average_replicate_scans("{}#{}:{}".format(os.path.basename(filename), G.graph['id'], n[0]), pls_scans, ppm, min_fraction, None, rsd_on_attr, block_size, ncpus)
                        if rsd_on_attr != "intensity":
                            pl_avg.add_attribute("rsd_{}_flag".format(rsd_on_attr), np.ones(pl_avg.full_size), flagged_only=False, is_flag=True)
                        else:
                            pl_avg.add_attribute("rsd_flag", np.ones(pl_avg.full_size), flagged_only=False, is_flag=True)
                    else:
                        pl_avg = average_replicate_scans("{}#{}:{}".format(os.path.basename(filename), G.graph['id'], n[0]), pls_scans, ppm, min_fraction, rsd_thres, rsd_on_attr, block_size, ncpus)

                    if exclusion_list is not None and len(exclusion_list) > 0:
                        pl_avg = filter_mz_ranges(pl_avg, mz_ranges, flag_name="exclusion_flag", flagged_only=False)

                    # add to full_scans to avoid redundant processing
                    if n[0] in ms1_headers and ms1_headers[n[0]] is None:
                        ms1_headers[n[0]] = pl_avg.copy()

                    pls_avg.append(pl_avg)
                    n_peaks, median_rsd = pl_avg.shape[0], np.nanmedian(pl_avg.get_attribute(rsd_label))

            if report is not None:
                out.write("{}\t{}\t{}\t{}\t{}\n".format(groups.index(G) + 1, n[0], nscans, n_peaks, median_rsd))

        if len(pls_avg) == 0:
            raise IOError("No peaks remaining after filtering. Remove file from Study (filelist).")

    if report is not None:
        out.close()

    return pls_avg


def create_spectral_trees(trees: Sequence[nx.OrderedDiGraph], peaklists: Sequence[PeakList]):

    """

    :param trees: list of NetworkX graphs
    :param peaklists: list of PeakList objects
    :return:
    :rtype: Sequence[nx.OrderedDiGraph]
    """

    spectral_trees = []

    headers = [pl.ID.split("#")[1] for pl in peaklists]

    for i, G in enumerate(trees):
        GG = nx.OrderedDiGraph()
        GG.graph["id"] = G.graph["id"]
        for edge in list(G.edges(data=True)):

            header_prec = "{}:{}".format(G.graph["id"], edge[0])
            if len(G.node[edge[0]]["scanids"]) == 0 or header_prec not in headers:
                if " ms " in header_prec:
                    warnings.warn("Cannot create a spectral tree without precursor from {}".format(header_prec))
                    break
                continue

            pl = peaklists[headers.index(header_prec)]
            mz_prec, intensity_prec = assign_precursor(pl, edge[1], tolerance=0.5)

            if not mz_prec:
                if " ms " in header_prec:
                    warnings.warn("Cannot create a spectral tree without precursor from {}".format(header_prec))
                    break
                continue
            else:
                mz_id_prec = "{}_{}_{}".format(round(mz_prec, 4), headers.index(header_prec), np.where(pl.mz == mz_prec)[0][0])
                GG.add_node(mz_id_prec, mz=mz_prec, intensity=intensity_prec, header=header_prec.split(":")[1], mslevel=G.node[edge[0]]["mslevel"], precursor=True)

            header_frag = "{}:{}".format(G.graph["id"], edge[1])
            if len(G.node[edge[1]]["scanids"]) == 0 or header_frag not in headers:
                continue

            pl_fragments = peaklists[headers.index("{}:{}".format(G.graph["id"], edge[1]))]
            for j, mz_frag in enumerate(pl_fragments.mz):
                mz_id_frag = "{}_{}_{}".format(round(mz_frag, 4), headers.index(header_frag), j)
                GG.add_node(mz_id_frag, mz=mz_frag, intensity=pl_fragments.intensity[j], header=header_frag.split(":")[1], mslevel=G.node[edge[1]]["mslevel"], precursor=False)
                GG.add_edge(mz_id_prec, mz_id_frag, mzdiff=round(mz_prec - mz_frag, 7), type="e")

        for node in nx.isolates(GG):
            GG.remove_node(node)

        spectral_trees.append(GG)
    return spectral_trees
