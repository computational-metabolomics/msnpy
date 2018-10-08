
import os
import re
import logging
import operator
import collections
import networkx as nx
import numpy as np
import h5py
from dimspy.portals import mzml_portal
from dimspy.portals import hdf5_portal
from dimspy.portals import thermo_raw_portal
from dimspy.process.replicate_processing import average_replicate_scans
from dimspy.process.peak_filters import filter_attr
from dimspy.process.peak_filters import filter_ringing
from dimspy.process.peak_filters import filter_mz_ranges
from filters import filter_ms1_by_injection_time, filter_by_replicates


def hdf5_peaklists_to_txt(filename, path_out, delimiter="\t"):

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
                        pl.add_attribute("window", pl.full_shape[0] * [sub_ids[i]], flagged_only=False, on_index=3)
                        str_out = pl.to_str(delimiter=delimiter)
                        if i > 0:
                            pk_out.write(str_out[str_out.index('\n'):])
                        else:
                            pk_out.write(str_out)
                        pl.drop_attribute("window")
    else:
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(delimiter=delimiter))
    return


def mz_tolerance(mz, tol, unit="ppm"):
    if unit.lower() == "ppm":
        return mz * (1 - (float(tol) * 0.000001)), mz * (1 + (float(tol) * 0.000001))
    elif unit.lower() == "da":
        return mz - float(tol), mz + float(tol)
    else:
        raise ValueError("Incorrect unit type (options: ppm or da)")


def create_graphs_from_scan_ids(scan_dependents, scan_events, injection_times):
    # Create Directed Graph from scan dependent relationships
    Gs = []
    G = nx.OrderedDiGraph()
    G.add_edges_from(sorted(list(scan_dependents), key=operator.itemgetter(0, 1)))
    for subgraph in nx.weakly_connected_component_subgraphs(G):

        edges = sorted(list(subgraph.edges()), key=operator.itemgetter(0, 1))
        nodes = sorted(subgraph.nodes())

        replicates_within = collections.OrderedDict()
        its = collections.OrderedDict()
        for n in nodes:
            replicates_within.setdefault(scan_events[n], []).append(n)
            its.setdefault(scan_events[n], []).append(injection_times[n])

        G = nx.OrderedDiGraph()
        G.add_nodes_from([(rw, {"scanids": replicates_within[rw], "mslevel": rw.count("@") + 1, "injectiontimes": its[rw]}) for rw in replicates_within])
        G.add_edges_from([(scan_events[e[0]], scan_events[e[1]]) for e in edges])

        h = list(G.nodes())[0]
        if h.count("@") > 0:
            print Warning("Scan dependent missing for {} {}".format(h, nodes))
        else:
            Gs.append(G)

    return Gs


def merge_ms1_scans(Gs):
    scan_ids = collections.OrderedDict()
    for G in Gs:
        edges = list(G.edges())
        scan_ids.setdefault(edges[0][0], []).extend(G.node[edges[0][0]]["scanids"])
    for G in Gs:
        edges = list(G.edges())
        G.node[edges[0][0]]["scanids"] = scan_ids[edges[0][0]]
    return Gs


def create_templates(G_all, nh):
    # Create a 'master' graph that include all the experimental trees
    # Loop through all the subgraphs/graphs

    G_temps = []

    for G in G_all:
        # Validate if the root node represents a scan event without fragmentation
        # Check if a subgraph, with a user defined number of nodes, exist in the list of templates
        # The nodes (scan events) are matched based on the order they have been collected
        if list(G.edges())[0:nh - 1] not in [list(g.edges())[0:nh - 1] for g in G_temps]:
            # Create a initial template with a particular number of nodes / edges
            G_temp = nx.OrderedDiGraph()
            G_temp.add_edges_from(list(G.edges())[0:nh - 1])
            for n in G_temp.nodes():
                G_temp.node[n]["mslevel"] = n.count("@") + 1
                G_temp.node[n]["scanids"] = list()
                G_temp.node[n]["template"] = True
            # Add template to list of templates
            G_temps.append(G_temp)
    return G_temps


def match_templates(G_all, G_temps, nh):
    g_sub = [list(g.edges())[0:nh - 1] for g in G_temps]
    for G in G_all:

        if g_sub.count(list(G.edges())[0:nh - 1]) > 1:
            raise IndexError("Too many templates match")

        i = g_sub.index(list(G.edges())[0:nh - 1])

        for e in G.edges():
            for j in range(0, 2):
                if e[j] not in G_temps[i].nodes():
                    G_temps[i].add_node(e[j], scanids=G.node[e[j]]["scanids"], mslevel=G.node[e[j]]["mslevel"])
                elif G.node[e[j]]["scanids"] not in G_temps[i].node[e[j]]["scanids"]:
                    for scan_id in G.node[e[j]]["scanids"]:
                        if scan_id not in G_temps[i].node[e[j]]["scanids"]:
                            G_temps[i].node[e[j]]["scanids"].append(scan_id)
            if e not in G_temps[i].edges():
                G_temps[i].add_edge(e[0], e[1])
    return G_temps


def match_templates_v2(G_all, G_temps):

    g_sub = [G.copy() for G in G_temps]
    
    for G in G_all:

        for Gt in g_sub:
            if G.subgraph(Gt.nodes()).number_of_edges() == Gt.number_of_edges():

                i = g_sub.index(Gt)

                for e in G.edges():
                    for j in range(0, 2):
                        if e[j] not in G_temps[i].nodes():
                            G_temps[i].add_node(e[j], scanids=G.node[e[j]]["scanids"], mslevel=G.node[e[j]]["mslevel"],
                                                template=False)
                        else:
                            for scan_id in G.node[e[j]]["scanids"]:
                                if scan_id not in G_temps[i].node[e[j]]["scanids"]:
                                    G_temps[i].node[e[j]]["scanids"].append(scan_id)
                    if e not in G_temps[i].edges():
                        G_temps[i].add_edge(e[0], e[1])
    return G_temps


def assign_precursor(peaklist, header_frag, tolerance=0.5):

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


def group_scans(filename, nh=2, min_replicates=1, report=None, max_injection_time=None, merge_ms1=False, split=False):

    if filename.lower().endswith(".mzml"):
        d = mzml_portal.Mzml(filename)
    elif filename.lower().endswith(".raw"):
        d = thermo_raw_portal.ThermoRaw(filename)
    else:
        raise IOError("Incorrect file format: {}".format(os.path.basename(filename)))

    si = d.scan_ids()
    sd = d.scan_dependents()
    sit = d.injection_times()

    graphs = create_graphs_from_scan_ids(sd, si, sit)
    if max_injection_time:
        graphs = [G for G in graphs if filter_ms1_by_injection_time(G, max_injection_time)]

    if split is False:
        templates = create_templates(graphs, nh)
        graphs_grouped = match_templates(graphs, templates, nh)
        graphs_grouped_v2 = match_templates_v2(graphs, templates)
        print graphs_grouped_v2 == graphs_grouped
        #raw_input()
    else:
        graphs_grouped = graphs

    if merge_ms1:
        graphs_grouped = merge_ms1_scans(graphs_grouped)

    graphs_grouped = [filter_by_replicates(G, min_replicates) for G in graphs_grouped]

    for i, G in enumerate(graphs_grouped):
        G.graph['id'] = i + 1

    if report is not None:
        out = open(report, "w")
        out.write("tree_id\tflag\ttemplate\tn\tincluded\n")

    valid_graphs = []
    for i, G in enumerate(graphs_grouped):

        e = list(G.edges(data=True))[0]  # MS1 to MS2
        flag = int((G.node[e[0]]["scanids"]) > 0 and len(G.node[e[1]]["scanids"]) > 0)

        if report is not None:
            scan_events = [n[0] for n in G.nodes(data=True) if len(n[1]["scanids"]) > 0]
            if not scan_events:
                scan_events = None
                n_scans = 0
            else:
                n_scans = len(scan_events)

            out.write("{}\t{}\t{}\t{}\t{}\n".format(G.graph['id'], flag, list(G.nodes())[0:nh], n_scans, scan_events))

        for n in G.nodes(data=True):
            if len(n[1]["scanids"]) == 0:
                G.remove_node(n[0])

        if flag:
            valid_graphs.append(G)

    if len(valid_graphs) == 0:
        raise IOError("No scan events remaining after filtering. Remove file from dataset or alter parameters.")

    return valid_graphs


def process_scans(filename, groups, function_noise, snr_thres, ppm, min_fraction=None, rsd_thres=None, normalise=False,
                  ringing_thres=None, exclusion_list={}, skip_ms1=False, report=None, block_size=2000, ncpus=None):
    print
    print os.path.basename(filename)

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
        print "Processing scans...."
        print "\n".join(map(str, [n[0] for n in nodes]))
        print
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
                    logging.warning("No scan data available for {}".format(n[0]))
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
    return pls_avg


def create_spectral_trees(trees, peaklists):

    spectral_trees = []

    headers = [pl.ID.split("#")[1] for pl in peaklists]

    for i, G in enumerate(trees):
        GG = nx.OrderedDiGraph()
        GG.graph["id"] = G.graph["id"]
        for edge in list(G.edges(data=True)):

            header_prec = "{}:{}".format(G.graph["id"], edge[0])
            if len(G.node[edge[0]]["scanids"]) == 0 or header_prec not in headers:
                if " ms " in header_prec:
                    logging.warning("Cannot create a spectral tree without precursor from {}".format(header_prec))
                    break
                continue

            pl = peaklists[headers.index(header_prec)]
            mz_prec, intensity_prec = assign_precursor(pl, edge[1], tolerance=0.5)

            if not mz_prec:
                if " ms " in header_prec:
                    logging.warning("Cannot create a spectral tree without precursor from {}".format(header_prec))
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

