#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import copy
import zipfile

from msnpy.processing import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class ProcessingTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        zip_ref = zipfile.ZipFile(to_test_data("mzml_DIMSn.zip"), 'r')
        zip_ref.extractall(to_test_results(""))
        zip_ref.close()

        cls.filename_01 = "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML"

        cls.groups = group_scans(to_test_results(cls.filename_01), 2, min_replicates=[1, 3, 3],
                                 max_injection_time=0.0, merge_ms1=False)

        cls.pls = process_scans(to_test_results(cls.filename_01), groups=cls.groups, function_noise="median",
                                snr_thres=3.0, ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                                report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)

        cls.trees = create_spectral_trees(cls.groups, cls.pls)

    def test_grouping(self):

        def _assert(node, event, scan_ids, ioninjectiontimes, mslevel, coltype, colenergy, template):
            self.assertEqual(node[0], event)
            self.assertListEqual(node[1]["scanids"], scan_ids)
            self.assertListEqual(node[1]["ioninjectiontimes"], ioninjectiontimes)
            self.assertEqual(node[1]["mslevel"], mslevel)
            self.assertEqual(node[1]["coltype"], coltype)
            self.assertEqual(node[1]["colenergy"], colenergy)
            self.assertEqual(node[1]["template"], template)


        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_01.txt"), max_injection_time=0.0,
                             merge_ms1=False)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [190.00-1200.00]', [1], [1.902963042259],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]', [3, 6, 9], [1000.0, 1000.0, 1000.0],
                2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [190.00-1200.00]', [502], [6.463267803192],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@cid35.00 [90.00-390.00]', [504, 507, 511],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_02.txt"), max_injection_time=0.0,
                             merge_ms1=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [190.00-1200.00]', [1, 12, 497, 502],
                [1.902963042259, 7.802415370941, 1000.0, 6.463267803192], 1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]', [3, 6, 9], [1000.0, 1000.0, 1000.0],
                2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [190.00-1200.00]', [1, 12, 497, 502],
                [1.902963042259, 7.802415370941, 1000.0, 6.463267803192], 1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@cid35.00 [90.00-390.00]', [504, 507, 511],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_03.txt"), max_injection_time=0.0,
                             merge_ms1=False, split=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [190.00-1200.00]', [1], [1.902963042259], 1, None, 0.0, False)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_03.txt"), max_injection_time=0.5,
                             merge_ms1=False)

        self.assertEqual(len(groups), 0)

    def test_processing(self):
        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             max_injection_time=0.0, merge_ms1=False)

        pls = process_scans(to_test_results(self.filename_01), groups=self.groups, function_noise="median", snr_thres=3.0,
                            ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                            report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)

        self.assertEqual(pls[0].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(pls[0].mz[0], 195.03663635253906)
        self.assertEqual(pls[0].mz[-1], 1197.8302001953125)
        self.assertEqual(list(pls[0].__dict__.keys()),
                         ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls[0].intensity[0], 103909.3046875)
        self.assertEqual(pls[0].intensity[-1], 93408.7578125)
        self.assertEqual(len(pls[0]), 223)
        self.assertEqual(pls[0].to_str()[0:52], "mz,intensity,intensity_norm,snr,present,fraction,rsd")
        self.assertEqual(len(pls[0].to_str()), 21723)

        self.assertEqual(pls[1].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI d Full ms2 376.34@hcd40.00 [50.00-390.00]")
        self.assertEqual(len(pls[1].mz), 0)

        self.assertEqual(pls[2].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI d Full ms2 376.34@hcd80.00 [50.00-390.00]")
        self.assertEqual(len(pls[2].mz), 0)

        self.assertEqual(pls[3].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#2:FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(len(pls[3].mz), 344)

        self.assertEqual(pls[4].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#2:FTMS + p NSI d Full ms2 554.55@hcd20.00 [50.00-565.00]")
        self.assertEqual(pls[4].mz[0], 206.2425791422526)
        self.assertEqual(pls[4].mz[-1], 554.5512084960938)
        self.assertEqual(pls[4].intensity[0], 80.20240529378255)
        self.assertEqual(pls[4].intensity[-1], 110562.52864583333)
        self.assertEqual(pls[4].intensity_norm[0], 0.0006135623314597699)
        self.assertEqual(pls[4].intensity_norm[-1], 0.83922192977622943)
        self.assertTupleEqual(pls[4].attributes,
                              ('mz', 'intensity', 'intensity_norm', 'snr', 'present', 'fraction', 'rsd_intensity_norm',
                               'occurrence', 'purity', 'snr_flag', 'fraction_flag', 'rsd_intensity_norm_flag'))
        self.assertEqual(pls[4].rsd_intensity_norm[0], 23.725594401714105)
        self.assertEqual(pls[4].rsd_intensity_norm[-1], 2.0593022319113317)
        self.assertEqual(len(pls[4]), 9)
        self.assertEqual(pls[4].to_str()[0:52], "mz,intensity,intensity_norm,snr,present,fraction,rsd")
        self.assertEqual(len(pls[4].to_str()), 2621)

    def test_create_spectral_trees(self):

        nodes = list(self.trees[0].nodes(data=True))
        edges = list(self.trees[0].edges(data=True))

        self.assertEqual(nodes, [])
        self.assertEqual(edges, [])

        nodes = list(self.trees[3].nodes(data=True))
        edges = list(self.trees[3].edges(data=True))

        self.assertEqual(nodes[0][0], "376.3435_9_177")
        self.assertListEqual(list(nodes[0][1].keys()), ['mz', 'intensity', 'header', 'mslevel', 'precursor'])
        self.assertEqual(nodes[0][1]["mz"], 376.3435363769531)
        self.assertEqual(nodes[0][1]["intensity"], 518535.78125)
        self.assertEqual(nodes[0][1]["header"], "FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(nodes[0][1]["mslevel"], 1)
        self.assertEqual(nodes[0][1]["precursor"], True)

        self.assertEqual(nodes[-1][0], "311.2957_11_2")
        self.assertEqual(nodes[-1][1]["mz"], 311.295654296875)
        self.assertEqual(nodes[-1][1]["intensity"], 253.14910888671875)
        self.assertEqual(nodes[-1][1]["header"], "FTMS + p NSI d Full ms3 376.34@cid35.00 326.27@cid35.00 [75.00-340.00]")
        self.assertEqual(nodes[-1][1]["mslevel"], 3)
        self.assertEqual(nodes[-1][1]["precursor"], False)

        self.assertEqual(edges[0][0], "376.3435_9_177")
        self.assertEqual(edges[0][1], "209.1182_10_0")
        self.assertListEqual(list(edges[0][2].keys()), ['mzdiff', 'type'])
        self.assertEqual(edges[0][2]["mzdiff"], 167.2253036)
        self.assertEqual(edges[0][2]["type"], "e")

        self.assertEqual(edges[-1][0], "326.2701_10_9")
        self.assertEqual(edges[-1][1], "311.2957_11_2")
        self.assertEqual(edges[-1][2]["mzdiff"], 14.9744771)
        self.assertEqual(edges[-1][2]["type"], "e")

    def test_create_graphs_from_scan_ids(self):
        sd = [[1, 2],[2, 3]]
        se = {1: "h1", 2: "h2", 3: "h3"}
        iit = {1: 0.5, 2: 0.5, 3: 0.5}
        graphs = create_graphs_from_scan_ids(scan_dependents=sd, scan_events=se, ion_injection_times=iit)
        nodes = list(graphs[0].nodes(data=True))
        edges = list(graphs[0].edges(data=True))
        self.assertEqual(nodes[0][0], 'h1')
        self.assertDictEqual(nodes[0][1], {'scanids': [1], 'mslevel': 1, 'coltype': None, 'colenergy': 0.0,
                                        'ioninjectiontimes': [0.5], 'flag': True})
        self.assertEqual(nodes[1][0], 'h2')
        self.assertDictEqual(nodes[1][1], {'scanids': [2], 'mslevel': 1, 'coltype': None, 'colenergy': 0.0,
                                        'ioninjectiontimes': [0.5], 'flag': True})

        self.assertEqual(edges[0][0], 'h1')
        self.assertEqual(edges[0][1], 'h2')
        self.assertDictEqual(edges[0][2], {})

    def test_assign_precursor(self):
        ap = assign_precursor(peaklist=self.pls[0], header_frag=self.pls[1].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (376.3430480957031, 2519428.5))

        ap = assign_precursor(peaklist=self.pls[0], header_frag=self.pls[2].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (376.3430480957031, 2519428.5))

        ap = assign_precursor(peaklist=self.pls[3], header_frag=self.pls[4].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (554.5517578125, 143097.390625))

        ap = assign_precursor(peaklist=self.pls[3], header_frag=self.pls[5].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (554.5517578125, 143097.390625))

        ap = assign_precursor(peaklist=self.pls[0], header_frag=self.pls[0].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (None, None))

    def test_mz_tol(self):
        mzt = mz_tol(mz=100.0, tol=5.0, unit="ppm")
        self.assertEqual(mzt, (99.9995, 100.0005))

        mpd = mz_pair_diff_tol(lower_mz=100.0, upper_mz=150.0, tol=5, unit="ppm")
        self.assertEqual(mpd, (49.99874999999999, 50.00125000000001))

    def test_merge_ms1_scans(self):
        graphs = merge_ms1_scans(graphs=copy.deepcopy(self.groups))
        nodes = list(graphs[0].nodes(data=True))
        self.assertListEqual(nodes[0][1]["scanids"], [1, 12, 497, 502])

    def test_templates(self):

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3], max_injection_time=0.0,
                             merge_ms1=False, remove=False)
        templ = create_templates(groups, 2)
        nodes = list(templ[0].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(nodes[0][1]["scanids"], [])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [])
        self.assertEqual(nodes[1][0], "FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]")
        self.assertEqual(nodes[1][1]["scanids"], [])
        self.assertEqual(nodes[1][1]["ioninjectiontimes"], [])
        nodes = list(templ[1].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(nodes[0][1]["scanids"], [])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [])

        gbt = group_by_template(graphs=groups, templates=templ)
        nodes = list(gbt[0].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(nodes[0][1]["scanids"], [1])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [1.902963042259])
        self.assertEqual(nodes[1][0], "FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]")
        self.assertEqual(nodes[1][1]["scanids"], [3, 6, 9])
        self.assertEqual(nodes[1][1]["ioninjectiontimes"], [1000.0, 1000.0, 1000.0])
        nodes = list(gbt[1].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI Full ms [190.00-1200.00]")
        self.assertEqual(nodes[0][1]["scanids"], [12])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [7.802415370941])
