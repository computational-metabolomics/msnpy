#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Ralf Weber
#
# This file is part of MSnPy.
#
# MSnPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MSnPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MSnPy.  If not, see <https://www.gnu.org/licenses/>.
#


import copy
import unittest

from msnpy.processing import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)


def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class ProcessingTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # zip_ref = zipfile.ZipFile(to_test_data("mzml_DIMSn.zip"), 'r')
        # zip_ref.extractall(to_test_results(""))
        # zip_ref.close()

        cls.filename_01 = to_test_data("A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML")
        cls.filename_02 = to_test_data("A06_Polar_Daph_WAX1_Phenyl_LCMS_Pos_DIMSn_subset.mzML")

        cls.groups_01 = group_scans(to_test_results(cls.filename_01), 2, min_replicates=[1, 3, 3],
                                    max_injection_time=0.0, merge_ms1=False)
        cls.groups_02 = group_scans(to_test_results(cls.filename_02), 2, min_replicates=[1, 3, 3],
                                    max_injection_time=0.0, merge_ms1=False)

        cls.pls_01 = process_scans(to_test_results(cls.filename_01), groups=cls.groups_01, function_noise="median",
                                   snr_thres=3.0, ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                                   report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)
        cls.pls_02 = process_scans(to_test_results(cls.filename_02), groups=cls.groups_02, function_noise="median",
                                   snr_thres=3.0, ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                                   report=to_test_results("processing_02.txt"), block_size=5000, ncpus=None)

        cls.trees_01 = create_spectral_trees(cls.groups_01, cls.pls_01)
        cls.trees_02 = create_spectral_trees(cls.groups_02, cls.pls_02)


    def test_grouping(self):

        def _assert(node, event, scan_ids, ioninjectiontimes, mslevel, coltype, colenergy, template):
            self.assertEqual(node[0], event)
            self.assertListEqual(node[1]["scanids"], scan_ids)
            self.assertListEqual(node[1]["ioninjectiontimes"], ioninjectiontimes)
            self.assertEqual(node[1]["mslevel"], mslevel)
            self.assertEqual(node[1]["coltype"], coltype)
            self.assertEqual(node[1]["colenergy"], colenergy)
            self.assertEqual(node[1]["template"], template)

        # -----------------------------------------------------
        # A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML
        # -----------------------------------------------------
        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_01_01.txt"), max_injection_time=0.0,
                             merge_ms1=False)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI d SIM ms [372.00-382.00]', [2], [369.747467041016],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]', [3, 6, 9], [1000.0, 1000.0, 1000.0],
                2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI d SIM ms [372.00-382.00]', [503], [573.602966308594],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@cid35.00 [90.00-390.00]', [504, 507, 511],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_01_02.txt"), max_injection_time=0.0,
                             merge_ms1=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI d SIM ms [372.00-382.00]', [2, 503],
                [369.747467041016, 573.602966308594], 1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]', [3, 6, 9], [1000.0, 1000.0, 1000.0],
                2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI d SIM ms [372.00-382.00]', [2, 503],
                [369.747467041016, 573.602966308594], 1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 376.34@cid35.00 [90.00-390.00]', [504, 507, 511],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_01_03.txt"), max_injection_time=0.0,
                             merge_ms1=False, split=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI d SIM ms [372.00-382.00]', [2], [369.747467041016], 1, None, 0.0, False)

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_01_04.txt"), max_injection_time=0.5,
                             merge_ms1=False)

        self.assertEqual(len(groups), 0)

        # -----------------------------------------------------
        # A06_Polar_Daph_WAX1_Phenyl_LCMS_Pos_DIMSn_subset.mzML
        # -----------------------------------------------------
        groups = group_scans(to_test_results(self.filename_02), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_02_01.txt"), max_injection_time=0.0,
                             merge_ms1=False)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [50.00-1000.00]', [1], [17.388236999512],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 137.05@hcd20.00 [50.00-150.00]', [2, 5, 8],
                [9.576781272888, 9.576781272888, 9.576781272888], 2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [50.00-1000.00]', [804], [19.371095657349],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 281.16@cid35.00 [65.00-295.00]', [805, 809, 813],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_02), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_02_02.txt"), max_injection_time=0.0,
                             merge_ms1=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [50.00-1000.00]', [1, 11, 21, 791, 804],
                [17.388236999512, 21.448980331421, 21.765810012817, 18.100595474243, 19.371095657349],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 137.05@hcd20.00 [50.00-150.00]', [2, 5, 8],
                [9.576781272888, 9.576781272888, 9.576781272888], 2, "hcd", 20.0, True)

        nodes = list(groups[-1].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [50.00-1000.00]', [1, 11, 21, 791, 804],
                [17.388236999512, 21.448980331421, 21.765810012817, 18.100595474243, 19.371095657349],
                1, None, 0.0, True)
        _assert(nodes[1], 'FTMS + p NSI d Full ms2 281.16@cid35.00 [65.00-295.00]', [805, 809, 813],
                [1000.0, 1000.0, 1000.0],
                2, "cid", 35.0, True)

        groups = group_scans(to_test_results(self.filename_02), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_02_03.txt"), max_injection_time=0.0,
                             merge_ms1=False, split=True)

        nodes = list(groups[0].nodes(data=True))
        _assert(nodes[0], 'FTMS + p NSI Full ms [50.00-1000.00]', [1], [17.388236999512], 1, None, 0.0, False)

        groups = group_scans(to_test_results(self.filename_02), 2, min_replicates=[1, 3, 3],
                             report=to_test_results("report_grouping_02_04.txt"), max_injection_time=0.5,
                             merge_ms1=False)

        self.assertEqual(len(groups), 0)


    def test_processing(self):
        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3],
                             max_injection_time=0.0, merge_ms1=False)

        pls = process_scans(to_test_results(self.filename_01), groups=self.groups_01, function_noise="median", snr_thres=3.0,
                            ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                            report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)

        self.assertEqual(pls[0].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI d SIM ms [372.00-382.00]")
        self.assertEqual(pls[0].mz[0], 372.1014404296875)
        self.assertEqual(pls[0].mz[-1], 379.3191223144531)
        self.assertEqual(list(pls[0].__dict__.keys()),
                         ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls[0].intensity[0], 6439.35986328125)
        self.assertEqual(pls[0].intensity[-1], 12307.2763671875)
        self.assertEqual(len(pls[0]), 12)
        self.assertEqual(pls[0].to_str()[0:52], "mz,intensity,intensity_norm,snr,present,fraction,rsd")
        self.assertEqual(len(pls[0].to_str()), 1312)

        self.assertEqual(pls[1].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI d Full ms2 376.34@hcd40.00 [50.00-390.00]")
        self.assertEqual(len(pls[1].mz), 0)

        self.assertEqual(pls[2].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#1:FTMS + p NSI d Full ms2 376.34@hcd80.00 [50.00-390.00]")
        self.assertEqual(len(pls[2].mz), 0)

        self.assertEqual(pls[3].ID,
                         "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML#2:FTMS + p NSI d SIM ms [550.00-560.00]")
        self.assertEqual(len(pls[3].mz), 54)

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

        nodes = list(self.trees_01[0].nodes(data=True))
        edges = list(self.trees_01[0].edges(data=True))

        self.assertEqual(nodes, [])
        self.assertEqual(edges, [])

        nodes = list(self.trees_01[3].nodes(data=True))
        edges = list(self.trees_01[3].edges(data=True))

        self.assertEqual(nodes[0][0], "538.1646_11_32")
        self.assertListEqual(list(nodes[0][1].keys()), ['mz', 'intensity', 'header', 'mslevel', 'precursor'])
        self.assertEqual(nodes[0][1]["mz"], 538.1646118164062)
        self.assertEqual(nodes[0][1]["intensity"], 66509.046875)
        self.assertEqual(nodes[0][1]["header"], "FTMS + p NSI d SIM ms [534.00-544.00]")
        self.assertEqual(nodes[0][1]["mslevel"], 1)
        self.assertEqual(nodes[0][1]["precursor"], True)

        self.assertEqual(nodes[-1][0], "538.5851_12_12")
        self.assertEqual(nodes[-1][1]["mz"], 538.5850626627604)
        self.assertEqual(nodes[-1][1]["intensity"], 79.27662404378255)
        self.assertEqual(nodes[-1][1]["header"], "FTMS + p NSI d Full ms2 538.17@cid35.00 [135.00-550.00]")
        self.assertEqual(nodes[-1][1]["mslevel"], 2)
        self.assertEqual(nodes[-1][1]["precursor"], False)

        self.assertEqual(edges[0][0], "538.1646_11_32")
        self.assertEqual(edges[0][1], "254.2855_12_0")
        self.assertListEqual(list(edges[0][2].keys()), ['mzdiff', 'type'])
        self.assertEqual(edges[0][2]["mzdiff"], 283.8791606)
        self.assertEqual(edges[0][2]["type"], "e")

        self.assertEqual(edges[-1][0], "538.1646_11_32")
        self.assertEqual(edges[-1][1], "538.5851_12_12")
        self.assertEqual(edges[-1][2]["mzdiff"], -0.4204508)
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

        ap = assign_precursor(peaklist=self.pls_01[0], header_frag=self.pls_01[1].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (376.3426208496094, 538519.8125))

        ap = assign_precursor(peaklist=self.pls_02[0], header_frag=self.pls_02[1].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (137.04583740234375, 11351725.0))

        ap = assign_precursor(peaklist=self.pls_02[0], header_frag=self.pls_02[2].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (137.04583740234375, 11351725.0))

        ap = assign_precursor(peaklist=self.pls_02[0], header_frag=self.pls_02[0].ID.split(":")[-1], tolerance=0.5)
        self.assertTupleEqual(ap, (None, None))

    def test_mz_tol(self):
        mzt = mz_tol(mz=100.0, tol=5.0, unit="ppm")
        self.assertEqual(mzt, (99.9995, 100.0005))

        mpd = mz_pair_diff_tol(lower_mz=100.0, upper_mz=150.0, tol=5, unit="ppm")
        self.assertEqual(mpd, (49.99874999999999, 50.00125000000001))

    def test_merge_ms1_scans(self):
        graphs = merge_ms1_scans(graphs=copy.deepcopy(self.groups_01))
        nodes = list(graphs[0].nodes(data=True))
        self.assertListEqual(nodes[0][1]["scanids"], [2, 503])

        graphs = merge_ms1_scans(graphs=copy.deepcopy(self.groups_02))
        nodes = list(graphs[0].nodes(data=True))
        self.assertListEqual(nodes[0][1]["scanids"], [1, 11, 21, 791, 804])

    def test_templates(self):

        groups = group_scans(to_test_results(self.filename_01), 2, min_replicates=[1, 3, 3], max_injection_time=0.0,
                             merge_ms1=False, remove=False)
        templ = create_templates(groups, 2)
        nodes = list(templ[0].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI d SIM ms [372.00-382.00]")
        self.assertEqual(nodes[0][1]["scanids"], [])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [])
        self.assertEqual(nodes[1][0], "FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]")
        self.assertEqual(nodes[1][1]["scanids"], [])
        self.assertEqual(nodes[1][1]["ioninjectiontimes"], [])
        nodes = list(templ[1].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI d SIM ms [550.00-560.00]")
        self.assertEqual(nodes[0][1]["scanids"], [])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [])

        gbt = group_by_template(graphs=groups, templates=templ)
        nodes = list(gbt[0].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI d SIM ms [372.00-382.00]")
        self.assertEqual(nodes[0][1]["scanids"], [2])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [369.747467041016])
        self.assertEqual(nodes[1][0], "FTMS + p NSI d Full ms2 376.34@hcd20.00 [50.00-390.00]")
        self.assertEqual(nodes[1][1]["scanids"], [3, 6, 9])
        self.assertEqual(nodes[1][1]["ioninjectiontimes"], [1000.0, 1000.0, 1000.0])
        nodes = list(gbt[1].nodes(data=True))
        self.assertEqual(nodes[0][0], "FTMS + p NSI d SIM ms [550.00-560.00]")
        self.assertEqual(nodes[0][1]["scanids"], [13])
        self.assertEqual(nodes[0][1]["ioninjectiontimes"], [882.302795410156])
