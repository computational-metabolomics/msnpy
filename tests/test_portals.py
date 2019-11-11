#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import zipfile

from msnpy.processing import *
from msnpy.portals import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)


def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class PortalsTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        filename_01 = to_test_data("A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML")
        cls.groups = group_scans(to_test_results(filename_01), 2, min_replicates=[1, 3, 3],
                             max_injection_time=0.0, merge_ms1=False)
        pls = process_scans(to_test_results(filename_01), groups=cls.groups, function_noise="median", snr_thres=3.0,
                            ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                            report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)
        cls.trees = create_spectral_trees(cls.groups, pls)

    def _compare_trees(self, G, GG):

        for i in range(len(G)):
            nodes = list(G[i].nodes(data=True))
            for j, node in enumerate((GG[i].nodes(data=True))):
                self.assertEqual(node[0], nodes[j][0])
                self.assertDictEqual(node[1], nodes[j][1])
            edges = list(G[i].edges(data=True))
            for j, edge in enumerate((GG[i].edges(data=True))):
                self.assertEqual(edge[0], edges[j][0])
                self.assertEqual(edge[1], edges[j][1])
                self.assertDictEqual(edge[2], edges[j][2])
            self.assertEqual(G[i].graph["id"], GG[i].graph["id"])

    def test_load_groups(self):

        save_groups(self.groups, to_test_results("groups.json"), format="json")
        groups_result = load_groups(to_test_results("groups.json"), format="json")
        self._compare_trees(self.groups, groups_result)

        save_groups(self.groups, to_test_results("groups.gml"), format="gml")
        groups_result = load_groups(to_test_results("groups.gml"), format="gml")
        self._compare_trees(self.groups, groups_result)

    def test_load_trees(self):

        save_trees(self.trees, to_test_results("spectral_trees.json"), format="json")
        trees_result = load_trees(to_test_results("spectral_trees.json"), format="json")
        self._compare_trees(self.trees, trees_result)

        save_trees(self.trees, to_test_results("spectral_trees.gml"), format="gml")
        trees_result = load_trees(to_test_results("spectral_trees.gml"), format="gml")
        self._compare_trees(self.trees, trees_result)
