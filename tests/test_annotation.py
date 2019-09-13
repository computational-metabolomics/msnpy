#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import copy
import zipfile

from msnpy.processing import *
from msnpy.annotation import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)


def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class AnnotationTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        zip_ref = zipfile.ZipFile(to_test_data("mzml_DIMSn.zip"), 'r')
        zip_ref.extractall(to_test_results(""))
        zip_ref.close()

        filename_01 = "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML"
        groups = group_scans(to_test_results(filename_01), 2, min_replicates=[1, 3, 3],
                             max_injection_time=0.0, merge_ms1=False)
        pls = process_scans(to_test_results(filename_01), groups=groups, function_noise="median", snr_thres=3.0,
                            ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                            report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)

        cls.trees = create_spectral_trees(groups, pls)

    def test_ApiMfdb(self):

        mfdb = ApiMfdb()

        # mz_d = mz_pair_diff_tol(lower_mz=100.000000, upper_mz=118.010565, tol=5, unit="ppm")
        # records = mfdb.select_mf(mz_d[0], mz_d[1], adducts=None, rules=False)
        # for record in records:
        #     print(record)

        mz_m = mz_tol(mz=71.03711 + (1.007825 - 0.0005486), tol=1.0, unit="ppm")
        records = mfdb.select_mf(mz_m[0], mz_m[1], adducts=["[M+H]+"], rules=True)
        self.assertDictEqual(records[0], {'mass': 72.0443904, 'atoms':
                                          {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'P': 0, 'S': 0},
                                          'adduct': '[M+H]+', 'DBE': 2, 'LEWIS': 1, 'SENIOR': 1, 'HC': 1, 'NOPSC': 1})

    def test_annotate_mf(self):
        pass
        #annotate_mf(spectral_trees: Sequence[nx.OrderedDiGraph], db_out: str, ppm: float,
        #                adducts: list = ["[M+H]+"],
        #                rules: bool = True, mf_db: str = "http://mfdb.bham.ac.uk", prefix_inp: str = ""):

    def test_mf_tree(self):
        pass
        # mf_tree(G: nx.OrderedDiGraph, path_db: str, max_mslevel: int, prefix: str)

    def test_filter_mf(self):
        pass
        # filter_mf(trees: Sequence[nx.OrderedDiGraph], path_db: str):

    def test_print_formula(self):
        self.assertEqual(print_formula({"C":6, "H":12, "N":0, "O":6, "P":0, "S":0}), "C6H12O6")
        self.assertEqual(print_formula({"C": 0, "H": 2, "N": 0, "O": 1, "P": 0, "S": 0}), "H2O")

    def test_rank_mf(self):
        pass
        # rank_mf(trees: Sequence[nx.OrderedDiGraph]):