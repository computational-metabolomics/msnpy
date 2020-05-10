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
from random import shuffle
from shutil import copyfile
from sqlite3 import connect

from pandas.testing import assert_frame_equal

from msnpy.annotation import *
from msnpy.processing import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)


def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class AnnotationTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # zip_ref = zipfile.ZipFile(to_test_data("mzml_DIMSn.zip"), 'r')
        # zip_ref.extractall(to_test_results(""))
        # zip_ref.close()

        filename_01 = to_test_data("pos_21-hydroxyprogesterone_subset.mzML")
        groups = group_scans(to_test_results(filename_01), 2, min_replicates=1,
                             max_injection_time=0.0, merge_ms1=False)
        pls = process_scans(to_test_results(filename_01), groups=groups, function_noise="median", snr_thres=3.0,
                            ppm=5.0, min_fraction=0.5, rsd_thres=30.0, normalise=True,
                            report=to_test_results("processing_01.txt"), block_size=5000, ncpus=None)

        cls.trees = create_spectral_trees(groups, pls)

        adducts = {"[M+H]+":  1.0072764}
        cls.db_single_adduct = to_test_results("test_mf_single_adduct.sqlite")
        cls.amf_single_adduct = annotate_mf(spectral_trees=cls.trees, db_out=cls.db_single_adduct, ppm=10.0, adducts=adducts, rules=True, mf_db="http://mfdb.bham.ac.uk")

        cls.db_single_adduct_filtered = cls.db_single_adduct.replace(".sqlite", "_filtered.sqlite")
        copyfile(cls.db_single_adduct, cls.db_single_adduct_filtered)

        cls.fmf = filter_mf(cls.trees, cls.db_single_adduct_filtered)

    def test_ApiMfdb(self):

        mfdb = ApiMfdb()

        # mz_d = mz_pair_diff_tol(lower_mz=100.000000, upper_mz=118.010565, tol=5, unit="ppm")
        # records = mfdb.select_mf(mz_d[0], mz_d[1], adducts=None, rules=False)
        # for record in records:
        #     print(record)

        mz_m = mz_tol(mz=71.03711 + (1.007825 - 0.0005486), tol=1.0, unit="ppm")
        records = mfdb.select_mf(mz_m[0], mz_m[1], adducts={"[M+H]+": 1.0072764}, rules=True)
        self.assertDictEqual(records[0], {'mass': 72.0443904, 'atoms':
                                          {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'P': 0, 'S': 0},
                                          'adduct': '[M+H]+', 'DBE': 2, 'LEWIS': 1, 'SENIOR': 1, 'HC': 1, 'NOPSC': 1})

    def test_annotate_mf(self):

        self.assertDictEqual(self.amf_single_adduct[0].nodes["331.2271_0_7"], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True})
        conn = connect(self.db_single_adduct)
        cursor = conn.cursor()
        cursor.execute("select * from MF_1")
        records = cursor.fetchall()

        self.assertTupleEqual(records[0], ('331.2271_0_7', 1, 14, 30, 6, 1, 0, 1, 3, 1, 1, 1, 1, '[M+H]+', 331.22745740000005, 331.22706604003906, -1.18, 1, '1', 0))
        self.assertTupleEqual(records[5], ('331.2271_0_7', 6, 21, 30, 0, 3, 0, 0, 7, 1, 1, 1, 1, '[M+H]+', 331.2267714, 331.22706604003906, 0.89, 1, '1', 0))
        self.assertTupleEqual(records[-1], ('313.2163_1_31__295.2056_3_0', 583, 0, 2, 0, 1, 0, 0, 0, 1, 1, 0, 1, None, 18.010565, 18.0106583, -5.18, None, '2_3', 0))

        cursor.execute("select count(*) from MF_1")
        self.assertEqual(cursor.fetchone()[0], 583)

        conn.close()

        e = 0.0005486
        adducts = {"[M+H]+": 1.007825 - e,
                   "[M+NH4]+": 18.034374 - e,
                   "[M+Na]+": 22.989770 - e
                   }
        db_out = to_test_results("test_mf_multiple_adducts.sqlite")

        amf = annotate_mf(spectral_trees=self.trees, db_out=db_out, ppm=10.0, adducts=adducts, rules=True, mf_db="http://mfdb.bham.ac.uk")
        self.assertDictEqual(amf[0].nodes["331.2271_0_7"], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True})

        conn = connect(db_out)
        cursor = conn.cursor()
        cursor.execute("select distinct ADDUCT from MF_1")
        records = cursor.fetchall()

        self.assertEqual(records[0][0], "[M+H]+")
        self.assertEqual(records[1][0], "[M+NH4]+")
        self.assertEqual(records[2][0], "[M+Na]+")

        cursor.execute("select count(*) from MF_1")
        self.assertEqual(cursor.fetchone()[0], 697)

        conn.close()

    def test_filter_mf(self):

        self.assertEqual(len(self.fmf), 6)

        self.assertDictEqual(self.fmf[0].nodes['331.2271_0_7'], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True, 'mf':
                                                           {'1': {'mass': 331.22745740000005, 'mf': 'C14H30N6OS', 'adduct': '[M+H]+'}}})

        self.assertDictEqual(self.fmf[5].nodes['331.2271_0_7'], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True, 'mf':
                                                           {'6': {'mass': 331.2267714, 'mf': 'C21H30O3', 'adduct': '[M+H]+'}}})

        self.assertEqual(self.fmf[0].number_of_nodes(), 9)
        self.assertEqual(self.fmf[0].number_of_edges(), 8)

        self.assertEqual(self.fmf[5].number_of_nodes(), 37)
        self.assertEqual(self.fmf[5].number_of_edges(), 36)

        conn = connect(self.db_single_adduct_filtered)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        records = cursor.fetchall()
        self.assertListEqual(records, [('MF_1',), ('EDGES_1',), ('MZ_PREC_FRAG_1',)])

    def test_mf_tree(self):
        mft = mf_tree(self.trees[0], self.db_single_adduct_filtered, max_mslevel=5, prefix="_1")

        self.assertEqual(len(mft), 6)

        self.assertDictEqual(mft[0].nodes['331.2271_0_7'], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True, 'mf':
                                                           {'1': {'mass': 331.22745740000005, 'mf': 'C14H30N6OS', 'adduct': '[M+H]+'}}})

        self.assertDictEqual(mft[5].nodes['331.2271_0_7'], {'mz': 331.22706604003906, 'intensity': 205494704.0, 'header': 'FTMS + c NSI Full ms [327.22-333.22]', 'mslevel': 1, 'precursor': True, 'mf':
                                                           {'6': {'mass': 331.2267714, 'mf': 'C21H30O3', 'adduct': '[M+H]+'}}})

        self.assertEqual(mft[0].number_of_nodes(), 9)
        self.assertEqual(mft[0].number_of_edges(), 8)

        self.assertEqual(mft[5].number_of_nodes(), 37)
        self.assertEqual(mft[5].number_of_edges(), 36)

    def test_print_formula(self):
        self.assertEqual(print_formula({"C": 6, "H": 12, "N": 0, "O": 6, "P": 0, "S": 0}), "C6H12O6")
        self.assertEqual(print_formula({"C": 0, "H": 2, "N": 0, "O": 1, "P": 0, "S": 0}), "H2O")

    def test_rank_mf(self):
        rmf = rank_mf(self.fmf)

        self.assertListEqual([('1_6', 1), ('1_2', 2), ('1_4', 2), ('1_5', 3), ('1_3', 4), ('1_1', 5)],
                             [(G.graph["id"], G.graph["rank"]) for G in rmf[0]])

        # rmf[1].to_csv(to_test_data("ranks_21-hydroxyprogesterone.txt"), sep="\t", index=False)
        rmf[1].to_csv(to_test_results("ranks_21-hydroxyprogesterone.txt"), sep="\t", index=False)

        assert_frame_equal(pd.read_csv(to_test_data("ranks_21-hydroxyprogesterone.txt"), sep="\t"),
                           pd.read_csv(to_test_results("ranks_21-hydroxyprogesterone.txt"), sep="\t"))

        rmf_filt = rank_mf(self.fmf, rank_threshold=2)
        self.assertListEqual([('1_6', 1), ('1_2', 2), ('1_4', 2)],
                             [(G.graph["id"], G.graph["rank"]) for G in rmf_filt[0]])

        # multiple groups and annotations
        fmf_copy = copy.deepcopy(self.fmf)
        for G in fmf_copy:
            G.graph["id"] = G.graph["id"].replace("1_", "2_")  # rename trees

        fmf_2 = self.fmf + fmf_copy
        shuffle(fmf_2)

        rmf_2 = rank_mf(fmf_2)
        self.assertListEqual([('1_6', 1), ('1_2', 2), ('1_4', 2), ('1_5', 3), ('1_3', 4), ('1_1', 5),
                              ('2_6', 1), ('2_2', 2), ('2_4', 2), ('2_5', 3), ('2_3', 4), ('2_1', 5)],
                             [(G.graph["id"], G.graph["rank"]) for G in rmf_2[0]])

        rmf_2 = rank_mf(fmf_2, rank_threshold=3)
        self.assertListEqual([('1_6', 1), ('1_2', 2), ('1_4', 2), ('1_5', 3),
                              ('2_6', 1), ('2_2', 2), ('2_4', 2), ('2_5', 3)],
                             [(G.graph["id"], G.graph["rank"]) for G in rmf_2[0]])