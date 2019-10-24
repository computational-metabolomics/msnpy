#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import os
import io
import unittest
import tempfile
import shutil
import numpy as np
from msnpy.convert import *
from dimspy.portals.hdf5_portal import load_peaklists_from_hdf5

def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class ConvertTestCase(unittest.TestCase):

    def test_tree2peaklist(self):

        non_merged_pls, merged_pls, ms1_precursor_pl = tree2peaklist(tree_pth=to_test_data('spectral_trees.json'),
                                                                     out_pth='',
                                                                     name='test',
                                                                     adjust_mz=False,
                                                                     merge=True,
                                                                     ppm=5)

        non_merged_original = load_peaklists_from_hdf5(
            to_test_data('convert', 'test_non_merged_pls.hdf5')
        )

        for i in range(0, len(non_merged_original)):
            self.assertTrue(np.array_equal(non_merged_pls[i].peaks,
                                           non_merged_original[i].peaks))

        merged_original = load_peaklists_from_hdf5(
            to_test_data('convert', 'test_merged_pls.hdf5')
        )
        for i in range(0, len(merged_original)):
            self.assertTrue(np.array_equal(merged_pls[i].peaks,
                                           merged_original[i].peaks))

        ms1_precursor_original = load_peaklists_from_hdf5(
            to_test_data('convert', 'test_ms1_precursors_pl.hdf5')
        )

        for i in range(0, len(ms1_precursor_original)):
            self.assertTrue(np.array_equal(ms1_precursor_pl[i].peaks,
                                           ms1_precursor_original[i].peaks))

    def test_peaklist2msp_non_merged(self):
        #test_out = to_test_data(os.path.join('convert',
        # 'test_non_merged.msp'))
        test_out = os.path.join(tempfile.mkdtemp(), 'test_non_merged.msp')
        pl = load_peaklists_from_hdf5(
            to_test_data('convert', 'test_non_merged_pls.hdf5')
        )
        peaklist2msp(pl,
                     test_out, #test_out
                     msp_type='massbank',
                     polarity='positive')
        self.assertListEqual(
            list(io.open(test_out)),
            list(io.open(to_test_data(os.path.join(
                'convert', 'test_non_merged.msp'))))
        )

    def test_peaklist2msp_merged(self):

        test_out = to_test_data('convert', 'test_merged.msp')
        #test_out = os.path.join(tempfile.mkdtemp(), 'test_merged.msp')
        peaklist2msp(load_peaklists_from_hdf5(to_test_data('convert', 'test_merged_pls.hdf5')),
                     test_out,
                     msp_type='massbank',
                     polarity='positive')
        self.assertListEqual(
            list(io.open(test_out)),
            list(io.open(to_test_data('convert', 'test_merged.msp'))))

    def test_peaklist2msp_ms1_precursors(self):
        test_out = to_test_data('convert', 'test_ms1_precursors.msp')
        #test_out = os.path.join(tempfile.mkdtemp(), 'test_ms1_precursors.msp')
        peaklist2msp(load_peaklists_from_hdf5(to_test_data(
                                    'convert', 'test_ms1_precursors_pl.hdf5')),
                     test_out,
                     msp_type='massbank',
                     polarity='positive',
                     include_ms1=True)
        self.assertListEqual(
            list(io.open(test_out)),
            list(io.open(to_test_data('convert','test_ms1_precursors.msp'))))