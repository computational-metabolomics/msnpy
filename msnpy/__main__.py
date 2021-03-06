# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import argparse
from msnpy.processing import group_scans
from msnpy.processing import process_scans
from msnpy.processing import create_spectral_trees
from msnpy.annotation import annotate_mf
from msnpy.annotation import filter_mf
from msnpy.annotation import rank_mf
from msnpy.portals import save_groups
from msnpy.portals import load_groups
from msnpy.portals import save_trees
from msnpy.portals import load_trees
from msnpy.convert import tree2peaklist
from msnpy.convert import peaklist2msp


from dimspy.portals import hdf5_portal
from msnpy import __version__


def map_delimiter(delimiter): # pragma: no cover
    seps = {"comma": ",", "tab": "\t"}
    if delimiter in seps:
        return seps[delimiter]
    else:
        return delimiter


def main(): # pragma: no cover

    print("Executing msnpy version %s." % __version__)

    parser = argparse.ArgumentParser(description='Python package to process and annotate MSn fragmentation data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_g = subparsers.add_parser('group-scans', help='Group fragmentation events and/or experiments.')
    parser_ps = subparsers.add_parser('process-scans', help='Process and filter scans.')

    parser_cst = subparsers.add_parser('create-spectral-trees', help='Create spectral trees from processed scan (fragmentation) data.')
    parser_ast = subparsers.add_parser('annotate-spectral-trees', help='Annotate and/or filter spectral trees.')
    parser_rst = subparsers.add_parser('rank-spectral-trees', help='Rank annotated spectral trees.')
    parser_cvst = subparsers.add_parser('convert-spectral-trees', help='Convert spectral trees to either dimspy peaklists, MSP files or both')

    ##################################################################
    # GROUP SCANS
    ##################################################################

    parser_g.add_argument('-i', '--input',
                          type=str, required=True,
                          help="Mzml or Thermo Scientific raw file")


    parser_g.add_argument('-o', '--output',
                          type=str, required=True,
                          help="")

    parser_g.add_argument('-u', '--report',
                          type=str, required=False, default=None,
                          help="Summary/Report of groups")

    parser_g.add_argument('-n', '--number-of-headers',
                          default=2, type=int, required=False,
                          help="")

    parser_g.add_argument('-r', '--min-replicates',
                          default=None, type=int, required=False,
                          help="")

    parser_g.add_argument('-t', '--max-injection-time',
                          default=None, type=float, required=False,
                          help="")

    parser_g.add_argument('-s', '--split',
                          action='store_true', required=False,
                          help="")

    parser_g.add_argument('-m', '--merge-ms1',
                          action='store_true', required=False,
                          help="")

    ##################################################################
    # PROCESS SCANS
    ##################################################################

    parser_ps.add_argument('-i', '--input',
                           type=str, required=True,
                           help="Mzml or Thermo Scientific raw file")

    parser_ps.add_argument('-g', '--groups',
                           type=str, required=True,
                           help="Json or gml file that includes the groups of scans")

    parser_ps.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects to.")

    parser_ps.add_argument('-m', '--function-noise',
                           choices=["median", "mean", "mad", "noise_packets"], required=True,
                           help="Select function to calculate noise.")

    parser_ps.add_argument('-s', '--snr-threshold',
                           default=3.0, type=float, required=True,
                           help="Signal-to-noise threshold")

    parser_ps.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_ps.add_argument('-a', '--min-fraction',
                           default=0.5, type=float, required=False,
                           help="Minimum fraction a peak has to be present. Use 0.0 to not apply this filter.")

    parser_ps.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Maximum threshold - relative standard deviation (Calculated for peaks that have been measured across a minimum of two scans).")

    parser_ps.add_argument('-n', '--normalise',
                           default=None, type=float, required=False,
                           help="Normalise scans by Total Ion Current (TIC)")

    parser_ps.add_argument('-e', '--exclusion-list', nargs='+',
                           default=None, required=False,
                           help="List of mz values to exclude from processing (e.g. from electrical noise)")

    parser_ps.add_argument('-r', '--ringing-threshold',
                           default=None, type=float, required=False,
                           help="Remove ringing artifacts.")

    parser_ps.add_argument('-u', '--report',
                           type=str, required=False, default=None,
                           help="Summary/Report of processed mass spectra")

    parser_ps.add_argument('-b', '--block-size',
                           default=5000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_ps.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")


    ##################################################################
    # CREATE SPECTRAL TREES
    ##################################################################

    parser_cst.add_argument('-i', '--input',
                            type=str, required=True,
                            help="HDF5 file (Peaklist objects) from step 'process-scans'.")

    parser_cst.add_argument('-g', '--groups',
                            type=str, required=True,
                            help="")

    parser_cst.add_argument('-o', '--output',
                            type=str, required=True,
                            help="")


    ##################################################################
    # ANNOTATE SPECTRAL TREES
    ##################################################################
    parser_ast.add_argument('-i', '--input',
                            type=str, required=True,
                            help="Json file containing spectral trees")

    parser_ast.add_argument('-p', '--ppm',
                            default=2.0, type=float, required=False,
                            help="Mass tolerance in Parts per million.")

    parser_ast.add_argument('-r', '--rules',
                            action='store_true', required=False,
                            help="")

    parser_ast.add_argument('-m', '--mf-db',
                            type=str, required=False, default="http://multiomics-int.cs.bham.ac.uk",
                            help="Molecular formulae database")

    parser_ast.add_argument('-d', '--output-db',
                            type=str, required=True,
                            help="Sqlite database file to store information regarding the annotations.")

    parser_ast.add_argument('-o', '--output-trees',
                            type=str, required=True,
                            help="Json file containing the annotated spectral trees.")

    parser_ast.add_argument('-a', '--adducts', nargs='+', required=True,
                            help="Adducts e.g. [M+H]+ [M+NH4]+ [M+Na]+ [M+(39K)]+",
                            default=['[M+H]+', '[M+Na]+', '[M+NH4]+'])

    parser_ast.add_argument('-f', '--filter',
                            action='store_true', required=False,
                            help="Filter the spectral tree annotations")



    #################################
    # RANK SPECTRAL TREES
    #################################
    parser_rst.add_argument('-i', '--input',
                            type=str, required=True,
                            help="Json file containing annotated spectral trees")

    parser_rst.add_argument('-o', '--output',
                            type=str, required=True,
                            help="Summary of the rankings")



    ##################################################################
    # CONVERT SPECTRA TREES - TO DIMSPY.PEAKLISTS AND MSP FILES
    ##################################################################
    parser_cvst.add_argument('-i', '--input',
                            type=str, required=True,
                            help="Json file containing annotated spectral trees or dimspy peaklist hdf5 file")

    parser_cvst.add_argument('-o', '--output',
                            type=str, required=True,
                            help="Out folder containing spectra")

    parser_cvst.add_argument('-x', '--input_type',
                             default="json", type=str, required=False,
                             help="If input is either a dimspy peaklist or a msnpy json")

    parser_cvst.add_argument('-n', '--name', type=str, required=False,
                           help="Name to use for suffixing files")

    parser_cvst.add_argument('-a', '--adjust_mz',
                             action='store_true', required=False,
                             help="Filter the spectral tree annotations")

    parser_cvst.add_argument('-m', '--merge',
                             action='store_true', required=False,
                             help="Filter the spectral tree annotations")

    parser_cvst.add_argument('-p', '--ppm',
                            default=5.0, type=float, required=False,
                            help="Mass tolerance in Parts per million.")

    parser_cvst.add_argument('-s', '--msp',
                             action='store_true', required=False,
                             help="Filter the spectral tree annotations")

    parser_cvst.add_argument('-t', '--msp_type',
                             default="massbank", type=str, required=False,
                             help="If MSP file is to be created what type (massbank, msp)")

    parser_cvst.add_argument('-z', '--polarity',
                           type=str, required=False, default='NA',
                           help="Polarity to add to the MSP file (positive or negative)")

    parser_cvst.add_argument('-y', '--ms1',
                             action='store_true', required=False,
                             help="Output ms1 spectra (creates spectra for the precursors in the MS1 spectra")


    args = parser.parse_args()

    print(args)

    if args.step == "group-scans":
        groups = group_scans(filename=args.input,
                             nh=args.number_of_headers,
                             min_replicates=args.min_replicates,
                             report=args.report,
                             max_injection_time=args.max_injection_time,
                             merge_ms1=args.merge_ms1,
                             split=args.split)

        save_groups(groups=groups,
                    filename=args.output,
                    format="json")

    if args.step == "process-scans":
        peaklists = process_scans(filename=args.input,
                                  groups=load_groups(args.groups, format="json") if args.groups else None,
                                  function_noise=args.function_noise,
                                  snr_thres=args.snr_threshold,
                                  ppm=args.ppm,
                                  min_fraction=args.min_fraction,
                                  rsd_thres=args.rsd_threshold,
                                  normalise=args.normalise,
                                  ringing_thres=args.ringing_threshold,
                                  exclusion_list=args.exclusion_list,
                                  report=args.report,
                                  block_size=args.block_size,
                                  ncpus=args.ncpus)

        hdf5_portal.save_peaklists_as_hdf5(peaklists, args.output)

    if args.step == "create-spectral-trees":
        groups = load_groups(args.groups, format="json")
        pls = hdf5_portal.load_peaklists_from_hdf5(args.input)
        spectral_trees = create_spectral_trees(groups, pls)
        save_trees(spectral_trees, args.output, format="json")

    if args.step == "annotate-spectral-trees":
        spectral_trees = load_trees(args.input, format="json")

        adducts = [a.replace('__ob__', '[').replace('__cb__', ']') for a in args.adducts]

        st = annotate_mf(spectral_trees=spectral_trees,
                         db_out=args.output_db,
                         ppm=args.ppm,
                         adducts=adducts,
                         rules=args.rules,
                         mf_db=args.mf_db)

        if args.filter:
            st = filter_mf(st, args.output_db)
        save_trees(st, args.output_trees, format="json")

    if args.step == "rank-spectral-trees":
        st = load_trees(args.input, format="json")
        ranks = rank_mf(st)
        ranks.to_csv(args.output, sep="\t", index=False)

    if args.step == "create-spectral-trees":
        groups = load_groups(args.groups, format="json")
        pls = hdf5_portal.load_peaklists_from_hdf5(args.input)
        spectral_trees = create_spectral_trees(groups, pls)
        save_trees(spectral_trees, args.output, format="json")

    if args.step == "convert-spectral-trees":
        print('converting trees to dimspy peaklists')
        if args.input_type=='json':
            non_merged_pls, merged_pls, ms1_precursor_pl = tree2peaklist(tree_pth=args.input,
                                                                         out_pth=args.output,
                                                                         name=args.name,
                                                                         adjust_mz=args.adjust_mz,
                                                                         merge=args.merge,
                                                                         ppm=args.ppm)
            if args.msp:
                print('Converting dimspy peaklists to MSP files')
                if non_merged_pls:
                    peaklist2msp(non_merged_pls,
                                 os.path.join(args.output, '{}_non_merged.msp'.format(args.name)),
                                 msp_type=args.msp_type,
                                 polarity=args.polarity)
                if merged_pls:
                    peaklist2msp(merged_pls,
                                 os.path.join(args.output, '{}_merged.msp'.format(args.name)),
                                 msp_type=args.msp_type,
                                 polarity=args.polarity)
                if ms1_precursor_pl:
                    peaklist2msp(ms1_precursor_pl,
                                 os.path.join(args.output, '{}_ms1_precursors.msp'.format(args.name)),
                                 msp_type=args.msp_type,
                                 polarity=args.polarity,
                                 include_ms1=True)
        else:
            pls = hdf5_portal.load_peaklists_from_hdf5(args.input)
            peaklist2msp(pls,
                         os.path.join(args.output, '{}.msp'.format(args.name)),
                         msp_type=args.msp_type,
                         polarity=args.polarity)



if __name__ == "__main__":
    main()
