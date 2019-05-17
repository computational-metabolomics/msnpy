#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
from six import iteritems
from msnpy.portals import load_trees
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_alignment import align_peaks
import itertools
import numpy as np
import re
import os


def tree2peaklist(tree_pth, out_pth, name='', adjust_mz=True, merge=True, ppm=5, ms1=True):
    ####################################################################################################################
    # Extract peaklists from msnpy
    ####################################################################################################################
    trees = load_trees(tree_pth)
    plsd = {}
    all_ms1_precursors = {}

    # get peaklist for each header
    for tree in trees:
        plsd[tree.graph['id']] = []

        # For each tree we look at each "header" e.g. the same mass spectrometry data (processed prior by dimspy-msnpy)
        # And create a peaklist for each header. (....probably a better way of doing this perhaps iterating through
        # the tree instead?). Anyway this seems to work OK.
        its = tree.nodes.items()
        # add id to tree values
        [i[1].update({'id': i[0]}) for i in its]
        tv = [i[1] for i in its]
        # requires sorting for itertools.groupby to work properly
        tv = sorted(tv, key=lambda i: i['header'])

        for header, group in itertools.groupby(tv, key=lambda x: x['header']):
            # get mz, intensity, mass, molecular formula, adduct

            mtch = re.search('.*Full ms .*', header)
            if mtch:
                # full scan
                continue

            precursor_detail_track = []

            mz = []
            intensity = []
            mass = []
            mf = []
            adduct = []

            metad = {'tree_id': tree.graph['id'], 'header': header, 'parent': {}}

            for d in list(group):

                # get precursor details for each level
                for n in tree.predecessors(d['id']):

                    pd = tree.nodes.get(n)

                    # check if we already have this precursor details
                    if pd['mslevel'] in precursor_detail_track:
                        continue

                    mf_details = pd['mf'].values()[0]
                    metad['parent'][pd['mslevel']] = {}
                    metad['parent'][pd['mslevel']]['mz'] = pd['mz']
                    metad['parent'][pd['mslevel']]['mass'] = mf_details['mass']
                    metad['parent'][pd['mslevel']]['adduct'] = mf_details['adduct']
                    metad['parent'][pd['mslevel']]['mf'] = mf_details['mf']
                    precursor_detail_track.append(pd['mslevel'])

                    if ms1:
                        if adjust_mz:
                            all_ms1_precursors[mf_details['mass']] = pd['intensity']
                        else:
                            all_ms1_precursors[pd['mz']] = pd['intensity']

                mz.append(d['mz'])
                intensity.append(d['intensity'])
                mf_details = d['mf'].values()[0]
                mass.append(mf_details['mass'])
                mf.append(mf_details['mf'])
                adduct.append(mf_details['adduct'])

            if adjust_mz:
                mza = mass
            else:
                mza = mz

            # create dimspy array object
            pl = PeakList(ID='{}: {}'.format(tree.graph['id'], header),
                          mz=mza,
                          intensity=intensity,
                          **metad)
            pl.add_attribute('mass', mass)
            pl.add_attribute('mz_original', mz)
            pl.add_attribute('mf', mf)
            pl.add_attribute('adduct', adduct)

            plsd[tree.graph['id']].append(pl)

    pls = [y for x in list(plsd.values()) for y in x]

    if out_pth:
        save_peaklists_as_hdf5(pls, os.path.join(out_pth, '{}_non_merged_pls.hdf5'.format(name)))


    # Merge
    if merge:
        merged_pls = []
        for (key, pls) in iteritems(plsd):
            merged_id = "<#>".join([pl.ID for pl in pls])
            pm = align_peaks(pls, ppm=ppm)
            plm = pm.to_peaklist(ID=merged_id)
            plm.metadata['parent'] = {1: pls[0].metadata['parent'][1]}

            merged_pls.append(plm)

        if out_pth:
            save_peaklists_as_hdf5(plms, os.path.join(out_pth, '{}_merged_pls.hdf5'.format(name)))
    else:
        merged_pls = ''

    if ms1:
        ms1_precursors_pl = PeakList(ID='ms1_precursors',
                                     mz=list(all_ms1_precursors.keys()),
                                     intensity=list(all_ms1_precursors.values()))
        save_peaklists_as_hdf5(ms1_precursors_pl, os.path.join(out_pth, '{}_ms1_precursors_pls.hdf5'.format(name)))
    else:
        ms1_precursors_pl = ''

    return pls, merged_pls, ms1_precursors_pl


def peaklist2msp(pls, out_pth, msp_type='massbank', polarity='positive', msnpy_annotations=True):

    msp_params = {}

    if msp_type == 'massbank':
        msp_params['name'] = 'RECORD_TITLE:'
        msp_params['polarity'] = 'AC$MASS_SPECTROMETRY: ION_MODE'
        msp_params['precursor_mz'] = 'MS$FOCUSED_ION: PRECURSOR_M/Z '
        msp_params['precursor_type'] = 'MS$FOCUSED_ION: PRECURSOR_TYPE'
        msp_params['num_peaks'] = 'PK$NUM_PEAK:'
        msp_params['cols'] = 'PK$PEAK: m/z int. rel.int.'
        msp_params['ms_level'] = 'AC$MASS_SPECTROMETRY: MS_TYPE '
        msp_params['resolution'] = 'AC$MASS_SPECTROMETRY: RESOLUTION '
        msp_params['fragmentation_mode'] = 'AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE'
        msp_params['collision_energy'] = 'AC$MASS_SPECTROMETRY: COLLISION_ENERGY'
        msp_params['mf'] = 'CH$FORMULA:'
    else:
        msp_params['name'] = 'NAME:'
        msp_params['polarity'] = 'POLARITY:'
        msp_params['precursor_mz'] = 'PRECURSOR_MZ:'
        msp_params['precursor_type'] = 'PRECURSOR_TYPE:'
        msp_params['num_peaks'] = 'Num Peaks:'
        msp_params['cols'] = ''
        msp_params['ms_level'] = 'MS_LEVEL:'
        msp_params['resolution'] = 'RESOLUTION:'
        msp_params['fragmentation_mode'] = 'FRAGMENTATION_MODE:'
        msp_params['mf'] = 'MOLECULAR_FORMULA:'

    with open(out_pth, "w+") as f:
        # Loop through peaklist
        for pl in pls:
            dt = pl.dtable[pl.flags]
            if dt.shape[0] == 0:
                continue

            if re.search('.*Full ms .*', pl.ID) and ms_level > 1:
                continue

            f.write('{} {}\n'.format(msp_params['name'], pl.ID))
            f.write('{} {}\n'.format(msp_params['polarity'], polarity))

            if msnpy_annotations:
                parent_metadata = pl.metadata['parent'][min(pl.metadata['parent'].keys())]

                f.write('{} {}\n'.format(msp_params['precursor_type'], parent_metadata['adduct']))
                f.write('{} {}\n'.format(msp_params['precursor_mz'], parent_metadata['mz']))
                f.write('{} {}\n'.format(msp_params['mf'], parent_metadata['mf']))

            else:
                mtch = re.search('.*Full ms(\d+).*', pl.ID)
                if mtch:
                    f.write('{} {}\n'.format(msp_params['ms_level'], mtch.group(1)))

                mtch = re.search('.*Full ms\d+ (.*) \[.*', pl.ID)
                if mtch:
                    dl = mtch.group(1).split(" ")
                    # get the last detail
                    detail = dl[-1]
                    mtch = re.search('(\d+.\d+)@(\D+)(.*)', detail)
                    if mtch:
                        f.write('{} {}\n'.format(msp_params['precursor_mz'], mtch.group(1)))


            mtch = re.findall('\d+.\d+@(\D+)(\d+.\d+)', pl.ID)
            if mtch:
                mtchz = list(zip(*mtch))
                f.write('{} {}\n'.format(msp_params['fragmentation_mode'], ', '.join(set(mtchz[0]))))
                f.write('{} {}\n'.format(msp_params['collision_energy'], ', '.join(set(mtchz[1]))))

            f.write('{} {}\n'.format(msp_params['num_peaks'], dt.shape[0]))
            if msp_params['cols']:
                f.write('{}\n'.format(msp_params['cols']))

            mz = dt['mz']
            intensity = dt['intensity']
            ra = dt['intensity'] / np.max(dt['intensity']) * 100

            if msp_type=='massbank':
                if 'mf' in dt.dtype.names:
                    mf = dt['mf']
                    adduct = dt['adduct']
                    mass = dt['mass']
                    f.write('PK$ANNOTATION: m/z tentative_formula formula_count adduct\n')
                else:
                    continue
                for i in range(0, len(mz)):
                    mzi = mz[i]
                    mfi = mf[i]
                    massi = mass[i]
                    adducti = adduct[i]
                    if msp_type == 'massbank':
                        f.write('{}\t{}\t{}\n'.format(mzi, mfi, massi, adducti))
                    else:
                        f.write('{}\t{}\n'.format(mzi, rai))

            for i in range(0, len(mz)):
                mzi = mz[i]
                intensityi = intensity[i]
                rai = ra[i]
                if msp_type == 'massbank':
                    f.write('{}\t{}\t{}\n'.format(mzi, intensityi, rai))
                else:
                    f.write('{}\t{}\n'.format(mzi, rai))
            f.write('\n')
