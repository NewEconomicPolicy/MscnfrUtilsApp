#-------------------------------------------------------------------------------
# Name:        merge_pt2.py
# Purpose:     Functions to replace merge_pt2
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'merge_pt2.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import isdir, split, join, isfile, splitext
from os import mkdir, remove
from datetime import timedelta
from time import time
from glob import glob
from pandas import DataFrame, read_csv
from mscnfr_utils_fns import update_progress_merge, update_progress_post, check_coord_bbox
from math import sqrt
from shape_funcs import format_bbox
'''
from geopandas import GeoDataFrame
from numpy import array
from scipy.spatial import cKDTree
from shapely.geometry import Point
'''
ngranularity = 120

timeElapsed = 5
hdrs_chan17 = ['longHW', 'LatHW', 'id', 'REF_DEPTH', 'SHARE', 'T_REF_BULK', 'T_CACO3', 'fcHW', 'pwpHW', 'socHW']

hdrs_chan11 = ['lon', 'lat', 'dumm02', 'dumm03', 'dumm04', 'lu_ar', 'lu_pa', 'lu_su', 'lu_bu', 'cch', 'dumm10',
               'dumm11', 'dumm12', 'dumm13', 'dumm14', 'dumm15', 'dumm16', 'dumm17', 'dumm18', 'dumm19', 'metelevn',
               'gridelevn', 'lucode', 'SWR']

def _hwsd_data_frame(hwsd_fname, num_set, max_lines, lgr):
    '''
    in this file lats run from North to South and lons run from West to East
    '''

    # process HWSD file
    # =================
    hwsd_dict = {'gran_lat_lon':[], 'lon':[], 'lat':[], 'id':[], 'ref_depth':[], 'share':[],
                                                't_ref_bulk':[], 't_caco3':[], 'fchw':[], 'pwphw':[], 'sochw':[]}
    with open(hwsd_fname,'r') as fobj:
        lines = fobj.readlines()

    num_total = len(lines)
    print('\nHave read {} lines from {}'.format(num_total, hwsd_fname))

    num_nodata = 0
    last_time = time()
    for nline, line in enumerate(lines):
        if nline >= max_lines:
            break

        srec = line.split()
        rec = [float(sval) for sval in srec]

        lon, lat = rec[:2]
        hwsd_dict['lon'].append(lon)
        hwsd_dict['lat'].append(lat)
        gran_lon = int(round((180.0 + lon)*ngranularity))
        gran_lat = int(round((90.0  - lat)*ngranularity))
        gran_lat_lon = '{}_{}'.format(gran_lat,gran_lon)
        hwsd_dict['gran_lat_lon'].append(gran_lat_lon)

        id, ref_depth, share, t_ref_bulk, t_caco3, fchw, pwphw, sochw = rec[2:]
        hwsd_dict['id'].append(id)
        hwsd_dict['ref_depth'].append(ref_depth)
        hwsd_dict['share'].append(share)
        hwsd_dict['t_ref_bulk'].append(t_ref_bulk)
        hwsd_dict['t_caco3'].append(t_caco3)
        hwsd_dict['fchw'].append(fchw)
        hwsd_dict['pwphw'].append(pwphw)
        hwsd_dict['sochw'].append(sochw)

        last_time = update_progress_post(last_time, num_nodata, nline, num_total, round(lat,3), round(lon,3), num_set)

    hwsd_df = DataFrame(hwsd_dict)
    hwsd_bbox = [hwsd_df['lon'].min(), hwsd_df['lat'].min(), hwsd_df['lon'].max(), hwsd_df['lat'].max()]

    mess = '\nCreated dataframe with {} records for {}'.format(len(hwsd_df), hwsd_fname)
    mess += '\n\tbounding box - ' + format_bbox(hwsd_bbox)
    print(mess)
    lgr.info(mess)

    return hwsd_df, hwsd_bbox

def _soil_data_frame(soil_fname, lgr, max_lines = 99999999999):
    '''
    in this file lats run from North to South and lons run from West to East
    '''
    soil_dict = {'gran_lat_lon':[], 'lon':[], 'lat':[], 'lu_ar':[], 'lu_pa':[], 'lu_su':[], 'lu_bu':[],
                                                    'cch':[], 'metelevn':[], 'gridelevn':[], 'lucode':[], 'SWR':[] }

    with open(soil_fname,'r') as fobj:
        lines = fobj.readlines()

    num_total = len(lines)
    print('\nHave read {} lines from {}'.format(len(lines), soil_fname))

    num_nodata = 0
    last_time = time()
    for nline, line in enumerate(lines):
        if nline >= max_lines:
            break

        srec = line.split()
        rec = [float(sval) for sval in srec]

        lon, lat = rec[:2]
        soil_dict['lon'].append(lon)
        soil_dict['lat'].append(lat)
        gran_lon = int(round((180.0 + lon)*ngranularity))
        gran_lat = int(round((90.0  - lat)*ngranularity))
        gran_lat_lon = '{}_{}'.format(gran_lat,gran_lon)
        soil_dict['gran_lat_lon'].append(gran_lat_lon)

        lu_ar, lu_pa, lu_su, lu_bu, cch = rec[5:10]
        metelevn, gridelevn, lucode, SWR = rec[20:]

        soil_dict['lu_ar'].append(lu_ar)
        soil_dict['lu_pa'].append(lu_pa)
        soil_dict['lu_su'].append(lu_su)
        soil_dict['lu_bu'].append(lu_bu)
        soil_dict['cch'].append(cch)
        soil_dict['metelevn'].append(metelevn)
        soil_dict['gridelevn'].append(gridelevn)
        soil_dict['lucode'].append(lucode)
        soil_dict['SWR'].append(SWR)

        last_time = update_progress_post(last_time, num_nodata, nline, num_total, round(lat,3), round(lon,3))


    # create set of points to assist with location
    # ============================================
    print('\nCreating point series with {} lines...'.format(nline))
    soil_dict['point'] = [(y, x) for y, x in zip(soil_dict['lat'], soil_dict['lon'])]
    soil_df = DataFrame(soil_dict)
    soil_bbox = [soil_df['lon'].min(), soil_df['lat'].min(), soil_df['lon'].max(), soil_df['lat'].max()]

    mess = '\nCreated dataframe with {} records for {}'.format(len(soil_df), soil_fname)
    mess += '\n\tbounding box - ' + format_bbox(soil_bbox)
    print(mess)
    lgr.info(mess)

    return soil_df, soil_bbox

def write_merge_pt2(form):

    '''
    read soildata_HWSD file and inputs
    '''
    lgr = form.lgr
    ignore_near = form.w_ignore_near.isChecked()
    dummy_run = form.w_dummy.isChecked()
    if dummy_run:
        output_flag = False
        out_fname = ''
    else:
        output_flag = True
    last_seq_flag = form.w_last_seq.isChecked()

    data_dir = 'E:\\MiscanforData\\MergePt2'
    fnames = glob(data_dir + '\\soildat*.txt')
    if len(fnames) == 0:
        print('No soil file names in ' + data_dir + ' cannot proceed')
        return

    max_lines = int(form.w_max_cells.text())
    soil_fname = fnames[0]

    # setup output staging area
    # =========================
    out_dir = data_dir + '\\outputs'
    if not isdir(out_dir):
        mkdir(out_dir)

    header = '  {:^10s}{:^10s}{:^12s}{:^8s}{:^10s}'.format('long','lat','id', 'REF_DEPTH', 'SHARE')
    header += '{:^9s}{:^9s}{:^9s}{:^9s}{:^9s}'.format('T_REF_BULK', 'T_CACO3', 'fc', 'pwp', 'soc')
    header += '{:^10s}{:^10s}{:^10s}{:^10s}'.format( 'lu_ar', 'lu_pa', 'lu_su', 'lu_bu')
    header += '{:^4s}{:^10s}{:^10s}{:^3s}'.format('cch', 'metelevn', 'gridelevn', 'SWR')

    # collect HWSD files
    # ==================
    data_dir += '\\inputs'
    hwsd_fnames = glob(data_dir + '\\soilHWSD*.txt')
    hwsd_fnames.reverse()
    if len(hwsd_fnames) == 0:
        print('No HWSD soil file names in ' + data_dir + ' cannot proceed')
        return

    # soildata HWSD file - 416674 lines for UK
    # ========================================
    soil_df, soil_bbox = _soil_data_frame(soil_fname, lgr)

    # main loop - read each sequence HWSD file
    # ========================================
    last_hwsd_fname = hwsd_fnames[-1]
    for hwsd_fname in hwsd_fnames:

        long_name, dummy = splitext(hwsd_fname)
        seq_no = long_name[-1]
        if last_seq_flag:
            if hwsd_fname != last_hwsd_fname:
                continue

        last_time = time()

        if output_flag:
            # open output file and write header
            # =================================
            out_fname = join(out_dir, 'HWSDlu_UKSEQ' + seq_no + '.txt')
            if isfile(out_fname):
                remove(out_fname)
                print('Deleted ' + out_fname)

            fobj_out = open(out_fname, 'w')
            fobj_out.write(header + '\n')

        # read HWSD sequence file
        # =======================
        dummy, hwsd_short_fname = split(hwsd_fname)
        hwsd_df, hwsd_bbox = _hwsd_data_frame(hwsd_fname, int(seq_no), max_lines, lgr)
        num_recs = len(hwsd_df)

        # step through HWSD sequence file
        # ===============================
        num_match = 0
        num_near = 0
        num_outside = 0
        num_written = 0
        start = time()
        for hindx, gran_lat_lon in enumerate(hwsd_df['gran_lat_lon']):

            lat = hwsd_df['lat'][hindx]
            lon = hwsd_df['lon'][hindx]
            rec_write = True

            if not check_coord_bbox(lat, lon, soil_bbox):
                lgr.info('Coordinate with lat: {}\tlon: {} is outside bounding box'.format(lat, lon))
                num_outside += 1
            else:
                soil_locat = soil_df.loc[(soil_df['gran_lat_lon'] == gran_lat_lon)]

                if soil_locat.empty:
                    if output_flag:
                        if ignore_near:
                            rec_write = False
                        else:
                            soil_df['cdist'] = [( sqrt((lat - pnt[0])**2 + (lon - pnt[1])**2) ) for pnt in soil_df['point']]
                            soil_indx  = soil_df['cdist'].idxmin()
                            lu_ar = soil_df['lu_ar'].values[soil_indx]
                            lu_pa = soil_df['lu_pa'].values[soil_indx]
                            lu_su = soil_df['lu_su'].values[soil_indx]
                            lu_bu = soil_df['lu_bu'].values[soil_indx]
                            cch   = soil_df['cch'].values[soil_indx]
                            metelevn  = soil_df['metelevn'].values[soil_indx]
                            gridelevn = soil_df['gridelevn'].values[soil_indx]
                            swr   = soil_df['SWR'].values[soil_indx]

                            soil_lat = soil_df['lat'].values[soil_indx]
                            soil_lon = soil_df['lon'].values[soil_indx]
                            lgr.info('No match for lat: {}\tlon: {}\tnearest found - lat: {}\tlon: {}'
                                                                            .format(lat, lon, soil_lat, soil_lon))
                    num_near += 1
                else:
                    if output_flag:
                        lu_ar = soil_locat['lu_ar'].values[0]
                        lu_pa = soil_locat['lu_pa'].values[0]
                        lu_su = soil_locat['lu_su'].values[0]
                        lu_bu = soil_locat['lu_bu'].values[0]
                        cch   = soil_locat['cch'].values[0]
                        metelevn  = soil_locat['metelevn'].values[0]
                        gridelevn = soil_locat['gridelevn'].values[0]
                        swr   = soil_locat['SWR'].values[0]

                    num_match += 1

                if output_flag and rec_write:
                    out_rec = '{:10.5f}{:10.5f}'.format(lon, lat)
                    out_rec += '{:12d}{:8.1f}{:10.4f}'.format(int(hwsd_df['id'][hindx]),
                                                                    hwsd_df['ref_depth'][hindx], hwsd_df['share'][hindx])

                    out_rec += '{:9.2f}{:9.2f}{:9.2f}{:9.2f}{:9.2f}'.format(hwsd_df['t_ref_bulk'][hindx],
                        hwsd_df['t_caco3'][hindx], hwsd_df['fchw'][hindx], hwsd_df['pwphw'][hindx], hwsd_df['sochw'][hindx])

                    out_rec += '{:10.7f}{:10.7f}{:10.7f}{:10.7f}'.format(lu_ar, lu_pa, lu_su, lu_bu)

                    out_rec += '{:4d}{:10.2f}{:10.2f}{:3d}'.format(int(cch), metelevn, gridelevn, int(swr))

                    fobj_out.write(out_rec + '\n')
                    num_written += 1

            last_time = update_progress_merge(last_time, num_outside, num_near, num_match, num_written, num_recs, int(seq_no))
            '''
            write(25,'(f10.5, f10.5, i12, f8.1, f10.4, f9.2, f9.2, f9.2, f9.2, f9.2, f10.7, f10.7, f10.7, f10.7, i4, f10.2, f10.2, i3)')
            longHW, LatHW, id, REF_DEPTH, SHARE, T_REF_BULK, T_CACO3, fcHW, pwpHW, socHW, lu_ar, lu_pa, lu_su, lu_bu, cch, metelevn, gridelevn, SWR
            '''
        end = time()
        time_taken = str(timedelta(0, int(end - start)))
        mess = '\tRead {} records from {}'.format(num_recs, hwsd_fname)
        mess += '\n\tWrote {} records to {}'.format(num_written, out_fname)
        mess += '\tnum outside: {}\tnear: {}\tmatched: {}\t{}'.format(num_outside, num_near, num_match, time_taken)
        lgr.info(mess)
        print('\n' + mess)

        if output_flag:
            fobj_out.close()

    print('Finished after processing {} HWSD files'.format(len(hwsd_fnames)))
    return
