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

import sys
from os.path import split, isfile
from time import time
from locale import format_string

ngranularity = 120

timeElapsed = 5
hdrs_chan17 = ['longHW', 'LatHW', 'id', 'REF_DEPTH', 'SHARE', 'T_REF_BULK', 'T_CACO3', 'fcHW', 'pwpHW', 'socHW']

hdrs_chan11 = ['lon', 'lat', 'dumm02', 'dumm03', 'dumm04', 'lu_ar', 'lu_pa', 'lu_su', 'lu_bu', 'cch', 'dumm10',
               'dumm11', 'dumm12', 'dumm13', 'dumm14', 'dumm15', 'dumm16', 'dumm17', 'dumm18', 'dumm19', 'metelevn',
               'gridelevn', 'lucode', 'SWR']

def _update_progress(last_time, nchecked, nbad, num_total):

    """Update progress bar."""
    this_time = time()
    if (this_time - last_time) > 3.0:
        sys.stdout.flush()
        nbad_str = format_string("%d", nbad, grouping=True)
        ncheck_str = format_string("%d", nchecked, grouping=True)
        percent_done = round(100.0*(1.0 - (num_total - nchecked - nbad)/num_total),1)
        sys.stdout.write('\rCells checked: {} bad: {} {:=5.1f}%'.format(ncheck_str, nbad_str, percent_done))
        return this_time

    return last_time


def _check_rec(rec):

    """ validate record """
    for sval in rec:
        try:
            float(sval)
        except ValueError as e:
            return False

    try:
        int(rec[14])
    except ValueError as e:
        return False

    try:
        int(rec[-1])
    except ValueError as e:
        return False

    return True

def check_soil_file(form):

    '''
    read soildata_HWSD file and inputs
    '''
    soil_fname = 'F:\\MiscanforData\\Problem 2020-04-17\\MISCANFOR-HYBRIDUKERCnewguiEVAP_RCP8.5HWSD\\MiscanforMay2012\\'
    soil_fname += 'input\\soil_data\\HWSDlu_UKSEQ1.txt'
    if not isfile(soil_fname):
        print('soil HWSD file ' + soil_fname + ' does not exist')
        return

    dumyy, short_fname = split(soil_fname)
    headers = list(['long', 'lat', 'id', 'ref_depth', 'share', 't_ref_bulk', 't_caco3', 'fc', 'pwp', 'soc', 'lu_ar', \
                    'lu_pa', 'lu_su', 'lu_bu', 'cch', 'metelevn', 'gridelevn', 'swr'])

    # convetional method
    # ================
    last_time = time()
    with open(soil_fname, 'r') as fobj:
        lines = fobj.readlines()

    num_recs = len(lines)
    bad_lines = []
    nbad = 0
    for indx, line in enumerate(lines[5:]):
        rec = line.split()
        nvals = len(rec)
        if nvals != 18:
            nbad += 1
            bad_lines.append(indx + 1)
        else:
            if not _check_rec(rec):
                nbad += 1
                bad_lines.append(indx + 1)

        last_time = _update_progress(last_time, indx, nbad, num_recs)
    '''
    # dataframe method
    # ================
    soil_df = read_csv(soil_fname, skiprows = 5, names = headers)
    num_recs = len(soil_df)

    # main loop - step through dataframe
    # ==================================
    last_time = time()
    for indx, gran_lat_lon in enumerate(soil_df['gran_lat_lon']):

        lat = soil_df['lat'][indx]
        lon = soil_df['lon'][indx]

        last_time = _update_progress(last_time, indx, num_recs)
    '''

    print('\nFinished after validating {}\tbad records: {}\tfrom total: {}'.format(short_fname, nbad, num_recs))
    return
