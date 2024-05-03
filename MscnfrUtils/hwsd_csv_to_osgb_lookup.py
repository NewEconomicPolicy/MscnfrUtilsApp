#-------------------------------------------------------------------------------
# Name:        ukcp18_fns.py
# Purpose:     Functions to create Miscanfor formatted metric files from UKCP18 RCP8.5 climate CSV files
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'ukcp18_fns.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import isdir, join, split, splitext
from os import mkdir
from time import time

from pandas import DataFrame, read_csv, read_excel

from mscnfr_utils_fns import update_progress_hwsd
from cvrtcoord import WGS84toOSGB36, OSGB36toWGS84
from mscnfr_osgb_fns import open_file_sets_osgb, write_mscnfr_out

ERROR_STR = '*** Error *** '

numSecsDay = 3600*24

HWSD_CSV_FNAME = 'E:\\GlobalEcosseData\\Hwsd_CSVs\\UK\\vault\\Wales_hwsd.csv'
HWSD_CSV_FNAME =  'E:\\GlobalEcosseData\\Hwsd_CSVs\\UK\\GBR_hwsd.csv'
LKUP_TBL_FN = 'G:\\MiscanforData\\All_outputs\\chess\\CHESS_hwsd_lkup_tble.csv'

AOI_HEADERS = ['gran_lat', 'gran_lon', 'mu_global', 'lat', 'lon']

ESTNG_MAX  = 655500
NRTHNG_MAX = 1056500

# ===================================================
#
def remove_cells_from_hwsd(form):
    """
    create dataframe from BritishGrid_HWSD.csv
    read log file comprising cells which have no CHESS data
    remove these cells from the dataframe and save it as a new CSV file.
    """
    brit_osgb_fn = form.settings['brit_osgb_fn']
    hwsd_df = read_csv(brit_osgb_fn)
    nrecs = len(hwsd_df)
    print('\nAOI HWSD file has {} records'.format(nrecs))

    # read log file
    # =============
    fltr_cells_fn = form.settings['fltr_cells_fn']
    with open(fltr_cells_fn, 'r') as fobj:
        fltr_lines = fobj.readlines()

    nrecs_rmv = len(fltr_lines) - 1
    print('\nWill remove {} cells from new HWSD file'.format(nrecs_rmv))

    eastings = []
    nrthings = []

    NREC_LIM = 999999999
    for nline, line in enumerate(fltr_lines[1:]):
        indices_str = line.split()[11]
        yindx, xindx = indices_str.split('/')
        easting = int(xindx) * 1000 - 500
        nrthing = int(yindx) * 1000 - 500

        eastings.append(easting)
        nrthings.append(nrthing)

    # remove cells
    # ============
    print('removed {} cells from dataframe'.format(len(nrthings)))

    out_dir = form.w_lbl_out_dir.text()
    if not isdir(out_dir):
        mkdir(out_dir)

    return

def _make_lookup_from_hwsd_xlsx(form, aoi_df):
    """
    reorder aoi_df and add index Series
    """
    HWSD_LKUP_HEADERS = ['cell_lat', 'cell_lon', 'northing', 'yindx', 'easting', 'xindx']
    indx_eastngs = []
    indx_nrthngs = []
    for ic, rec in enumerate(aoi_df.values):
        eastng = float(rec[0])
        nrthng = float(rec[1])
        indx_eastngs.append(int(eastng/1000))
        indx_nrthngs.append(int(nrthng/1000))
        # last_time = update_progress_hwsd(last_time, ic, nrecs)

    lkup_df = DataFrame(columns=HWSD_LKUP_HEADERS)

    lkup_df['cell_lat'] = aoi_df['latitude']
    lkup_df['cell_lon'] = aoi_df['longitude']

    lkup_df['northing'] = aoi_df['BNG_Y']
    lkup_df['yindx'] = indx_nrthngs

    lkup_df['easting'] = aoi_df['BNG_X']
    lkup_df['xindx'] = indx_eastngs

    return lkup_df

def write_lookup_from_hwsd_xlsx(form):
    """
    brit_osgb_fn is a CSV file of BritishGrid_HWSD.xlsx
    conversion to CSV speeds up pandas read_csv compared with read_excel
    """
    col_names = ['BNG_X', 'BNG_Y', 'longitude', 'latitude']
    brit_osgb_fn = form.settings['brit_osgb_fn']
    hwsd_df = read_csv(brit_osgb_fn, usecols=col_names)
    nrecs = len(hwsd_df)
    print('\nAOI HWSD file has {} records'.format(nrecs))

    out_dir = form.w_lbl_out_dir.text()
    if not isdir(out_dir):
        mkdir(out_dir)

    lkup_df = _make_lookup_from_hwsd_xlsx(form, hwsd_df)

    out_fn = join(out_dir, 'CHESS_hwsd_lkup_tble.csv')
    lkup_df.to_csv(out_fn, index=False, header=True)
    print('\nWrote ' + out_fn)

    return

def write_osgb_meteogrid(form):
    """

    """

    # TODO: write sorted list
    lkup_df = read_csv(LKUP_TBL_FN)
    lkup_fltr_df = lkup_df[lkup_df.duplicated(['nrthng_indx', 'estng_indx'])]
    ncells = len(lkup_fltr_df)
    print('Lookup dataframe created from ' + LKUP_TBL_FN + ' with length {}'.format(ncells))

    var_names = []
    out_dir = split(LKUP_TBL_FN)[0]
    miscan_fobjs, writers = open_file_sets_osgb(var_names + ['meteogrid'], out_dir)

    pettmp = {}
    nyears = 100
    for indx in range(ncells):
        east_indx = lkup_fltr_df['estng_indx'].values[indx]
        easting = east_indx * 1000 - 500
        nrth_indx = lkup_fltr_df['nrthng_indx'].values[indx]
        nrthing = nrth_indx * 1000 - 500
        pettmp['meteogrid'] = list([easting, nrthing])
        write_mscnfr_out(pettmp, writers, nyears, meteogrid_flag = True)

        lon = lkup_fltr_df['lon'].values[indx]
        lat = lkup_fltr_df['lat'].values[indx]
        lon2, lat2 = OSGB36toWGS84(easting, nrthing)
        pass

    miscan_fobjs['meteogrid'].close()
    print('Created ' + miscan_fobjs['meteogrid'].name)

    return

# ===================================================
#
def _make_lookup_from_hwsd_csv(form, aoi_df, nrecs):
    """

    """
    last_time = time()

    indx_eastngs = []
    indx_nrthngs = []
    for ic, rec in enumerate(aoi_df.values):
        lat = float(rec[3])
        lon = float(rec[4])
        eastng, nrthng = WGS84toOSGB36(lon, lat)
        indx_eastngs.append(int(eastng/1000))
        indx_nrthngs.append(int(nrthng/1000))
        last_time = update_progress_hwsd(last_time, ic, nrecs)

    '''
    indx_nrthngs_s = Series(indx_nrthngs)
    lkup_df = aoi_df.assign(nrthng_indx = indx_nrthngs_s)
    indx_easthngs_s = Series(indx_eastngs)
    aoi_df = lkup_df.assign(eastng_indx = indx_easthngs_s)
    '''
    aoi_df['nrthng_indx'] = indx_nrthngs
    aoi_df['estng_indx'] = indx_eastngs

    return aoi_df

def write_lookup_from_hwsd_csv(form):
    """

    """
    # hwsd_df = HWSD_mu_globals_csv(form, HWSD_CSV_FNAME)
    aoi_df = read_csv(HWSD_CSV_FNAME, sep=',', names=AOI_HEADERS)
    nrecs = len(aoi_df)
    print('\nAOI HWSD file has {} records'.format(nrecs))

    out_dir = join(form.w_lbl_out_dir.text(), 'chess')
    if not isdir(out_dir):
        mkdir(out_dir)

    aoi_df = _make_lookup_from_hwsd_csv(form, aoi_df, nrecs)

    lkup_fltr_df = aoi_df[aoi_df.duplicated(['nrthng_indx', 'estng_indx'])]     # compress dataframe

    out_fn = join(out_dir, 'CHESS_hwsd_lkup_tble.csv')
    print('\nWriting ' + out_fn)
    lkup_fltr_df.to_csv(out_fn, index=False, header=True)

    return