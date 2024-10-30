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

from os.path import isdir, split, exists, join, isfile
from os import mkdir, remove, makedirs
from copy import copy
from time import time
from glob import glob
from calendar import monthrange

from warnings import filterwarnings

from math import exp
from pandas import DataFrame, read_csv
from netCDF4 import Dataset
from numpy.ma.core import MaskedConstant
from numpy.ma import masked as MaskedConstant, isMaskedArray

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

from mscnfr_utils_fns import update_progress_chess, open_file_sets_cru
from mscnfr_osgb_fns import open_file_sets_osgb

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '

KELVIN_TO_CENTIGRADE = -273.15
TAS_METRICS = ['tasmax', 'tasmin', 'tas']   # temperature at surface

START_YEAR, END_YEAR = 1961, 2080
GRID_SIZE = 1000

numSecsDay = 3600*24

#RQRD_METRICS = ['precip']
#SCENARIOS = ['rcp26']
#REALISATIONS = ['01']

RQRD_METRICS = ['precip', 'hurs', 'huss', 'psurf', 'rsds', 'sfcWind', 'tasmax', 'tasmin', 'tas']  # 9
# RQRD_METRICS = ['hurs', 'huss', 'psurf', 'tas']
SCENARIOS = ['rcp60', 'rcp45', 'rcp26', 'rcp85']
REALISATIONS = ['01', '04', '06', '15']

NSEARCH_PTS = 1     # do not exceed maximum of 9 in function _fetch_valid_hist_rcp_data

MSCNFR_FLD_WDTH = 8
YR_BLOCK_LNGTH = MSCNFR_FLD_WDTH*12
FRMT_STR = '{:' + str(MSCNFR_FLD_WDTH) + 'd}'

filterwarnings("error")

def _fetch_valid_hist_rcp_data(lggr, vals_hist_all, vals_rcp_all, metric, pettmp, yindx, xindx, lat, lon):
    """
    simplified version of _fetch_valid_hist_rcp_data_orig which searches in neighbouring cells if no data in
    cell specified by yindx, xindx
    """
    vals_hist = None
    vals_rcp = None

    valid_data_flag = True

    if not metric == 'hurs':        # no hurs historic dataset
        vals_hist = vals_hist_all[:, yindx, xindx]
        if vals_hist[0] is MaskedConstant:
            valid_data_flag = False

    if not metric == 'psurf':       # no psurf rcp dataset
        vals_rcp = vals_rcp_all[:, yindx, xindx]
        if vals_rcp[0] is MaskedConstant:
            valid_data_flag = False

    pettmp[metric] = []
    if vals_hist is not None:
        try:
            pettmp[metric] = [float(val) for val in vals_hist]
        except UserWarning as warn:
            valid_data_flag = False

    if vals_rcp is not None:
        try:
            pettmp[metric] += [float(val) for val in vals_rcp]
            pettmp[metric] += [float(vals_rcp[-1])]     # fudge due to CHESS incompleteness
        except UserWarning as warn:
            valid_data_flag = False

    if not valid_data_flag:
        mess = WARN_STR + 'No valid data for metric {}\tyindx, xindx: {}/{}'.format(metric, yindx, xindx)
        mess += '\tlat, lon: {}/{}'.format(lat, lon)
        lggr.info(mess)

    return valid_data_flag

def _make_meteo_csvs_from_chess(lggr, out_dir, nc_fnames_hist, nc_fnames_rcp, meteo_df, max_cells, rqrd_metrics):
    """
    read y/x indices from lookup table
    """
    nrth_min = meteo_df['northing'].min()
    nrth_max = meteo_df['northing'].max()
    east_min = meteo_df['easting'].min()
    east_max = meteo_df['easting'].max()

    num_cells = len(meteo_df)
    max_cells = min(max_cells, num_cells)

    metric_names = set(list(nc_fnames_rcp.keys()) + list(nc_fnames_hist.keys()))    # required for update_progress_chess
    n_metrics = len(metric_names)

    # ===========================================================================
    miscan_fobjs, writers = open_file_sets_osgb(rqrd_metrics + ['meteogrid'], out_dir, START_YEAR, END_YEAR,
                                                                GRID_SIZE, nrth_min, nrth_max, east_min, east_max)
    vals_hist_all = {}
    vals_rcp_all = {}
    nrqrd = len(rqrd_metrics)
    for metric in rqrd_metrics:
        '''
        historic weather
        ================
        historic finishes at last month of 2017 so we start RCP on first month of 2018
        '''
        if metric == 'hurs':    # no historic data for relative humidity
            t2 = time()
            vals_hist_all[metric] = None
        else:
            t1 = time()
            nc_dset = Dataset(nc_fnames_hist[metric])
            vals_hist_all[metric] = nc_dset.variables[metric][:240, :, :]
            nc_dset.close()

            t2 = time()
            print('Read {} in {} secs'.format(nc_fnames_hist[metric], int(t2 - t1)))
        '''
        future weather
        ==============
        start RCP on first month of 2018 since historic finishes at last month of 2017 
        '''
        if metric == 'precip':
            var_name = 'pr'
        else:
            var_name = metric

        # future weather
        # ==============
        if metric == 'psurf':  # no RCP data for surface air pressure
            vals_rcp_all[metric] = None
        else:
            nc_dset_rcp = Dataset(nc_fnames_rcp[metric])
            vals_rcp_all[metric] = nc_dset_rcp.variables[var_name][1:, :, :]
            nc_dset_rcp.close()

            t3 = time()
            print('Read {} in {} secs'.format(split(nc_fnames_rcp[metric])[1], int(t3 - t2)))

    # main loop
    # =========
    yindx_max = vals_rcp_all[metric].shape[1] - 1
    xindx_max = vals_rcp_all[metric].shape[2] - 1
    num_total = max_cells*len(metric_names)
    num_set = 0     # a set is a set of metric values for each grid point
    icount = 0      # counter for each value retrieved for all metrics
    last_time = time()
    nbad_cells = 0
    bad_cells_list = []
    nmasked = 0
    for lat, lon, yindx, xindx in zip(meteo_df['cell_lat'], meteo_df['cell_lon'], meteo_df['yindx'], meteo_df['xindx']):

        if _out_of_bounds_check(yindx, yindx_max, xindx, xindx_max, lat, lon):
            continue

        pettmp = {}
        for metric in rqrd_metrics:

            valid_data_flag = _fetch_valid_hist_rcp_data(lggr, vals_hist_all[metric], vals_rcp_all[metric],
                                                                            metric, pettmp, yindx, xindx, lat, lon)
            if not valid_data_flag:
                nmasked += 1
                break

            last_time = update_progress_chess(last_time, num_set, icount, nbad_cells, num_total,
                                                                                            round(lat,4), round(lon,4))
            icount += 1

        if not valid_data_flag:
            continue

        # set flag for relative humidity
        # ==============================
        hurs_flag = True
        for metric in ['huss', 'psurf', 'tas']:
            if metric not in rqrd_metrics:
                hurs_flag = False
                break

        # stanza for relative humidity
        # ============================
        if hurs_flag:
            _relative_humidity_calc(pettmp)  # create historic hurs data from tas, psurf and huss
            metric = 'hurs'
            vals_rcp = vals_rcp_all[metric][:, yindx, xindx]
            try:
                pettmp[metric] += [float(val) for val in vals_rcp]
            except UserWarning as warn:
                valid_data_flag = False

        # write data for this grid point
        # ==============================
        # TODO: improve checking
        if len(pettmp) == nrqrd:
            _write_out_chess(writers, pettmp, yindx, xindx, lat, lon, bad_cells_list)
            num_set += 1
            if num_set >= max_cells:
                break
        else:
            mess = 'No CHESS data for HWSD grid cell lat: {}\tlong: {}\tnorthing: {}\teasting: {}'.format(lat, lon,
                                                                                yindx*1000 + 500, xindx*1000 + 500)              #                                                            metric,
            lggr.info(mess)
            nbad_cells += 1

    # close all file objects
    # ======================
    for var in miscan_fobjs:
        miscan_fobjs[var].close()

    print('\nFinished writing {} datasets comprising {} grid cells\n\tto: {}\n'.format(n_metrics, num_set, out_dir))
    if nmasked > 0:
        print(WARN_STR + '{} HWSD grid cells have no CHESS data - see logfile for detail'.format(nmasked))

    return nbad_cells

def read_chess_nc_write_meteo_csv(form):
    """
    process NCs for required metrics only
    """
    meteogrid_only = form.w_dummy.isChecked()       # also used in merge_pt2.py
    if meteogrid_only:
        rqrd_metrics = RQRD_METRICS[:1]
    else:
        rqrd_metrics = []
        for metric in RQRD_METRICS:
            rqrd_metrics.append(metric)

    nc_fnames_hist = {}
    for nc_fname in glob(form.settings['chess_hist'] + '/*.nc'):
        short_fn = split(nc_fname)[1]
        for metric in rqrd_metrics:
            if short_fn.rfind(metric) == 0:
                nc_fnames_hist[metric] = nc_fname
                break

    nfiles = len(nc_fnames_hist)
    if nfiles == 0:
        print('No output files required - nothing to do')
        return
    else:
        print('Will process {} NC files'.format(nfiles))

    out_dir = join(form.w_lbl_out_dir.text(), 'chess_exp')
    if not isdir(out_dir):
        mkdir(out_dir)

    # OSGB file is CHESS_hwsd_lkup_tble.csv derived from BritishGrid_HWSD.xlsx
    # ========================================================================
    meteo_lat_lon_osgb_fn = form.settings['lat_lon_osgb']
    if not isfile(meteo_lat_lon_osgb_fn):
        print(ERROR_STR + 'Look up file ' + meteo_lat_lon_osgb_fn + ' must exist')
        return
    try:
        meteo_df = read_csv(meteo_lat_lon_osgb_fn)
    except PermissionError as err:
        print('Could not read ' + meteo_lat_lon_osgb_fn + ' due to: ' + str(err))
        return
    num_cells = len(meteo_df)
    num_cells_str = format_string("%d", num_cells, grouping=True)
    print('Coordinates mapping file ' + meteo_lat_lon_osgb_fn + ' has ' + num_cells_str + ' + records')

    # =======================================
    try:
        max_cells = int(form.w_max_cells.text())
    except ValueError as err:
        max_cells = 10
    for scenario in SCENARIOS:
        for realisation in REALISATIONS:
            nc_fnames_rcp = _get_rcp_nc_files(form.settings, scenario, realisation, rqrd_metrics)
            out_dir_plus = join(out_dir, scenario, realisation)
            if not isdir(out_dir_plus):
                makedirs(out_dir_plus)
                print('Created ' + out_dir_plus)

            nbad_cells = _make_meteo_csvs_from_chess(form.lgr, out_dir_plus, nc_fnames_hist, nc_fnames_rcp,
                                                                meteo_df, max_cells, rqrd_metrics)
            mess = '\tN bad cells: {}'.format(nbad_cells)
            print('Completed met files for scenario ' + scenario + ' realisation ' + realisation + mess)

    return

def _make_bad_cell_mess(bad_cells_list, metric, yindx, xindx, lat, lon, mess):
    """

    """
    mess += 'metric {}\tyindx, xindx: {}/{}\tlat, lon: {}/{}'.format(metric, yindx, xindx, lat, lon)
    bad_cells_list.append(mess)
    return

def _relative_humidity_calc(pettmp, tref = 273.16):
    """
    purpose: use Clausius-Clapeyron equation to derive historic hurs data since RCP data for hurs is extant
                 thanks to Jon Mccalmont <jon.mccalmont@abdn.ac.uk> for equation

        spec_hum:   specific humidity (huss kg kg-1)    huss  - hist and rcp available    1440 values 120 yrs
        tmean_K:    mean temperature (tas, K)           tas   - hist and rcp available    1440 values 120 yrs
        psurf:      pressure (psurf, Pa)                psurf - hist only    available     684 values  57 yrs

    """
    hurs = []
    for tmean_K, psurf, spec_hum in zip(pettmp['tas'], pettmp['psurf'], pettmp['huss']):

        # Clausius-Clapeyron equation
        # ===========================
        rltv_humid = 0.263 * psurf * spec_hum * (exp((17.67 * (tmean_K - tref)) / (tmean_K - 29.65))) ** -1

        # for values greater than 100 or less than 0
        # ===========================================
        if rltv_humid > 100:
            rltv_humid = 100
        elif rltv_humid < 0:
            rltv_humid = 0

        hurs.append(rltv_humid)

    pettmp['hurs'] = hurs

    return

def _write_out_chess(writers, pettmp, yindx, xindx, lat, lon, bad_cells_list):
    """
    for this coordinate write each variable to a separate file
    metrics are passed as a list of real values which we convert to integers after times by 10
        huss    Near-surface (1.5 m above surface) specific humidity, eg 0.003511
                Specific humidity is the ratio of water vapor mass to total moist air parcel mass.
        precip  Precipitation flux kg m-2 s-1 to mm
        rlds    Surface downwelling longwave radiation W m-2, eg 286.7
        rsds    Surface downwelling shortwave radiation W m-2, eg 24.7
        sfcWind Near-surface (10 m above surface) wind speed m s-1, eg 4.5
        tas     Near-surface (1.5 m above surface) air temperature    K, eg 277.0
    """

    # stanza to create mean temperature
    # =================================
    metrics = list(pettmp.keys())
    if 'tasmean' in metrics:
        pettmp['tasmean'] = []
        pettmp['tasrange'] = []
        for tmax, tmin in zip(pettmp['tasmax'], pettmp['tasmin']):
            pettmp['tasmean'].append((tmax + tmin) / 2.0)
            pettmp['tasrange'].append(tmax - tmin)

    # main loop
    # =========
    integrity_flag = True
    newlist = {}
    recs = {}
    for metric in metrics:
        if metric == 'huss':
            continue

        out_rec = pettmp[metric]

        if metric == 'precip':

            # convert rainfall which has units kg m-2 s-1 to mm
            # =================================================
            yr = START_YEAR - 1
            for imnth, val in enumerate(out_rec):
                this_mnth = imnth % 12
                if this_mnth == 0:
                    yr += 1
                ndays = monthrange(yr, this_mnth + 1)[1] #
                out_rec[imnth] = val * numSecsDay * ndays

        if metric in TAS_METRICS:
            # convert temperature from K to C degrees
            # =======================================
            for imnth, val in enumerate(out_rec):
                out_rec[imnth] = val + KELVIN_TO_CENTIGRADE

        # convert values to integers multiplied by 10
        # ===========================================
        try:
            newlist[metric] = [FRMT_STR.format(int(10.0 * val)) for val in out_rec]
        except ValueError as err:
            _make_bad_cell_mess(bad_cells_list, metric, yindx, xindx, lat, lon, str(END_YEAR) + ' ' + str(err))
            integrity_flag = False
            break

        try:
            recs[metric] = ''.join(newlist[metric])
        except TypeError or KeyError as err:
            _make_bad_cell_mess(bad_cells_list, metric, yindx, xindx, lat, lon, str(END_YEAR) + ' ' + str(err))
            integrity_flag = False
            break

    if integrity_flag:
        for metric in metrics:
            if metric == 'huss':
                continue

            ntsteps = len(newlist[metric])
            for indx in range(0, ntsteps, 12):
                indx_str = indx*MSCNFR_FLD_WDTH
                writers[metric].writerow(list([recs[metric][indx_str:indx_str + YR_BLOCK_LNGTH]]))

        easting = xindx * 1000 - 500
        nrthing = yindx * 1000 - 500
        writers['meteogrid'].writerow(list(['{:10d}{:10d}'.format(easting, nrthing)]))

    return

def _out_of_bounds_check(yindx, yindx_max, xindx, xindx_max, lat, lon):
    """

    """
    if yindx > yindx_max or xindx > xindx_max:
        mess = ERROR_STR + 'meteo dataframe index out of bounds - lat, yindx, max: '
        mess += '{} {} {}\tlon, xindx, max: {} {} {}'.format(lat, yindx, yindx_max, lon, xindx, xindx_max)
        print('\n' + mess)
        out_of_bounds_flag = True
    else:
        out_of_bounds_flag = False

    return  out_of_bounds_flag

def _get_rcp_nc_files(settings, scenario, realisation, rqrd_metrics):
    """

    """
    wthr_dir = join(settings['chess_rcp'], scenario + '_bias-corrected', realisation)
    all_nc_fns = glob(wthr_dir + '/*.nc')

    # list required NCs only
    # ======================
    nc_fnames = {}
    for nc_fname in all_nc_fns:
        short_fn = split(nc_fname)[1]
        for metric in rqrd_metrics:
            if metric == 'precip':
                srch_str = 'pr'
            else:
                srch_str = metric

            nloc = short_fn.rfind(srch_str + '_')
            if nloc >= 0:
                nc_fnames[metric] = nc_fname
                break

    return nc_fnames

def _fetch_valid_hist_rcp_data_orig(lggr, vals_hist_all, vals_rcp_all, metric, pettmp, yindx, xindx, lat, lon):
    """
    try initial location and if no valid data, then look in neighbouring cells using the lookup table
    defined by ydsp and xdsp
    """
    ydsp = [0,  0, 1, 0, -1, -1,  1, 1, -1]
    xdsp = [0, -1, 0, 1,  0, -1, -1, 1,  1]

    vals_hist = None
    vals_rcp = None

    for ic in range(NSEARCH_PTS):
        valid_data_flag = True
        yindex = yindx + ydsp[ic]
        xindex = xindx + xdsp[ic]

        if not metric == 'hurs':        # no hurs historic dataset
            vals_hist = vals_hist_all[:, yindex, xindex]
            if vals_hist[0] is MaskedConstant:
                valid_data_flag = False

        if not metric == 'psurf':       # no psurf rcp dataset
            vals_rcp = vals_rcp_all[:, yindex, xindex]
            if vals_rcp[0] is MaskedConstant:
                valid_data_flag = False

        if valid_data_flag:
            break

    pettmp[metric] = []
    if vals_hist is not None:
        try:
            pettmp[metric] = [float(val) for val in vals_hist]
        except UserWarning as warn:
            valid_data_flag = False

    if vals_rcp is not None:
        try:
            pettmp[metric] += [float(val) for val in vals_rcp]
        except UserWarning as warn:
            valid_data_flag = False

    if not valid_data_flag:
        mess = WARN_STR + 'No valid data for metric {}\tyindx, xindx: {}/{}'.format(metric, yindx, xindx)
        mess += '\tlat, lon: {}/{}'.format(lat, lon)
        lggr.info(mess)

    return valid_data_flag