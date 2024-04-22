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
from time import time
from glob import glob
from calendar import monthrange

from math import exp
from pandas import DataFrame, read_csv
from netCDF4 import Dataset
from numpy import array, int32
from numpy.ma import masked as MaskedConstant

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

from mscnfr_utils_fns import update_progress_chess, open_file_sets_cru
from mscnfr_osgb_fns import open_file_sets_osgb

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '

RQRD_METRICS = ['precip', 'hurs', 'huss', 'psurf', 'rsds', 'sfcWind', 'tasmax', 'tasmin', 'tas']

KELVIN_TO_CENTIGRADE = -273.15
TAS_METRICS = ['tasmax', 'tasmin', 'tas']   # temperature at surface

START_YEAR, END_YEAR = 1961, 2080
GRID_RESOL = 0.0083333

numSecsDay = 3600*24

SCENARIOS = ['rcp60', 'rcp45', 'rcp26', 'rcp85']
SCENARIOS = ['rcp26']
REALISATIONS = ['01', '04', '06', '15']
REALISATIONS = ['01']

MSCNFR_FLD_WDTH = 8
YR_BLOCK_LNGTH = MSCNFR_FLD_WDTH*12
FRMT_STR = '{:' + str(MSCNFR_FLD_WDTH) + 'd}'

def _make_mess(bad_cells, metric, yindx, xindx, lat, lon, mess):
    """

    """
    mess += 'metric {}\tyindx, xindx: {}/{}\tlat, lon: {}/{}'.format(metric, yindx, xindx, lat, lon)
    bad_cells.append(mess)
    return

def _relative_humidity_calc(pettmp, tref = 273.16):
    """
    thanks to Jon Mccalmont <jon.mccalmont@abdn.ac.uk> for equation

    spec_hum:   specific humidity (huss kg kg-1)
    tmean_K:    mean temperature (tas, K)
    psurf:      pressure (psurf, Pa),

    """
    for metric in ['huss', 'psurf', 'tas']:
        if metric not in pettmp:
            return

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

def _make_meteo_csvs_from_chess(out_dir, nc_fnames, nc_fnames_rcp, meteo_df, max_cells, rqrd_metrics):
    """
    read y/x indices from lookup table
    """
    num_cells = len(meteo_df)
    num_cells_str = format_string("%d", num_cells, grouping=True)
    print('Coordinates mapping file has ' + num_cells_str + ' + records')

    max_cells = min(max_cells, num_cells)

    metric_names = nc_fnames.keys()
    n_metrics = len(metric_names)

    # ===========================================================================
    miscan_fobjs, writers = open_file_sets_osgb(rqrd_metrics + ['meteogrid'], out_dir, START_YEAR, END_YEAR)

    vals_hist_all = {}
    vals_rcp_all = {}
    for metric in rqrd_metrics:

        if metric == 'hurs':
            t2 = time()
        else:
            # historic weather
            # ================
            t1 = time()
            nc_dset = Dataset(nc_fnames[metric])
            vals_hist_all[metric] = nc_dset.variables[metric][:-1, :, :]
            nc_dset.close()

            t2 = time()
            print('Read {} in {} secs'.format(nc_fnames[metric], int(t2 - t1)))

        # historic finishes at last month of 2017 so we start RCP on first month of 2018
        # ==============================================================================
        if metric == 'precip':
            var_name = 'pr'
        elif metric == 'psurf':
            continue
        else:
            var_name = metric

        # future weather
        # ==============
        nc_dset_rcp = Dataset(nc_fnames_rcp[metric])
        vals_rcp_all[metric] = nc_dset_rcp.variables[var_name][444:, :, :]
        nc_dset_rcp.close()

        t3 = time()
        print('Read {} in {} secs'.format(split(nc_fnames_rcp[metric])[1], int(t3 - t2)))

    # main loop
    # =========
    yindx_max = vals_rcp_all[metric].shape[1] - 1
    xindx_max = vals_rcp_all[metric].shape[2] - 1
    num_total = max_cells*len(metric_names)
    num_set = 0     # a set is is a set of metric values for each grid point
    icount = 0      # counter for each value retrieved for all metrics
    last_time = time()
    bad_cells = []
    for lat, lon, yindx, xindx in zip(meteo_df['cell_lat'], meteo_df['cell_lon'], meteo_df['yindx'], meteo_df['xindx']):
        if yindx > yindx_max or xindx > xindx_max:
            mess = ERROR_STR + 'meteo dataframe index out of bounds - lat, yindx, max: '
            mess += '{} {} {}\tlon, xindx, max: {} {} {}'.format(lat, yindx, yindx_max, lon, xindx, xindx_max)
            print('\n' + mess)
            continue

        pettmp = {}
        for metric in rqrd_metrics:
            if metric == 'hurs':
                continue

            vals_hist = vals_hist_all[metric][:, yindx, xindx]
            pettmp[metric] = [float(val) for val in vals_hist]

            if metric != 'psurf':
                vals_rcp = vals_rcp_all[metric][:, yindx, xindx]
                try:
                    pettmp[metric] = pettmp[metric] + [float(val) for val in vals_rcp]
                except UserWarning as warning:
                    _make_mess(bad_cells, metric, yindx, xindx, lat, lon, WARN_STR + warning)

            last_time = update_progress_chess(last_time, num_set, icount, len(bad_cells), num_total,
                                                                                            round(lat,4), round(lon,4))
            icount += 1

        # stanza for relative humidity
        # ============================
        _relative_humidity_calc(pettmp)
        metric = 'hurs'
        if metric in pettmp:
            vals_rcp = vals_rcp_all[metric][:, yindx, xindx]
            pettmp[metric] = pettmp[metric] + [float(val) for val in vals_rcp]

        # write data for this grid point
        # ==============================
        _write_out_chess(writers, pettmp, yindx, xindx, lat, lon, bad_cells)

        num_set += 1
        if num_set >= max_cells:
            break

    # close all file objects
    # ======================
    for var in miscan_fobjs:
        miscan_fobjs[var].close()

    print('\n\nFinished writing {} datasets comprising {} grid cells\tto: {}\n'.format(n_metrics, num_set, out_dir))

    return bad_cells

def read_chess_nc_write_meteo_csv(form):
    """
    process NCs for required metrics only
    """
    meteogrid_only = True
    if meteogrid_only:
        rqrd_metrics = RQRD_METRICS[:1]

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

    #
    # =========================================================================
    meteo_lat_lon_osgb_fn = form.settings['lat_lon_osgb']
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

            bad_cells = _make_meteo_csvs_from_chess(out_dir_plus, nc_fnames_hist, nc_fnames_rcp,
                                                                                meteo_df, max_cells, rqrd_metrics)
            mess = ' N bad cells: '.format(len(bad_cells))
            print('Completed met files for scenario ' + scenario + ' realisation ' + realisation + mess)

    return

def _write_out_chess(writers, pettmp, yindx, xindx, lat, lon, bad_cells):
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
            _make_mess(bad_cells, metric, yindx, xindx, lat, lon, str(END_YEAR) + ' ' + str(err))
            integrity_flag = False
            break

        try:
            recs[metric] = ''.join(newlist[metric])
        except TypeError or KeyError as err:
            _make_mess(bad_cells, metric, yindx, xindx, lat, lon, END_YEAR + err)
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
