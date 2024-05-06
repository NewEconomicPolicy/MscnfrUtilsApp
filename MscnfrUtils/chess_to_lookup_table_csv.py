#-------------------------------------------------------------------------------
# Name:        chess_to_lookup_table_csv.py
# Purpose:
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'chess_to_lookup_table_csv.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import isdir, split, exists, join, isfile
from os import mkdir, remove, makedirs
from time import time
from glob import glob
from copy import copy

from pandas import DataFrame
from netCDF4 import Dataset
from numpy import array, int32
from numpy.ma import masked as MaskedConstant

from mscnfr_utils_fns import update_progress_post
from cvrtcoord import WGS84toOSGB36, OSGB36toWGS84

ERROR_STR = '*** Error *** '

numSecsDay = 3600*24

# ========================
#
def _write_check_csv(out_dir, chess_nc_fn, settings, max_cells):
    """

    """
    out_fn = join(out_dir, 'cvrtcoord.csv')

    nc_dset = Dataset(chess_nc_fn, mode='r')
    if 'x' not in nc_dset.variables:
        print(ERROR_STR + 'variable x must be in ' + chess_nc_fn)
        return None

    eastngs = [int32(x) for x in array(nc_dset.variables['x'])]
    nrthngs = [int32(y) for y in array(nc_dset.variables['y'])]
    num_total = len(eastngs) * len(nrthngs)

    lats = []
    y_coords = []
    y_indices = []

    lons = []
    x_coords = []
    x_indices = []

    eastngs_calc = []
    nrthngs_calc = []
    lons_calc = []
    lats_calc = []

    ic = 0
    num_vals = 0
    num_nans = 0
    last_time = time()
    varname = 'huss'
    for yndx, nrthng in enumerate(nrthngs):

        for xndx, eastng in enumerate(eastngs):

            lat = nc_dset.variables['lat'][yndx][xndx]
            lon = nc_dset.variables['lon'][yndx][xndx]
            val = nc_dset.variables[varname][0][yndx][xndx]
            if val is MaskedConstant:
                num_nans += 1
            else:
                lats.append(lat)
                y_coords.append(int(nc_dset.variables['y'][yndx]))
                y_indices.append(yndx)

                lons.append(lon)
                x_coords.append(int(nc_dset.variables['x'][xndx]))
                x_indices.append(xndx)

                eastng_calc, nrthng_calc = WGS84toOSGB36(lon, lat)
                eastngs_calc.append(eastng_calc)
                nrthngs_calc.append(nrthng_calc)

                lon_calc, lat_calc = OSGB36toWGS84(eastng, nrthng)
                lons_calc.append(lon_calc)
                lats_calc.append(lat_calc)

                num_vals += 1

            ic += 1

            last_time = update_progress_post(last_time, num_nans, num_vals, num_total, round(lat, 4), round(lon, 4))

            if num_vals > max_cells:
                break

        if num_vals >= max_cells:
            print('\nnumber of vals: {}\texceeds requested: {}'.format(num_vals, max_cells))
            break

    nc_dset.close()

    meteo_df = DataFrame()
    meteo_df['cell_lat'] = lats
    meteo_df['cell_lon'] = lons
    meteo_df['northing'] = y_coords
    meteo_df['yindx'] = y_indices
    meteo_df['easting'] = x_coords
    meteo_df['xindx'] = x_indices

    meteo_df['lats_calc'] = lats_calc
    meteo_df['lons_calc'] = lons_calc
    meteo_df['nrthngs_calc'] = nrthngs_calc
    meteo_df['eastngs_calc'] = eastngs_calc

    meteo_df = meteo_df.sort_values(by=['cell_lat', 'cell_lon'], ascending=[False, True])

    try:
        meteo_df.to_csv(out_fn, index=False, header=True)
    except PermissionError as err:
        print('Could not create ' + out_fn + ' due to: ' + str(err))

    print('\nWrote {} lat/lons to {}'.format(num_vals, out_fn))

    return

# ========================

def write_check_csv(form):
    """

    """
    func_name =  __prog__ + ' write_check_csv'

    # preliminary checks
    # ==================
    nc_dir = join(form.settings['weather_dir'],'CHESS_historic\\Monthly')
    nc_fnames = glob(nc_dir + '/*.nc')
    if len(nc_fnames) == 0:
        print('No NC files found')
        return

    out_dir = join(form.w_lbl_out_dir.text(), 'chess')
    if not isdir(out_dir):
        mkdir(out_dir)

    max_cells = int(form.w_max_cells.text())

    # =======================================
    meteo_df = _write_check_csv(out_dir, nc_fnames[0], form.settings, max_cells)

    return

def _make_lookup_csv_from_chess(out_dir, chess_nc_fn, meteo_lat_lon_osgb_fn, max_cells):
    """
    Create OSGB lookup table
    read a CHESS file and write a lookup table comprising Lat/Lons, Northings/Eastings and X/Y indices to a CSV file
    """
    nc_dset = Dataset(chess_nc_fn, mode='r')
    if 'x' not in nc_dset.variables:
        print(ERROR_STR + 'variable x must be in ' + chess_nc_fn)
        return None

    fnd_flag = False
    for varname in ['hurs', 'pr', 'precip']:
        if varname in nc_dset.variables:
            fnd_flag = True
            break

    if not fnd_flag:
        print(ERROR_STR + 'variable pr or precip must be in ' + chess_nc_fn)
        return None

    eastngs = [int32(x) for x in array(nc_dset.variables['x'])]
    nrthngs = [int32(y) for y in array(nc_dset.variables['y'])]
    num_total = len(eastngs) * len(nrthngs)

    lats = []
    y_coords = []
    y_indices = []

    lons = []
    x_coords = []
    x_indices = []

    ic = 0
    num_vals = 0
    num_nans = 0
    last_time = time()
    for yndx, nrthng in enumerate(nrthngs):

        for xndx, eastng in enumerate(eastngs):

            lat = nc_dset.variables['lat'][yndx][xndx]
            lon = nc_dset.variables['lon'][yndx][xndx]
            val = nc_dset.variables[varname][0][yndx][xndx]
            if val is MaskedConstant:
                num_nans += 1
            else:
                lats.append(lat)
                lons.append(lon)

                y_coords.append(int(nc_dset.variables['y'][yndx]))
                y_indices.append(yndx)

                x_coords.append(int(nc_dset.variables['x'][xndx]))
                x_indices.append(xndx)

                num_vals += 1

            ic += 1

            last_time = update_progress_post(last_time, num_nans, num_vals, num_total, round(lat,4), round(lon,4))

            if num_vals > max_cells:
                break

        if num_vals > max_cells:
            print('\nnumber of vals: {}\texceeds requested: {}'.format(num_vals, max_cells))
            break

    nc_dset.close()

    meteo_df = DataFrame()
    meteo_df['cell_lat'] = lats
    meteo_df['cell_lon'] = lons
    meteo_df['northing'] = y_coords
    meteo_df['yindx'] = y_indices
    meteo_df['easting'] = x_coords
    meteo_df['xindx'] = x_indices
    meteo_df = meteo_df.sort_values(by=['cell_lat', 'cell_lon'], ascending=[False, True])

    try:
        meteo_df.to_csv(meteo_lat_lon_osgb_fn, index=False, header=True)
    except PermissionError as err:
        print('Could not create ' + meteo_lat_lon_osgb_fn + ' due to: ' + str(err))

    print('\nWrote {} lat/lons to {}'.format(num_vals, meteo_lat_lon_osgb_fn))

    return

def write_lookup_csv(form):
    """

    """

    # preliminary checks
    # ==================
    nc_dir = form.w_lbl_ncdir.text()
    nc_fnames = glob(nc_dir + '/*.nc')
    if len(nc_fnames) == 0:
        print('No NC files found')
        return

    out_dir = join(form.w_lbl_out_dir.text(), 'chess')
    if not isdir(out_dir):
        mkdir(out_dir)

    max_cells = int(form.w_max_cells.text())

    # =======================================
    chess_nc_fn = nc_fnames[0]
    meteo_df = _make_lookup_csv_from_chess(out_dir, chess_nc_fn, form.settings['lat_lon_osgb'], max_cells)

    return
