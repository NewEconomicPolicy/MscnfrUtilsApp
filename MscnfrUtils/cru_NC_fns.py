#-------------------------------------------------------------------------------
# Name:        cru_NC_fns.py
# Purpose:     Functions to create Miscanfor formatted metric files from CRU CSV files
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'cru_NC_fns.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import split, join, isdir
from os import mkdir
import sys
from time import time
from netCDF4 import Dataset
from numpy.ma import is_masked
from math import floor

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

from mscnfr_utils_fns import open_file_sets_cru, update_progress_post, write_mscnfr_out

timeElapsed = 5
resol = 0.5   # degrees
out_delim = '\t'
numSecsDay = 3600*24
NLONS = 720

def fetch_cru_data(var_names, nc_dsets, lat_indx, lon_indx):

    '''
    check each metric and if data is not present then return
    '''

    pettmp = {}

    for var_name in var_names:

        # check first 10 values
        # =====================
        data_flag = True
        for timindx in range(10):
            val = nc_dsets[var_name].variables[var_name][timindx,lat_indx,lon_indx]

            if is_masked(val):
                return None

        # collect readings for all time values
        pettmp[var_name] = nc_dsets[var_name].variables[var_name][:, lat_indx, lon_indx]

    return pettmp

def generate_cru_csv_files(form):
    '''

    '''
    func_name =  __prog__ + ' generate_csv_files'

    if len(form.settings['cru_fnames']) == 0:
        print('No CRU datasets to process')
        return

    print('In function ' + func_name)
    latitudeZeroIndex = -89.75

    # retrieve latitude ranges etc
    # ============================
    lat_min = float(form.w_min_lat.text())
    min_lat_indx = floor((lat_min - latitudeZeroIndex)/resol)
    lat_max = float(form.w_max_lat.text())
    max_lat_indx = floor((lat_max - latitudeZeroIndex)/resol)
    print('Latitude minimum and maximum indices: {} {} '.format(min_lat_indx, max_lat_indx))
    output_dir = join(form.w_lbl_out_dir.text(), 'cru')
    if not isdir(output_dir):
        mkdir(output_dir)

    # for user feedback
    nlats = max_lat_indx - min_lat_indx + 1
    num_total = nlats*NLONS
    max_cells = int(form.w_max_cells.text())

    # open all ncfiles
    # ================
    nc_dsets = {}
    var_names = []
    for long_fname in form.settings['cru_fnames']:
        dummy, fname = split(long_fname)
        fcomponents = fname.split('.')
        var_name = fcomponents[-3]
        var_names.append(var_name)
        nc_dsets[var_name] = Dataset(long_fname, mode='r')

        print('Have opened NetCDF file for metric: ' + var_name + ' long_name: ' + nc_dsets[var_name].variables[var_name].long_name)

    num_metrics = len(var_names)
    print('Will generate {} csv files consisting of metrics and a meteogrid file consisting of grid coordinates'\
                                                                                                .format(num_metrics))
    miscan_fobjs, writers = open_file_sets_cru(var_names + ['elevation','meteogrid'], output_dir, lat_min, lat_max)

    last_time = time()
    start_time = time()
    num_nodata = 0; num_with_data = 0
    num_time_steps = len(nc_dsets[var_name].variables['time'])
    print('Each grid cell for each metric should have {} time steps corresponding to {} years'\
                                                                        .format(num_time_steps, num_time_steps/12))
    # for each location, where there is data, build set of data
    # =========================================================
    for lon_indx in range(NLONS):
        lon = nc_dsets[var_name].variables['lon'][lon_indx]

        for lat_indx in range(min_lat_indx, max_lat_indx + 1):

            lat = nc_dsets[var_name].variables['lat'][lat_indx]
            last_time = update_progress_post(last_time, num_nodata, num_with_data, num_total, lat, lon)
            pettmp = fetch_cru_data(var_names, nc_dsets, lat_indx, lon_indx)

            # write data
            # ==========
            if pettmp == None:
                num_nodata += 1
            else:
                num_with_data += 1
                pettmp['meteogrid'] = list([lon, lat])
                write_mscnfr_out(pettmp, writers, num_time_steps)
                if num_with_data >= max_cells:
                    break

            if num_with_data >= max_cells:
                last_time = update_progress_post(last_time, start_time, num_nodata, num_with_data, num_total, lat, lon)
                break

    # close netCDF and csv files
    # ==========================
    for var_name in var_names:
        nc_dsets[var_name].close()
        miscan_fobjs[var_name].close()

    print('\nAll done...')

    return

def check_csv_file(form):
    '''
    check CSV file where each record comprises lat,lon, mu_gloobal
    '''
    csv_file = form.settings['csv_file']
    print('Reading ' + csv_file)
    fobj = open(csv_file)
    lines = fobj.readlines()
    fobj.close()

    # delete headers
    del lines[0]

    latitudes = set()
    longitudes = set()
    mu_globals = set()

    print()
    previous_time = time()
    for ic, line in enumerate(lines):
        rec = line.rstrip('\n').split(',')
        latitudes.add(float(rec[0]))
        longitudes.add(float(rec[1]))
        mu_globals.add(int(rec[2]))
        this_time = time()
        if this_time - previous_time > timeElapsed:
            sys.stdout.write('\rNumber of lines read: ' + format_string("%d", ic, grouping=True))
            previous_time = this_time

    min_lat = min(latitudes)
    max_lat = max(latitudes)
    min_lon = min(longitudes)
    max_lon = max(longitudes)

    print('\nRead {} lines\tMin/Max latitude: {} {}\tMin/Max longitude: {} {}\tnumber of unique mu globals: {}'
          .format(len(lines), min_lat, max_lat, min_lon, max_lon, len(mu_globals)))

    return
