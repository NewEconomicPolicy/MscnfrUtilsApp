#-------------------------------------------------------------------------------
# Name:        hadley_NC_fns.py
# Purpose:     Functions to create Miscanfor formatted metric files from Hadley netCDF climate files
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'hadley_NC_fns.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import isdir, exists, join
from os import remove, mkdir
from time import time
from netCDF4 import Dataset
from numpy.ma import is_masked
from math import floor
import csv
from copy import copy

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

from mscnfr_utils_fns import update_progress_post

ERROR_STR = '*** Error *** '

timeElapsed = 2
NLONS = 720
RESOL = 0.5   # degrees
OUT_DELIM = '\t'
numSecsDay = 3600*24

SHORT_FNAMES = dict({'rlds': 'radiation_long', 'tasmin': 'tempmin', 'tasmax': 'tempmax', 'pr': 'precip',
                     'rsds': 'radiation_short'})
#  SHORT_FNAMES = dict ( {'monthly_wind':'wind_speed'})

def generate_hadley_csv_files(form):
    """

    """
    func_name =  __prog__ + ' generate_hadley_csv_files'

    print('In function ' + func_name)
    if form.settings['hadley_year_sets'] is None:
        print('No Hadley datasets to process')
        return

    latitudeZeroIndex = 89.75
    hadley_sets_mapped = form.settings['hadley_sets_mapped']
    hadley_year_sets   = form.settings['hadley_year_sets']
    first_year_set = hadley_year_sets[0]

    metric_names = SHORT_FNAMES.keys()
    metric_names_out = list(metric_names)
    metric_names_out.append('tasmean')
    metric_names_out.append('tasrange')

    # retrieve latitude ranges etc
    # ============================
    lat_min = float(form.w_min_lat.text())
    max_lat_indx = floor((latitudeZeroIndex - lat_min)/RESOL)
    lat_max = float(form.w_max_lat.text())
    min_lat_indx = floor((latitudeZeroIndex - lat_max)/RESOL)
    print('Latitude minimum and maximum indices: {} {} '.format(min_lat_indx, max_lat_indx))
    output_dir = join(form.w_lbl_out_dir.text(), 'hadley')
    if not isdir(output_dir):
        mkdir(output_dir)

    # for user feedback
    # =================
    nlats = max_lat_indx - min_lat_indx + 1
    num_total = nlats*NLONS
    max_cells = int(form.w_max_cells.text())

    print('Will generate {} csv files consisting of {} and a meteogrid file of grid coordinates'\
                                                .format(len(metric_names), ' '.join(sorted(SHORT_FNAMES.values()))))
    miscan_fobjs, writers = _open_file_sets_hadley(metric_names_out, output_dir, lat_min, lat_max,
                                                   start_year = form.settings['hadley_start_year'],
                                                   stop_year = form.settings['hadley_end_year'])
    if writers == None:
        return

    # open all ten nc year sets for each of the five metrics
    # ======================================================
    nc_dsets = {}
    nc_variables = {}
    for metric in metric_names:
        nc_dsets[metric] = {}
        nc_variables[metric] = {}
        mapped_name = hadley_sets_mapped[metric]

        for num_yr_set, year_set in enumerate(hadley_year_sets):
            try:
                nc_dsets[metric][year_set] = Dataset(form.settings['hadley_fnames'][metric][num_yr_set], mode='r')
            except IndexError as err:
                print(ERROR_STR + 'metric: {}\t# year_set: {}\t'.format(metric, num_yr_set) + str(err))
                return

            nc_variables[metric][year_set] = nc_dsets[metric][year_set].variables[mapped_name]

        mess = 'Have opened NetCDF files for metric: ' + metric
        if 'long_name' in nc_variables[metric][year_set].ncattrs():
            print(mess + ' long_name: ' + nc_variables[metric][year_set].getncattr('long_name'))
        else:
            print(mess)

    # lats and longs are same for all NC files
    # ========================================
    latitudes = nc_dsets[metric][year_set].variables['lat'][:]
    longitudes = nc_dsets[metric][year_set].variables['lon'][:]

    # collect days in each month - since time variable starts with January assume 31 for first value
    # ==============================================================================================
    ave_days_per_mnth = list([])
    for year_set in hadley_year_sets:
        time_var = nc_dsets[metric][year_set].variables['time'][:]
        prev_day = time_var[0]
        ave_days_per_mnth.append(time_var[1] - prev_day)
        for day in time_var[1:]:
            ave_days_per_mnth.append(day - prev_day)
            prev_day = day

    num_time_steps = len(ave_days_per_mnth)
    print('Each grid cell has {} monthly time steps corresponding to {} years and {} year sets'\
                                               .format(num_time_steps, int(num_time_steps/12), len(hadley_year_sets)))

    # main loop: step through 720 longitudes from -179.75 to 179.75 degrees
    # for each lat/long location, where there is data, build set of data
    # ==================================================================
    last_time = time()
    num_nodata = 0
    num_with_data = 0
    pettmp = {}

    for lon_indx in range(NLONS):
        lon = longitudes[lon_indx]

        # collect readings for all time and latitude values for each metric and year set
        # ==============================================================================
        for metric in metric_names:
            mapped_name = hadley_sets_mapped[metric]
            pettmp[metric] = {}
            for year_set in hadley_year_sets:
                pettmp[metric][year_set] = nc_dsets[metric][year_set].variables[mapped_name][:, :, lon_indx]

        # step through latitudes from the south proceding northwards
        # ==========================================================
        for lat_indx in range(max_lat_indx, min_lat_indx - 1, -1):
            lat = latitudes[lat_indx]

            # check that there is data for this grid point
            # ============================================
            write_records_flag = True
            for metric in metric_names:
                val = pettmp[metric][first_year_set][0][lat_indx]
                if is_masked(val):
                    num_nodata += 1
                    write_records_flag = False
                    break

            # append data for subsequent year sets and write
            # ==============================================
            if write_records_flag:
                num_with_data += 1
                meteogrid = list([lon, lat])
                out_pettmp = {}
                for metric in metric_names:
                    out_pettmp[metric] = []
                    for year_set in hadley_year_sets:
                        out_pettmp[metric] += list(pettmp[metric][year_set][:, lat_indx])

                _write_out_hadley(writers, out_pettmp, metric_names_out, ave_days_per_mnth, meteogrid, num_time_steps)

                last_time = update_progress_post(last_time, num_nodata, num_with_data, num_total, lat, lon)

            # break out of latitudes loop
            # ===========================
            if num_with_data >= max_cells:
                break

        # break out of longitudes loop
        # ============================
        if num_with_data >= max_cells:
            break

    # close all files
    # ===============
    for metric in metric_names:
        miscan_fobjs[metric].close()

        for year_set in form.settings['hadley_year_sets']:
            nc_dsets[metric][year_set].close()

    print('\nHave closed all CSV and netCDF files - finished...')

    return

def _write_out_hadley(writers, out_pettmp, metric_names, ave_days_per_mnth, meteogrid, num_time_steps):
    """
    write each variable to a separate file
    """
    
    # meteogrid is passed as a list comprising lat/long - only write to meteogrid.txt once
    # ====================================================================================
    newlist = ['{:10.4f}'.format(val) for val in meteogrid]
    rec = ''.join(newlist)
    writers['meteogrid'].writerow(list([rec]))

    # stanza to create mean temperature
    if 'tasmean' in out_pettmp.keys():
        out_pettmp['tasmean'] = []
        out_pettmp['tasrange'] = []
        for tmax, tmin in zip(out_pettmp['tasmax'],out_pettmp['tasmin']):
            out_pettmp['tasmean'].append((tmax+tmin)/2.0)
            out_pettmp['tasrange'].append(tmax-tmin)

    for var_name in metric_names:

        if var_name not in out_pettmp.keys():
            continue
        out_rec = out_pettmp[var_name]

        # metrics are passed as a list of real values which we convert to integers after times by 10
        # ==========================================================================================
        if var_name == 'pr':
            ic = 0
            for val, ave_days in zip(out_rec, ave_days_per_mnth):
                out_rec[ic] = val*numSecsDay*ave_days
                ic += 1

        newlist = ['{:8d}'.format(int(10.0*val)) for val in out_rec]
        for indx in range(0, num_time_steps, 12):
            rec = ''.join(newlist[indx:indx + 12])
            writers[var_name].writerow(list([rec]))

    return

def _open_file_sets_hadley(var_names, out_folder, lat_min, lat_max,
                           start_year = 1901, stop_year = 2016, out_suff = '.txt', remove_flag = True):
    """

    """
    if not isdir(out_folder):
        print(out_folder + ' does not exist - please reselect output folder')
        return None, None

    # common header except for meteogrid
    # ==================================
    header_recs = []
    header_recs.append('GridSize    0.5000000')
    header_recs.append('LongMin     -180.0000')
    header_recs.append('LongMax      180.0000 ')
    header_recs.append('LatMin      ' + str(lat_min))
    header_recs.append('LatMax      ' + str(lat_max))
    header_recs.append('StartYear   ' + str(start_year))
    header_recs.append('StopYear    ' + str(stop_year))

    miscan_fobjs = {}; writers = {}

    short_fnames_out = copy(SHORT_FNAMES)
    short_fnames_out['tasmean'] = 'tempmean'
    short_fnames_out['tasrange'] = 'temprange'

    # for each file write header records
    # ==================================
    for var_name in var_names:
        file_name = join(out_folder, short_fnames_out[var_name] + out_suff)
        if remove_flag:
            if exists(file_name):
                try:
                    remove(file_name)
                except (PermissionError) as err:
                    print('Could not remove file ' + file_name)
                    return None, None

        miscan_fobjs[var_name] = open(file_name, 'w', newline='')
        # writers[var_name] = csv.writer(miscan_fobjs[var_name], delimiter = out_delim)
        writers[var_name] = csv.writer(miscan_fobjs[var_name])
        for header_rec in header_recs:
            writers[var_name].writerow(list([header_rec]))

    # add two less lines for the grid coordinates file
    meteogrid = 'meteogrid'
    file_name = join(out_folder, meteogrid + out_suff)
    if remove_flag:
        if exists(file_name):
            try:
                remove(file_name)
            except (PermissionError) as err:
                print('Could not remove file ' + file_name)
                return None, None

    miscan_fobjs[meteogrid] = open(file_name, 'w', newline='')
    # writers[meteogrid] = csv.writer(miscan_fobjs[meteogrid], delimiter = out_delim)
    writers[meteogrid] = csv.writer(miscan_fobjs[meteogrid])
    for header_rec in header_recs[0:5]:
        writers[meteogrid].writerow(list([header_rec]))

    return miscan_fobjs, writers
