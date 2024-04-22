#-------------------------------------------------------------------------------
# Name:        mscnfr_osgb_fns.py
# Purpose:     Diverese functions
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = ' mscnfr_osgb_fns.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import isdir, split, exists, join, isfile, splitext
from os import mkdir, remove
from glob import glob
import csv
import sys
from math import cos
from time import time

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

timeElapsed = 2

def write_mscnfr_out(pettmp, writers, num_time_steps, meteogrid_flag = True):
    """
    write each variable to a separate file
    """
    for var_name in pettmp.keys():

        # meteogrid is passed as a list comprising lat/long
        # =================================================
        if var_name == 'meteogrid':
            if meteogrid_flag:
                newlist = ['{:10d}'.format(val) for val in pettmp[var_name]]
                rec = ''.join(newlist)
                writers[var_name].writerow(list([rec]))
        else:
            # other metrics are passed as an ndarray which we convert to an integers after times by 10
            # ========================================================================================
            newlist = ['{:8d}'.format(int(10.0*val)) for val in pettmp[var_name]]
            for indx in range(0, num_time_steps, 12):
                rec = ''.join(newlist[indx:indx + 12])
                writers[var_name].writerow(list([rec]))

    return

def open_file_sets_osgb(var_names, out_folder, start_year = 1901, stop_year = 2019, grid_size = 1000,
        nrth_min = 0, nrth_max = 1056500, east_min = 0, east_max = 655500, out_suff = '.txt', remove_flag = True):

    if not isdir(out_folder):
        print(out_folder + ' does not exist - please reselect output folder')
        return None, None

    header_recs = []
    header_recs.append('GridSize    ' + str(round(grid_size,3)))
    header_recs.append('EastMin     ' + str(east_min))
    header_recs.append('EastMax     ' + str(east_max))
    header_recs.append('NorthMin    ' + str(nrth_min))
    header_recs.append('NorthMax    ' + str(nrth_max))
    header_recs.append('StartYear   ' + str(start_year))
    header_recs.append('StopYear    ' + str(stop_year))

    '''
    short_fnames = dict ( {'humidity': 'UKCP18_RHumidity','radiat_short': 'UKCP18_RadShortWaveNet', 'precip': 'precip',
                           'tempmin': 'Tmin', 'tempmax': 'Tmax', 'wind': 'Wind',
                            'cld': 'cloud','dtr': 'temprange', 'tmp': 'temperature', 'pre': 'precip', 'pet': 'pet'})
    '''
    miscan_fobjs = {}; writers = {}

    # for each file write header records
    # ==================================
    for var_name in var_names:
        file_name = join(out_folder, var_name + out_suff)
        if remove_flag:
            if exists(file_name):
                remove(file_name)

        miscan_fobjs[var_name] = open(file_name, 'w', newline='')
        writers[var_name] = csv.writer(miscan_fobjs[var_name])

        if var_name == 'meteogrid':
            hdr_recs = header_recs[0:5]
        else:
            hdr_recs = header_recs

        for header_rec in hdr_recs:
            writers[var_name].writerow(list([header_rec]))

    return miscan_fobjs, writers
