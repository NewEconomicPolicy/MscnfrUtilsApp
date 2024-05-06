#-------------------------------------------------------------------------------
# Name:        mscnfr_utils_fns.py
# Purpose:     Diverese functions
# Author:      Mike Martin
# Created:     22/03/2020
# Description:
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'mscnfr_utils_fns.py'
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
OUT_DELIM = ','
HADLEY_SETS_MAPPED = {'pr': 'prAdjust', 'rlds': 'rldsAdjust', 'rsds': 'rsdsAdjust', 'tasmin': 'tasmin',
                      'tasmax': 'tasmax', 'monthly_wind': 'windAdjust'}
WARN_STR = '*** Warning *** '
'''
from geopandas import GeoDataFrame
from numpy import array
from scipy.spatial import cKDTree
from shapely.geometry import Point

def ckdnearest(gdA, gdB):
    nA = array(list(zip(gdA.geometry.x, gdA.geometry.y)) )
    nB = array(list(zip(gdB.geometry.x, gdB.geometry.y)) )
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=1)
    gdf = concat(
        [gdA.reset_index(drop=True), gdB.loc[idx, gdB.columns != 'geometry'].reset_index(drop=True),
                                                                        Series(dist, name='dist')], axis=1)
    return gdf

'''
def check_coord_bbox(lat, lon, bbox):
    """
    this function returns true or false depending on whether lat/lon lies within specicified bounding box
    """
    func_name =  __prog__ + ' check_coord_bbox'

    lon_ll, lat_ll, lon_ur, lat_ur = bbox
    if (lat_ur >= lat and lat >= lat_ll) and (lon_ur >= lon and lon >= lon_ll):
        return True
    else:
        return False

def get_block_ref(date_rec):

    yr, mnth, day = [val for val in date_rec.split('-')]
    block_ref = yr + '_' + mnth

    return block_ref

def find_grid_resolution(grid_cell_km, latitude):
    """
    determnie approximate cell size in degrees for a given latitude
    """
    radius_earth = 6371.0 # km average earth radius at equator
    pi = 3.14159265
    deg2rad = pi/180.0

    km_per_deg = 2*pi*cos(latitude*deg2rad)*radius_earth/360
    cell_resol = grid_cell_km/km_per_deg

    return cell_resol

def write_mscnfr_out(pettmp, writers, num_time_steps, meteogrid_flag = True):
    """
    write each variable to a separate file
    """
    for var_name in pettmp.keys():

        # meteogrid is passed as a list comprising lat/long
        # =================================================
        if var_name == 'meteogrid':
            if meteogrid_flag:
                newlist = ['{:10.4f}'.format(val) for val in pettmp[var_name]]
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

def open_file_sets_cru(var_names, out_folder, lat_min, lat_max, lon_min = -180.0, lon_max = 180, grid_size = 0.5,
            start_year = 1901, stop_year = 2019, out_suff = '.txt', remove_flag = True):

    if not isdir(out_folder):
        print(out_folder + ' does not exist - please reselect output folder')
        return None, None

    header_recs = []
    header_recs.append('GridSize    ' + str(round(grid_size,3)))
    header_recs.append('LongMin     ' + str(lon_min))
    header_recs.append('LongMax     ' + str(lon_max))
    header_recs.append('LatMin      ' + str(lat_min))
    header_recs.append('LatMax      ' + str(lat_max))
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

def identify_ukcp18_dirs(form, prefix = 'Perturbations: '):
    """
    metric_dirs are typically:
        precip, Tmax, Tmin, UKCP18_RHumidity, UKCP18_RadShortWaveNet and Wind
    Typical file:
        subset_2020-02-14T16-13-39_HadREM3-GA705-r001i1p00000.csv
    each monthly time step is defined is a block of 65 columns (Eastings) by 112 rows (Northings)
    """
    ukcp18_dir = form.w_lbl_ukcp18.text()
    if not isdir(ukcp18_dir):
        return prefix

    metric_dirs = glob(ukcp18_dir + '/*')
    perturbs = []

    # since all paths have similar sets of files only the first one is used
    # =====================================================================
    for dirname in metric_dirs:
        perturbs = []
        pstr = ''
        metric_fnames  = glob(dirname + '/subset*.csv')
        for fname in metric_fnames:
            dummy, short_fname = split(fname)
            root_name, dummy = splitext(short_fname)
            ensemble = root_name.split('-')[-1]
            perturb = ensemble.lstrip('r001i1')
            perturbs.append(perturb)
            pstr += perturb + ', '
        break

    form.settings['perturbs'] = perturbs
    return prefix + pstr.rstrip(', ')

def reformat_csv_files(form):
    """

    """
    func_name =  __prog__ + ' reformat_csv_files'
    print('In function ' + func_name)

    max_num_cells = int(form.w_max_cells.text())

    num_created = 0
    root_dir = join(form.settings['weather_dir'], 'Hwsd_CSVs')
    dump_dir = form.settings['mscnfr_dump_dir']
    country_list = list(['North_Korea','South_Korea','Taiwan','Japan','China'])

    for country in country_list:

        out_fname = join(dump_dir,country + '.csv')
        if isfile(out_fname):
            remove(out_fname)

        if not isfile(out_fname):
            print(WARN_STR + out_fname + ' does not exist')
            continue

        fname_obj = open(out_fname, 'w', newline='')
        fname_writer = csv.writer(fname_obj, delimiter = OUT_DELIM)
        print('\nOpened ' + out_fname)
        nwrites = 0

        dir_name = join(root_dir,country)
        csv_files = glob( dir_name + '\\*.csv')
        for csv_file in csv_files:
            nreads = 0
            fninp_obj = open(csv_file, 'r', newline='')
            fninp_reader = csv.reader(fninp_obj, delimiter = OUT_DELIM)

            for row_in in fninp_reader:
                # reorganise row
                row_out = list([ row_in[3], row_in[4], row_in[2]])
                fname_writer.writerow( row_out )
                nreads += 1
                nwrites += 1
                if nreads > max_num_cells:
                    break

            print('Closed ' + csv_file + ' after {} reads'.format(nreads))
            fninp_obj.close()

        print('Closed ' + out_fname + ' after {} writes'.format(nwrites))
        fname_obj.close()
        num_created += 1

    return num_created

def identify_nc_files(form, nc_dir):
    """
    conversion_flag can have three states: None, Hadley and Cru
    CRU files are located in CRU 4_04:
                    cru_ts4.04.1901.2019.dtr.dat.nc   diurnal temperature range
    Hadley files are located in CRU 4_04:
                    pr_bced_1960_1999_hadgem2-es_rcp2p6_2071-2080.nc   Bias-Corrected Precipitation
    CRU and Hadley files use Lat/Longs
    """
    if not isdir(nc_dir):
        print('NetCDF directory ' + nc_dir + ' does not exist')
        form.settings['hadley_fnames'] = None
        form.settings['cru_fnames'] = []
        form.settings['hadley_year_sets'] = None
        return


    form.settings['hadley_sets_mapped'] = HADLEY_SETS_MAPPED
    hadley_metrics = HADLEY_SETS_MAPPED.keys()
    form.settings['hadley_metrics'] = hadley_metrics

    # stanza to filter out cloud cover and PET
    cru_fnames = []
    fnames = glob(nc_dir + '/cru*.nc')
    for long_fname in fnames:
        dummy, fname = split(long_fname)
        fcomponents = fname.split('.')
        var_name = fcomponents[-3]
        cru_fnames.append(long_fname)       # was: if var_name == 'cld' or var_name == 'pet':

    nfiles_convert = len(cru_fnames)

    # ===================
    year_sets = []
    hadley_fnames = {}
    if nfiles_convert > 0:
        conversion_flag = 'CRU'
        suffix = ' NetCDF files'
    else:
        # create dictionary of nc files with each element consisting
        # of an ordered set of nc files for that metric
        # =============================================
        for metric in hadley_metrics:
            hadley_fnames[metric] = glob(nc_dir + '/' + metric + '*hadgem*.nc')
        nfiles_convert = len(hadley_fnames['pr'])
        if nfiles_convert == 0:
            print('No netCDF files to convert in ' + nc_dir)
            form.w_combo18.setCurrentText('None')
            return
        else:
            conversion_flag = 'Hadley'
            suffix = ' NetCDF file sets'
            # extract year sets
            for fname in hadley_fnames['pr']:
                root_name, extens = splitext(fname)
                year_sets.append(root_name.split('_')[-1])

            start_year, end_year = year_sets[0].split('-')
            form.settings['hadley_start_year'] = int(start_year)
            start_year, end_year = year_sets[-1].split('-')
            form.settings['hadley_end_year'] = int(end_year)

    form.w_lbl_num_ncs.setText('Will convert '+ str(nfiles_convert) + ' ' + conversion_flag + suffix)
    form.settings['hadley_fnames'] = hadley_fnames
    form.settings['cru_fnames'] = cru_fnames
    form.settings['hadley_year_sets'] = year_sets
    form.w_combo18.setCurrentText(conversion_flag)

def update_progress_hwsd(last_time, ic, nrecs):
    """
    Update progress bar
    """
    this_time = time()
    if (this_time - last_time) > timeElapsed:
        sys.stdout.flush()
        percent_done = round(100.0*(1.0 - (nrecs - ic)/nrecs), 2)
        out_str = '\rProcessed: {:=6.2f}% of HWSD csv file'.format(percent_done)
        sys.stdout.write(out_str)
        return this_time

    return last_time

def update_progress_chess(last_time, num_set, num_atoms, nbad_cells, num_total, lat, lon):
    """
    Update progress bar
    """
    this_time = time()
    if (this_time - last_time) > timeElapsed:
        sys.stdout.flush()

        num_set_str = format_string("%d", num_set, grouping=True)
        num_bad_str = format_string("%d", nbad_cells, grouping=True)
        num_atoms_str = format_string("%d", num_atoms, grouping=True)
        percent_done = round(100.0*(1.0 - (num_total - num_set - nbad_cells)/num_total), 4)
        out_str = '\rN good cells: {:10s}\trejected: {:10s}\tN atoms: {:10s}\t{:=6.2f}%\tlat/lon: {:=7.3f} {:=7.3f}'\
                                .format(num_set_str, num_bad_str, num_atoms_str, percent_done, lat, lon)
        sys.stdout.write(out_str)
        return this_time

    return last_time

def update_progress_post(last_time, num_nodata, num_with_data, num_total, lat, lon, num_set = 1):
    """
    Update progress bar
    """
    this_time = time()
    if (this_time - last_time) > timeElapsed:
        sys.stdout.flush()
        num_with_data_str = format_string("%d", num_with_data, grouping=True)
        num_nodata_str = format_string("%d", num_nodata, grouping=True)
        percent_done = round(100.0*(1.0 - (num_total - num_nodata - num_with_data)/num_total),2)
        out_str = '\rSet: {:=2d} Cells with data: {:10s} none: {:10s} {:=5.1f}%\tlat/lon: {:=7.3f} {:=7.3f}'\
                                .format(int(num_set), num_with_data_str, num_nodata_str, percent_done, lat, lon)
        sys.stdout.write(out_str)
        return this_time

    return last_time

def update_progress_ukcp18(last_time, perturb, metric, nblock_cntr, ntotal_blocks):

    """Update progress bar."""
    this_time = time()
    if (this_time - last_time) > 3.0:
        sys.stdout.flush()
        nblock_cntr_str = format_string("%d", nblock_cntr, grouping=True)
        percent_done = round(100.0*(1.0 - (ntotal_blocks - nblock_cntr)/ntotal_blocks),1)
        sys.stdout.write('\rPerturbation: {}\tmetric: {:12s}\tblock counter: {}\tcomplete: {:=5.1f}%'
                .format(perturb, metric, nblock_cntr_str, percent_done))
        return this_time

    return last_time

def update_progress_merge(last_time, num_outside, num_near, num_match, num_written, num_recs, seq_no):

    """Update progress bar."""
    this_time = time()
    if (this_time - last_time) > 3.0:
        sys.stdout.flush()
        num_written_str = format_string("%d", num_written, grouping=True)
        num_out_str     = format_string("%d", num_outside, grouping=True)
        num_near_str    = format_string("%d", num_near, grouping=True)
        num_match_str   = format_string("%d", num_match, grouping=True)
        percent_done = round(100.0*(1.0 - (num_recs - num_outside - num_near - num_match)/num_recs),1)
        sys.stdout.write('\rSet: {:=2} Cells written: {}\toutside: {}\tnear: {}\tmatched: {} {:=5.1f}%'
                .format(seq_no, num_written_str, num_out_str, num_near_str, num_match_str, percent_done))
        return this_time

    return last_time