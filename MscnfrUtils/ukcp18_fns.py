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
from os import mkdir, remove
from time import time
from calendar import monthrange
from glob import glob

from csv import writer
from pandas import read_csv, DataFrame
from netCDF4 import Dataset
from numpy import ma, float64

from locale import setlocale, format_string, LC_ALL
setlocale(LC_ALL, '')

from mscnfr_utils_fns import open_file_sets_cru, update_progress_ukcp18, get_block_ref, update_progress_post,\
                                                                                            find_grid_resolution
HWSD_FNAME = 'E:\\GlobalEcosseData\\Hwsd_CSVs\\UK\\UK_hwsd.csv'

NGRANULARITY = 120
NLINES_PER_BLOCK = 114
MAX_PERTURBS  = 1000
MAX_PERTURBS  = 1
MAX_METRICS   = 1000
READ_FLAG = True

METRIC_LOOKUP = {'precip': 'precip', 'Tmax': 'tempmax', 'Tmin': 'tempmin',
                        'UKCP18_RadShortWaveNet': 'radiat_short', 'UKCP18_RHumidity' : 'humidity', 'Wind': 'wind'}
METRICS = list(METRIC_LOOKUP.values())

def _write_ukcp18_elevations(out_dir, meteo_df):
    """
    write miscanfor file sets
    """
    lat_min = meteo_df['cell_lat'].min()
    lat_max = meteo_df['cell_lat'].max()
    lon_min = meteo_df['cell_lon'].min()
    lon_max = meteo_df['cell_lon'].max()

    header_recs = []
    header_recs.append('GridSize    0.5000000')
    header_recs.append('LongMin     ' + str(lon_min))
    header_recs.append('LongMax     ' + str(lon_max))
    header_recs.append('LatMin      ' + str(lat_min))
    header_recs.append('LatMax      ' + str(lat_max))

    file_name = join(out_dir, 'elevation.txt')
    if exists(file_name):
        remove(file_name)

    fobj = open(file_name, 'w', newline='')
    elev_writer = writer(fobj)
    for header_rec in header_recs:
        elev_writer.writerow(list([header_rec]))

    for lat, lon, elevation in zip(meteo_df['cell_lat'].values, meteo_df['cell_lon'].values, meteo_df['elevations'].values):
        rec = '{:10.4f}{:10.4f}{:10.1f}'.format(lat, lon, elevation)
        elev_writer.writerow(list([rec]))

    fobj.close()

    print('wrote ' + file_name)

    return

def _fetch_ukcp18_elevations(uk_coords_dir, out_dir, meteo_df):
    """
    create Miscanfor compatible file of elevations
    """
    func_name =  __prog__ + ' _write_ukcp18_elevations'

    nc_file = join(uk_coords_dir, 'ELEVN_orog_land-rcm_uk_12km_osgb.nc')
    if not isfile(nc_file):
        print('Elevation NC file: ' + nc_file + ' does not exist')
        return

    nc_dset = Dataset(nc_file, mode='r')
    elev_northings = list(nc_dset.variables['projection_y_coordinate'][:])
    elev_eastings  = list(nc_dset.variables['projection_x_coordinate'][:])

    bad_locat = 0
    elevations = []
    for easting, northing in zip(meteo_df['easting'], meteo_df['northing']):
        bad_locat_flag = False
        try:
            xindx = elev_eastings.index(easting)
        except ValueError as e:
            print('Easting: {} not found in evelvations NC file')
            bad_locat_flag = True
        try:
            yindx = elev_northings.index(northing)
        except ValueError as e:
            print('Northing: {} not found in evelvations NC file')
            bad_locat_flag = True

        if bad_locat_flag:
            bad_locat += 1
        if bad_locat > 10:
            print('Too many bad locations - aborting')
            break

        elevation = nc_dset.variables['surface_altitude'][yindx, xindx]
        if ma.is_masked(elevation):
            elevation = float64(0.0)
        elevations.append(elevation)

    nc_dset.close()

    return elevations

def _filter_block(meteo_df, block_df, block_ref, metric):
    """

    """
    year, mnth = block_ref.split('_')
    ndays = monthrange(int(year), int(mnth))[1]

    val_list = []
    for easting, northing in zip(meteo_df['easting'], meteo_df['northing']):
        val = block_df[str(easting)][str(northing)]
        if metric == 'precip':
            val_list.append(ndays*val)
        else:
            val_list.append(val)

    return val_list

def generate_ukcp18_csv_files(form):
    """
    write each variable to a separate file
    """
    func_name =  __prog__ + ' generate_ukcp18_csv_files'

    # preliminary checks
    # ==================
    ukcp18_dir = form.w_lbl_ukcp18.text()
    if not isdir(ukcp18_dir) :
        print('UK CP18 path: ' + ukcp18_dir + ' must exist')
        return

    metric_dirs = glob(ukcp18_dir + '/*')
    perturbs = form.settings['perturbs']
    if len(perturbs) == 0 or perturbs is None:
        print('UK CP18 conversion must have perturbations')
        return

    out_dir = join(form.w_lbl_out_dir.text(), 'ukcp18')
    if not isdir(out_dir):
        mkdir(out_dir)

    max_blocks = int(form.w_max_cells.text())

    # essential preliminaries
    # =======================
    dirname = metric_dirs[0]
    root_dir = split(split(dirname)[0])[0]
    uk_coords_dir = join(root_dir, 'UK_coords')
    cnvrtd_coords_fname = join(uk_coords_dir, 'UK_conversion_coords.txt')
    if not isfile(cnvrtd_coords_fname):
        print('File: ' + cnvrtd_coords_fname + ' must exist - please check')
        return

    meteo_fname = join(uk_coords_dir,'meteogrid_extend.txt')
    if not isfile(meteo_fname):
        print('File: ' + meteo_fname + ' must exist - please check')
        return

    perturb = perturbs[0]
    this_block, nblocks, start_year, end_year, northing_keys, easting_keys = \
                                                                        _fetch_nrthng_and_eastng_keys(dirname, perturb)

    meteo_df = _generate_meteogrid(HWSD_FNAME, cnvrtd_coords_fname, meteo_fname, northing_keys, easting_keys)
    meteo_df['elevation'] = _fetch_ukcp18_elevations(uk_coords_dir, out_dir, meteo_df)

    # set limits
    # ==========
    nperturbs = min(MAX_PERTURBS, len(perturbs))
    nmetrics  = min(MAX_METRICS, len(metric_dirs))
    ntotal_blocks = min(max_blocks, nblocks)*nperturbs*nmetrics

    # each block is a timestep
    # ========================
    last_time  = time()
    result_dfs = {}
    nblock_cntr = 0
    for ipert, perturb in enumerate(perturbs):
        print('\nProcessing {} blocks from perturbation: {}'.format(nblocks, perturb))

        for imetr, dirname in enumerate(metric_dirs):

            dummy, metric_fold = split(dirname)
            metric = METRIC_LOOKUP[metric_fold]
            '''
            if metric != 'precip':
                continue
            '''
            fname_csv = glob(dirname + '\\*' + perturb + '.csv')[0]
            dummy, short_fname = split(fname_csv)

            # retrieve entire file content
            # ============================
            with open(fname_csv, 'r') as fobj:
                data = fobj.readlines()

            # first 15 lines are header records
            # =================================
            header_block = data[0:15]

            result_dfs[metric] = DataFrame()
            result_dfs[metric]['cell_lat'] = meteo_df['cell_lat']
            result_dfs[metric]['cell_lon'] = meteo_df['cell_lon']

            nlines = len(data)
            hblocks = (nlines - 15)/NLINES_PER_BLOCK
            nblocks = int(hblocks)
            if nblocks - hblocks != 0.0:
                print('*** Warning *** check length of file ' + short_fname + '\tfile length: {}, should be: 136815'.
                                                                                                        format(nlines))
            # read each block
            # ===============
            indx1 = 15
            indx2 = indx1 + NLINES_PER_BLOCK
            for nblock in range(nblocks):
                nblock_cntr += 1
                new_block = data[indx1:indx2]
                block_ref = get_block_ref(new_block[0])
                del new_block[0:2] # discard date and line of eastings

                this_block = [val.strip() for val in new_block]     # remove \n
                del new_block

                # use row of eastings an index
                # ============================
                recs = []
                for line in this_block:
                    rec = [float(sval) for sval in line.split(',')]
                    northing = str(int(rec[0]))
                    recs.append([northing] + rec[1:])

                block_df = DataFrame(recs, columns = easting_keys)
                block_df.set_index('northing', inplace = True)
                block_df = block_df.copy()  # attempt to address PerformanceWarning: DataFrame is highly fragmented

                result_dfs[metric][block_ref] = _filter_block(meteo_df, block_df, block_ref, metric)

                indx1 += NLINES_PER_BLOCK
                indx2 = indx1 + NLINES_PER_BLOCK

                last_time = update_progress_ukcp18(last_time, perturb, metric, nblock_cntr, ntotal_blocks)
                if nblock + 1 >= max_blocks:
                    print('\nExtraction complete after {} blocks read'.format(nblock))
                    break

            if imetr + 1 >= MAX_METRICS:
                break

        # write miscanfor file sets
        # =========================
        _write_ukcp18_results(out_dir, perturb, meteo_df, result_dfs, start_year, end_year)

        if ipert + 1 >= MAX_PERTURBS :
            break

    print('\nCreation of Miscanfor datasets complete after processing {} perturbations and {} metrics'
                                                                                        .format(ipert + 1, imetr + 1))
    return

def _write_ukcp18_results(out_dir, perturb, meteo_df, result_dfs, start_year, end_year):
    """
    write miscanfor file sets
    """
    lat_min = meteo_df['cell_lat'].min()
    lat_max = meteo_df['cell_lat'].max()
    lon_min = meteo_df['cell_lon'].min()
    lon_max = meteo_df['cell_lon'].max()
    lat_mid = lat_min + (lat_max - lat_min)/2
    grid_resol = find_grid_resolution(12, lat_mid)   # cell size in kms and latitude UK mid-point

    out_dir_pert = join(out_dir, perturb)
    if not isdir(out_dir_pert):
        mkdir(out_dir_pert)

    miscan_fobjs, writers = open_file_sets_cru(METRICS + ['elevation','meteogrid'], out_dir_pert,
                                                lat_min, lat_max, lon_min, lon_max, grid_resol, start_year, end_year)
    # meteogrid and elevation
    # =======================
    for lat, lon, elevation in zip(meteo_df['cell_lat'], meteo_df['cell_lon'], meteo_df['elevation']):

        writers['meteogrid'].writerow(list(['{:10.4f}{:10.4f}'.format(lon, lat)]))
        writers['elevation'].writerow(list(['{:8d}'.format(int(10*elevation))]))

    # all the rest
    # ============
    for metric in METRICS:
        if metric in result_dfs:
            for rec_result in result_dfs[metric].values:
                nblocks = len(rec_result) - 2
                newlist = ['{:8d}'.format(int(round(10.0*val))) for val in rec_result[2:]]
                for indx in range(0, nblocks, 12):
                    rec = ''.join(newlist[indx:indx + 12])
                    writers[metric].writerow(list([rec]))

    # close all file objects
    # ======================
    for var in miscan_fobjs:
        miscan_fobjs[var].close()

    return

def _fetch_nrthng_and_eastng_keys(dirname, perturb):
    """

    """
    fname_csv = glob(dirname + '\\*' + perturb + '.csv')[0]
    dummy, short_fname = split(fname_csv)

    # retrieve entire file contents
    # =============================
    with open(fname_csv, 'r') as fobj:
        data = fobj.readlines()

    # first 15 lines comprise header records
    header_block = data[0:15]

    nlines = len(data)
    hblocks = (nlines - 15)/NLINES_PER_BLOCK
    nblocks = int(hblocks)
    if nblocks - hblocks != 0.0:
        print('*** Warning *** check length of file ' + short_fname + '\tfile length: {}, should be: 136815'.
                                                                                                format(nlines))

    dummy, metric = split(dirname)
    print('Processing {} blocks from perturbation: {}\tCSV file: {}\t\tmetric: {}'
                                                                    .format(nblocks, perturb, short_fname, metric))

    indx1 = 15; indx2 = indx1 + NLINES_PER_BLOCK

    new_block = data[indx1:indx2]
    this_block = [val.strip() for val in new_block]     # remove \n
    del new_block

    this_date = this_block[0]   # read date then delete line
    del this_block[0]

    start_year, start_month, start_day = [int(val) for val in this_date.split('-')]
    if start_month == 12:
        start_year += 1
    end_year = int(start_year + nblocks/12 -1)

    # create headers
    # ==============
    hdr_str = this_block[0]
    hdrs = hdr_str.split(',')[1:]
    eastings = [str(int(float(sval))) for sval in hdrs]	# crude but only needs to be done once
    easting_keys = ['northing'] + eastings
    del this_block[0] # delete eastings

    # convert first field of rest of lines
    # ====================================
    northing_keys = []
    for line in this_block:
        northing = line.split(',')[0]
        northing_keys.append(northing)

    return this_block, nblocks, start_year, end_year, northing_keys, easting_keys

def _generate_meteogrid(hwsd_fname, cnvrtd_coords_fname, meteo_fname, northing_keys, easting_keys,
                        read_flag = READ_FLAG):
    """
    result should run east to west, then north to south

    could be done more easily if a reliable NC file already exists:
            nc_dset = Dataset(nc_file, mode='r')
            lat = nc_dset.variables['grid_latitude'][indx_nrth,indx_east]
            lon = nc_dset.variables['grid_longitude'][indx_nrth,indx_east]
            nc_dset.close()
    """

    # might already exist
    # ===================
    if read_flag:
        if isfile(meteo_fname):
            meteo_df = read_csv(meteo_fname)
            print('\nRead {} lat/lons from {}'.format(len(meteo_df), meteo_fname))
            return meteo_df.copy()
        else:
            print('\nCould not find ' + meteo_fname + ' - will create')

    if not isfile(cnvrtd_coords_fname):
        print('File ' + cnvrtd_coords_fname + ' does not exist, cannot continue')
        return None

    cnvrt_hdrs = ['easting', 'northing', 'grid_ref', 'lat', 'lon']
    cnvrt_df = read_csv(cnvrtd_coords_fname, sep = ',', names = cnvrt_hdrs)
    num_total = len(cnvrt_df)

    hwsd_hdrs  = ['gran_lat', 'gran_lon', 'mu_global', 'lat', 'lon']
    hwsd_frame = read_csv(hwsd_fname, sep = ',', names = hwsd_hdrs)
    print('Read HWSD AOI file: ' + hwsd_fname + ' comprising {} grid cells'.format(len(hwsd_fname)))

    eastings_list  = [float(val) for val in easting_keys[1:]]
    northings_list = [float(val) for val in northing_keys]

    lats = []; lons = []; nrthngs = []; eastngs = []; cell_lons = []; cell_lats = []
    nwrite = 0
    num_nodata = 0
    last_time  = time()
    for indx_east, easting in enumerate(eastings_list):
        for indx_nrth, northing in enumerate(northings_list):

            locat = cnvrt_df.loc[(cnvrt_df['easting'] == easting) & (cnvrt_df['northing'] == northing)]
            if locat.empty:
                num_nodata += 1
            else:
                lat = locat['lat'].values[0]
                lon = locat['lon'].values[0]
                gran_lon = int(round((180.0 + lon)*NGRANULARITY))
                gran_lat = int(round((90.0  - lat)*NGRANULARITY))
                hwsd_locat = hwsd_frame.loc[(hwsd_frame['gran_lat'] == gran_lat) & (hwsd_frame['gran_lon'] == gran_lon)]
                if not hwsd_locat.empty:
                    lats.append(round(lat,4))
                    lons.append(round(lon,4))
                    eastngs.append(int(easting))
                    nrthngs.append(int(northing))
                    cell_lon = gran_lon/NGRANULARITY - 180.0
                    cell_lons.append(round(cell_lon,4))
                    cell_lat = 90.0 - gran_lat/NGRANULARITY
                    cell_lats.append(round(cell_lat,4))
                    nwrite += 1

            last_time = update_progress_post(last_time, num_nodata, nwrite, num_total,
                                                                                round(lat,3), round(lon,3))
    meteo_df = DataFrame()
    meteo_df['lat']  = lats
    meteo_df['lon']  = lons
    meteo_df['cell_lat'] = cell_lats
    meteo_df['cell_lon'] = cell_lons
    meteo_df['northing'] = nrthngs
    meteo_df['easting']  = eastngs
    meteo_df = meteo_df.sort_values(by=['cell_lat', 'cell_lon'], ascending = [False, True])

    try:
        meteo_df.to_csv (meteo_fname, index = False, header = True)
    except PermissionError as err:
        print('Could not create ' + meteo_fname + ' due to: ' + str(err) )

    print('\nWrote {} lat/lons to {}'.format(nwrite, meteo_fname))
    return meteo_df
