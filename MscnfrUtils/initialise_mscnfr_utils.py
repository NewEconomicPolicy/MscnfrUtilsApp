"""
#-------------------------------------------------------------------------------
# Name:        initialise_funcs.py
# Purpose:     script to read read and write the setup and configuration files
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
# setup file in CWD:            reform_hwsd_setup.txt
# config file in config_dir:    hwsd_soils_modis_config.txt
# MODIS land use file name:     <current working directory>\Docs\MODIS_land_uses.json
"""

__prog__ = 'initialise_mscnfr_utils.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import isdir, split, exists, join, isfile
from os import makedirs, mkdir, getcwd
from json import JSONDecodeError, load as json_load, dump as json_dump
from time import time, sleep
import sys
from set_up_logging import set_up_logging
from mscnfr_utils_fns import identify_nc_files, identify_ukcp18_dirs

PROGRAM_ID = 'mscnfr_utils'
ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '
sleepTime = 3

def initiation(form):
    '''
    this function is called to initiate the programme to process non-GUI settings.
    '''

    # retrieve settings
    # =================
    form.settings = _read_setup_file(PROGRAM_ID)
    set_up_logging(form, PROGRAM_ID)
    print('\nWill write logfile to: ' + form.lgr.handlers[0].baseFilename)

    return

def _read_setup_file(program_id):
    """
    read settings used for program from the setup file, if it exists,
    or create setup file using default values if file does not exist
    """
    func_name =  __prog__ +  ' _read_setup_file'

    fname_setup = program_id + '_setup.json'
    setup_file = join(getcwd(), fname_setup)

    if exists(setup_file):
        try:
            with open(setup_file, 'r') as fsetup:
                setup = json_load(fsetup)
        except (OSError, IOError) as err:
                print(err)
                sleep(sleepTime)
                sys.exit(0)
    else:
        setup = _write_default_setup_file(setup_file)

    settings = setup['setup']
    settings_list = ['chess_dir', 'config_dir', 'fname_png', 'log_dir', 'root_dir']
    for key in settings_list:
        if key not in settings:
            print(ERROR_STR + 'setting {} is required in setup file {} '.format(key, setup_file))
            sleep(sleepTime)
            exit(0)

    fname_png  = settings['fname_png']
    if not isfile(fname_png):
        print('Image file ' + fname_png + ' does not exist')

    settings['config_file'] = join(settings['config_dir'], program_id + '_config.json')
    settings['conversions'] = list(['CRU', 'Hadley', 'UK CP18'])
    settings['perturbs'] = None
    settings['hadley_year_sets'] = None
    settings['hadley_sets_mapped'] = None
    settings['cru_fnames'] = []

    return settings

def _write_default_setup_file(setup_file):
    """
    stanza if setup_file needs to be created
    """

    root_dir = 'E:\\Miscanfor'
    if not isdir(root_dir):
        root_dir = 'C:\\Miscanfor'

    log_dir    = root_dir + '\\logs'
    if not isdir(log_dir):
        makedirs(log_dir)

    nc_dir   = root_dir + '\\netCDF_4.01_met'
    if not isdir(nc_dir):
        makedirs(nc_dir)

    _default_setup = {
        'setup': {
            'root_dir'   : root_dir,
            'fname_png'  : join(root_dir + '\\Images', 'mxg_field.jpg'),
            'log_dir'    : log_dir,
            'config_dir' : root_dir + '\\config'
        }
    }
    # if setup file does not exist then create it...
    with open(setup_file, 'w') as fsetup:
        json_dump(_default_setup, fsetup, indent=2, sort_keys=True)
        return _default_setup

def read_config_file(form):
    """
    read widget settings used in the previous programme session from the config file, if it exists,
    or create config file using default settings if config file does not exist
    """
    func_name =  __prog__ +  ' read_config_file'

    config_file = form.settings['config_file']
    if exists(config_file):
        try:
            with open(config_file, 'r') as fconfig:
                config = json_load(fconfig)
        except (OSError, IOError, JSONDecodeError) as err:         # does not catch all errors
            print(str(err))
            return False
    else:
        config = _write_default_config_file(config_file)

    expectation = {'fnames':['csv_file', 'nc_dir', 'out_dir', 'ukcp18_dir'],
                   'run_settings':['min_lat', 'max_lat', 'max_cells', 'ignore_near', 'dummy_run', 'last_seq']}

    integrity_flag, mess = _check_two_level_json_file(config, expectation)
    if not integrity_flag:
            print(ERROR_STR + mess + ' in configuration file: ' + config_file)
            sleep(sleepTime)
            sys.exit(0)

    grp = 'fnames'
    csv_file = config[grp]['csv_file']
    out_dir = config[grp]['out_dir']
    nc_dir = config[grp]['nc_dir']
    ukcp18_dir = config[grp]['ukcp18_dir']

    grp = 'run_settings'
    min_lat     = config[grp]['min_lat']
    max_lat     = config[grp]['max_lat']
    max_cells   = config[grp]['max_cells']
    ignore_near = config[grp]['ignore_near']
    dummy_run   = config[grp]['dummy_run']
    last_seq    = config[grp]['last_seq']

    # set check boxes
    # ===============
    if ignore_near:
        form.w_ignore_near.setCheckState(2)
    else:
        form.w_ignore_near.setCheckState(0)

    if dummy_run:
        form.w_dummy.setCheckState(2)
    else:
        form.w_dummy.setCheckState(0)

    if last_seq:
        form.w_last_seq.setCheckState(2)
    else:
        form.w_last_seq.setCheckState(0)

    # retrieve the set of NC files
    # ============================
    form.w_lbl_ncdir.setText(nc_dir)
    identify_nc_files(form, nc_dir)

    form.lbl01.setText(csv_file)
    form.w_lbl_ukcp18.setText(ukcp18_dir)
    form.w_lbl_out_dir.setText(out_dir)
    form.w_min_lat.setText(min_lat)
    form.w_max_lat.setText(max_lat)
    form.w_max_cells.setText(max_cells)
    form.w_perturbs.setText(identify_ukcp18_dirs(form))

    _check_chess_wthr(form.settings)
    if form.settings['chess_hist'] is None:
        form.w_lookup.setEnabled(False)
        form.w_chess_csv.setEnabled(False)

    return True

def _check_chess_wthr(sttngs):
    """
    validate CHESS data installation
    """
    rcp_dir = join(sttngs['chess_dir'], 'CHESS_RCPs', 'Monthly')
    hist_dir = join(sttngs['chess_dir'], 'CHESS_historic', 'Monthly')
    lookup_dir = join(split(hist_dir)[0], 'lookup_table')
    # lat_lon_osgb_fn = join(lookup_dir, 'meteo_lat_lon_osgb.csv')
    lat_lon_osgb_fn = join(lookup_dir, 'CHESS_hwsd_lkup_tble.csv')
    brit_osgb_fn = join(lookup_dir, 'BritishGrid_HWSD.csv')

    if not isdir(hist_dir) or not isdir(rcp_dir) or not isdir(lookup_dir) or not isfile(brit_osgb_fn):
        mess = WARN_STR + 'Check CHESS installation:\n\thistory dir: ' + hist_dir + '\n\tRCP dir: ' + rcp_dir
        print(mess + '\n\tlookup dir: ' + lookup_dir + '\n\tLookup OSGB file: ' + brit_osgb_fn)
        hist_dir, rcp_dir, lookup_dir, lat_lon_osgb_fn, brit_osgb_fn = 5*[None]

    sttngs['chess_hist'] = hist_dir
    sttngs['chess_rcp'] = rcp_dir
    sttngs['chess_lkup'] = lookup_dir
    sttngs['lat_lon_osgb'] = lat_lon_osgb_fn
    sttngs['brit_osgb_fn'] = brit_osgb_fn

    return

def write_config_file(form):
    """

    """
    config_file = form.settings['config_file']
    config = {
        'fnames' : {
            'csv_file':form.lbl01.text(),
            'out_dir' : form.w_lbl_out_dir.text(),
            'nc_dir'  : form.w_lbl_ncdir.text(),
            'ukcp18_dir': form.w_lbl_ukcp18.text()
        },
        'run_settings': {
            'min_lat'   : form.w_min_lat.text(), # set limits to restrict the search when traversing the NetCDF files
            'max_lat'   : form.w_max_lat.text(),
            'max_cells' : form.w_max_cells.text(),
            'ignore_near' : form.w_ignore_near.isChecked(),
            'dummy_run'   : form.w_dummy.isChecked(),
            'last_seq'    : form.w_last_seq.isChecked()
        }
    }
    with open(config_file, 'w') as fconfig:
        json_dump(config, fconfig, indent=2, sort_keys=True)
        print('Updated ' + config_file)

    return True

def _check_two_level_json_file(config, expectation):
    """

    """
    integrity_flag = True
    mess = ''
    for grp in expectation:
        if grp in config:
            for attrib in expectation[grp]:
                if attrib not in config[grp]:
                    mess += 'attribute ' + attrib + ' not present'
                    integrity_flag = False
                    break
        else:
            mess += 'group ' + grp + ' not present'
            integrity_flag = False
            break

    return  integrity_flag, mess

def _write_default_config_file(config_file):
    """
    stub
    """
    config_dir, dummy = split(config_file)
    if not isdir(config_dir):
        mkdir(config_dir)

    _default_config = {
         'fnames' : {
            'csv_file' : '',
            'out_dir'  : '',
            'nc_dir'   : '',
            'ukcp18_dir': ''
        },
        'run_settings': {
            'min_lat'   : "-56.0",        # set limits to restrict the search when traversing the NetCDF files
            'max_lat'   : "84.0",
            'max_cells' : "999999999",
            'ignore_near' : True,
            'dummy_run'   : True,
            'last_seq'    : True
        }
    }
    # if config file does not exist then create it...
    with open(config_file, 'w') as fconfig:
        json_dump(_default_config, fconfig, indent=2, sort_keys=True)
        return _default_config

    return True
