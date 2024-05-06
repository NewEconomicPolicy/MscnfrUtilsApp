#-------------------------------------------------------------------------------
# Name:        MiscanforUtilsGUI.py
# Purpose:     GUI to identify HWSD cells for selected types and write a CSV of valid HWSD cells
# Author:      Mike Martin
# Created:     16/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#

__prog__ = 'MiscanforUtilsGUI.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from os.path import normpath
import sys

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QLabel, QWidget, QApplication, QHBoxLayout, QVBoxLayout, QGridLayout, QLineEdit,
                                                                     QComboBox, QPushButton, QCheckBox, QFileDialog)

from initialise_mscnfr_utils import read_config_file, initiation, write_config_file
from cru_NC_fns import generate_cru_csv_files, check_csv_file
from validate import check_soil_file
from ukcp18_fns import generate_ukcp18_csv_files
from hadley_NC_fns import generate_hadley_csv_files
from mscnfr_utils_fns import reformat_csv_files, identify_ukcp18_dirs, identify_nc_files
from hwsd_csv_to_osgb_lookup import write_lookup_from_hwsd_xlsx, write_osgb_meteogrid, remove_cells_from_hwsd

from merge_pt2 import write_merge_pt2
from chess_to_lookup_table_csv import write_lookup_csv, write_check_csv
from chess_to_mscnfr_txt_frmt import read_chess_nc_write_meteo_csv

STD_FLD_SIZE_60 = 60
STD_FLD_SIZE_80 = 80
STD_FLD_SIZE_120 = 120

class Form(QWidget):

    def __init__(self, parent=None):

        super(Form, self).__init__(parent)

        self.version = 'NetCDF'
        initiation(self)
        # define two vertical boxes, in LH vertical box put the painter and in RH put the grid
        # define horizon box to put LH and RH vertical boxes in
        hbox = QHBoxLayout()
        hbox.setSpacing(10)

        # left hand vertical box consists of png image
        # ============================================
        lh_vbox = QVBoxLayout()

        # LH vertical box contains image only
        lbl20 = QLabel()
        pixmap = QPixmap(self.settings['fname_png'])
        lbl20.setPixmap(pixmap)

        lh_vbox.addWidget(lbl20)

        # add LH vertical box to horizontal box
        hbox.addLayout(lh_vbox)

        # right hand box consists of combo boxes, labels and buttons
        # ==========================================================
        rh_vbox = QVBoxLayout()

        # The layout is done with the QGridLayout
        grid = QGridLayout()
        grid.setSpacing(10)	# set spacing between widgets

        rh_vbox.addLayout(grid)     # add grid to RH vertical box
        hbox.addLayout(rh_vbox)     # vertical box goes into horizontal box
        self.setLayout(hbox)        # the horizontal box fits inside the window

        self.setGeometry(300, 300, 300, 250)   # posx, posy, width, height
        self.setWindowTitle(self.version + ' - Generate CSV of HWSD mu_globals based on soil types and land use')

        # ======
        irow = 1
        w_csv_file = QPushButton("CSV file to check")
        helpText = 'Option to select CSV text file containing the data to be checked'
        w_csv_file.setEnabled(False)
        w_csv_file.setToolTip(helpText)
        grid.addWidget(w_csv_file, irow, 0)
        w_csv_file.clicked.connect(self.checkCsvFileClicked)

        lbl01 = QLabel()
        grid.addWidget(lbl01, irow, 1, 1, 5)        # row, column, rowSpan, columnSpan
        self.lbl01 = lbl01

        # =============
        irow += 2
        w_ukcp18_dir = QPushButton("UK CP18 location")
        helpText = 'Select directory for the UK CP18 project under which there are 5 sub-directories\n'
        helpText += 'for each of 6 metrics. Under each sub-directory there are sets of 12 CSV files,\n'
        helpText += 'each one describing a realisation or perturbation of the RCP 8.5 scenario e.g. p01113'
        w_ukcp18_dir.setToolTip(helpText)
        grid.addWidget(w_ukcp18_dir, irow, 0)  # Provinces hwsd directory
        w_ukcp18_dir.clicked.connect(self.fetchUkCp18Dir)

        w_lbl_ukcp18 = QLabel()
        grid.addWidget(w_lbl_ukcp18, irow, 1, 1, 5)
        self.w_lbl_ukcp18 = w_lbl_ukcp18

        irow += 1
        lbl04 = QLabel()
        grid.addWidget(lbl04, irow, 0, 1, 6)
        self.w_perturbs = lbl04

        irow += 1
        grid.addWidget(QLabel(), irow, 0)   # spacer

        # ======
        irow += 1
        w_nc_dir = QPushButton("Location of NC files")
        helpText = 'Select a directory where NetCDF files are located'
        w_nc_dir.setToolTip(helpText)
        grid.addWidget(w_nc_dir, irow, 0)  # Provinces hwsd directory
        w_nc_dir.clicked.connect(self.fetchNCsDir)

        w_lbl_ncdir = QLabel()
        grid.addWidget(w_lbl_ncdir, irow, 1, 1, 5)
        self.w_lbl_ncdir = w_lbl_ncdir

        # num NC files
        # ============
        irow += 1
        w_lbl_num_ncs = QLabel('')
        grid.addWidget(w_lbl_num_ncs, irow, 1, 1, 4)
        self.w_lbl_num_ncs = w_lbl_num_ncs
        
        # ==============================
        irow += 1
        lbl08a = QLabel('Min latitude:')
        lbl08a.setAlignment(Qt.AlignRight)
        helpText = 'any data less than this latitude will be ignored'
        lbl08a.setToolTip(helpText)
        grid.addWidget(lbl08a, irow, 0)

        w_min_lat = QLineEdit()
        w_min_lat.setFixedWidth(STD_FLD_SIZE_60)
        grid.addWidget(w_min_lat, irow, 1)
        self.w_min_lat = w_min_lat

        lbl08b = QLabel('Max latitude:')
        lbl08b.setAlignment(Qt.AlignRight)
        helpText = 'any data greater than this latitude will be ignored'
        lbl08b.setToolTip(helpText)
        grid.addWidget(lbl08b, irow, 2)

        w_max_lat = QLineEdit()
        w_max_lat.setFixedWidth(STD_FLD_SIZE_60)
        grid.addWidget(w_max_lat, irow, 3)
        self.w_max_lat = w_max_lat

        irow += 1
        grid.addWidget(QLabel(), irow, 0)   # spacer

        # =============
        irow += 1
        w_out_dir = QPushButton('Output Directory')
        helpText = 'Directory into which CSV files will be written'
        w_out_dir.setToolTip(helpText)
        grid.addWidget(w_out_dir, irow, 0)
        w_out_dir.clicked.connect(self.fetchOutDir)
        self.w_out_dir = w_out_dir

        w_lbl_out_dir = QLabel('')     # Output Directory
        grid.addWidget(w_lbl_out_dir, irow, 1, 1, 5)
        self.w_lbl_out_dir = w_lbl_out_dir

        irow += 1
        grid.addWidget(QLabel(), irow, 0)   # spacer

        # user can restrict length of run
        # ===============================
        irow += 1
        lbl08c = QLabel('Max grid cells:')
        lbl08c.setAlignment(Qt.AlignRight)
        helpText = 'Max data blocks or grid cells with data'
        lbl08c.setToolTip(helpText)
        grid.addWidget(lbl08c, irow, 0)

        w_max_cells = QLineEdit()
        w_max_cells.setFixedWidth(STD_FLD_SIZE_60)
        helpText = 'No output will be written'
        w_max_cells.setToolTip(helpText)
        grid.addWidget(w_max_cells, irow, 1)
        self.w_max_cells = w_max_cells

        # =======
        w_dummy = QCheckBox('Dummy run')
        w_dummy.setChecked(True)
        helpText = 'Applies to writeMergePt2: No output will be written\n\t'
        helpText += 'and CSVs from CHESS: write meteogrid and one metric file only'
        w_dummy.setToolTip(helpText)
        grid.addWidget(w_dummy, irow, 2)
        self.w_dummy = w_dummy

        w_last_seq = QCheckBox('Last sequence')
        w_last_seq.setChecked(True)
        helpText = 'Process last sequence from the set of HWSD sequeces only'
        w_last_seq.setToolTip(helpText)
        w_last_seq.setFixedWidth(STD_FLD_SIZE_120)
        grid.addWidget(w_last_seq, irow, 3)
        self.w_last_seq = w_last_seq

        w_ignore_near =  QCheckBox('Ignore near')
        w_ignore_near.setChecked(True)
        helpText = 'Skip cells where is not an exact match for a cell coordinate from sequence file\n' \
                                                        + ' and cell coordinates from soil data files'
        w_ignore_near.setToolTip(helpText)
        w_ignore_near.setFixedWidth(STD_FLD_SIZE_120)
        grid.addWidget(w_ignore_near, irow, 4)
        self.w_ignore_near = w_ignore_near

        # =======
        irow += 1
        helpText = 'Set type of conversion:\n\tCRU: \n\tHadley: \n\tUK CP18: '
        lbl18 = QLabel('CSV conversion:')
        lbl18.setAlignment(Qt.AlignRight)
        lbl18.setToolTip(helpText)
        grid.addWidget(lbl18, irow, 0)

        combo18 = QComboBox()
        for conversion in self.settings['conversions']:
            combo18.addItem(conversion)
        combo18.setToolTip(helpText)
        combo18.setCurrentIndex(2)
        combo18.setFixedWidth(STD_FLD_SIZE_80)
        combo18.setToolTip(helpText)
        grid.addWidget(combo18, irow, 1)
        self.w_combo18 = combo18

        irow += 1
        grid.addWidget(QLabel(), irow, 0)   # spacer

        # first command line
        # ==================
        irow += 1
        w_create_file = QPushButton("Generate CSV files")
        helpText = 'Generate a new Miscanfor CSV formatted metric and meteogrid files from NetCDF or CSV files'
        helpText += '\n\tuses selection from CSV conversion drop down'
        w_create_file.setToolTip(helpText)
        w_create_file.setFixedWidth(STD_FLD_SIZE_120)
        w_create_file.clicked.connect(self.genCsvFilesClicked)
        grid.addWidget(w_create_file, irow, 1)

        w_write_mrge = QPushButton("writeMergePt2")
        helpText = 'Run writeMergePt2 soil data preparation program which downscales and merges the\n' + \
                            'landuse from IGBP & HWSD soil data with the lat, long and soil data parameters\n' + \
                            ' output from HWSD_extract13.f90'
        w_write_mrge.setToolTip(helpText)
        w_write_mrge.setFixedWidth(STD_FLD_SIZE_120)
        w_write_mrge.clicked.connect(self.writeMergePt2)
        grid.addWidget(w_write_mrge, irow, 2)

        w_check_csv = QPushButton("Check CSV file")
        helpText = 'Check extents of very large CSV file'
        w_check_csv.setToolTip(helpText)
        w_check_csv.clicked.connect(self.checkCsvFileClicked)
        w_check_csv.setEnabled(False)
        grid.addWidget(w_check_csv, irow, 3)

        w_reform_csv = QPushButton("Reformat CSVs")
        helpText = 'Read list of HWSD files and write lat, long, mu_global files for Miscanfor'
        helpText += '\n\tuses module mscnfr_utils_fns'
        w_reform_csv.setToolTip(helpText)
        w_reform_csv.clicked.connect(self.reformatCsvFilesClicked)
        grid.addWidget(w_reform_csv, irow, 4)

        w_save = QPushButton("Save", self)
        w_save.clicked.connect(self.saveClicked)
        grid.addWidget(w_save, irow, 5)

        # second command line
        # ===================
        irow += 1

        w_check_soil = QPushButton("Check soil file")
        helpText = 'Check CSV soil file integrity'
        w_check_soil.setToolTip(helpText)
        w_check_soil.clicked.connect(self.checkSoilFileClicked)
        # w_check_soil.setEnabled(False)
        grid.addWidget(w_check_soil, irow, 3)

        w_chck_ukcp18 = QPushButton("Check UKCP18 outputs")
        helpText = 'Check UKCP18 outputs'
        w_chck_ukcp18.setToolTip(helpText)
        w_chck_ukcp18.clicked.connect(self.checkUkCp18Clicked)
        w_chck_ukcp18.setEnabled(False)
        grid.addWidget(w_chck_ukcp18, irow, 4)

        # third command line for CHESS
        # ============================
        irow += 1
        w_cvrtcoord = QPushButton("Check OSGB to WGS84")
        helpText = 'Test Hannah Fry Python module to convert OSGB eastings and northings '
        helpText += '\nto the World Geodetic System (WGS) lon/lat coordinate system and vice-versa.'
        helpText += '\nCreates an output file: cvrtcoord.csv'
        w_cvrtcoord.setToolTip(helpText)
        w_cvrtcoord.setFixedWidth(STD_FLD_SIZE_120)
        w_cvrtcoord.clicked.connect(self.checkCvrtCoord)
        grid.addWidget(w_cvrtcoord, irow, 0)

        w_lookup = QPushButton("Lookup from CHESS")
        helpText = 'Write meteo_lat_lon_osgb.csv lookup and meteogrid files from CHESS NetCDF'
        w_lookup.setToolTip(helpText)
        w_lookup.setFixedWidth(STD_FLD_SIZE_120)
        w_lookup.clicked.connect(self.genLookupFromChess)
        # w_lookup.setEnabled(False)
        grid.addWidget(w_lookup, irow, 1)

        w_lkup_hwsd = QPushButton("Lookup from HWSD")
        helpText = 'Write CHESS_hwsd_lkup_tble.csv based HWSD csv file and compress so there are no duplicated OSGB indices'
        w_lkup_hwsd.setToolTip(helpText)
        w_lkup_hwsd.setFixedWidth(STD_FLD_SIZE_120)
        w_lkup_hwsd.clicked.connect(self.genLookupFromHwsd)
        grid.addWidget(w_lkup_hwsd, irow, 2)

        w_osgb_csv = QPushButton("Test meteogrid")
        helpText = 'Compress CHESS_hwsd_lkup_tble.csv'
        helpText = 'Write test meteogrid'
        # w_osgb_csv.setEnabled(False)
        w_osgb_csv.setToolTip(helpText)
        w_osgb_csv.setFixedWidth(STD_FLD_SIZE_120)
        w_osgb_csv.clicked.connect(self.genOsgbMeteogrid)
        grid.addWidget(w_osgb_csv, irow, 3)

        w_chess_csv = QPushButton("CSVs from CHESS")
        helpText = 'Write Miscanfor CSV formatted metric files from CHESS NetCDF using meteo_lat_lon_osgb.csv lookup file'
        w_chess_csv.setToolTip(helpText)
        w_chess_csv.setFixedWidth(STD_FLD_SIZE_120)
        w_chess_csv.clicked.connect(self.genCsvFromChess)
        grid.addWidget(w_chess_csv, irow, 4)

        w_exit = QPushButton("Exit", self)
        w_exit.clicked.connect(self.exitClicked)
        grid.addWidget(w_exit, irow, 6)

        # fourth command line for CHESS
        # =============================
        irow += 1
        w_filter_hwsd = QPushButton("Filter HWSD")
        helpText = 'Read BritishGrid_HWSD.csv and create copy with identified cells rmoved'
        w_filter_hwsd.setToolTip(helpText)
        w_filter_hwsd.setFixedWidth(STD_FLD_SIZE_120)
        w_filter_hwsd.clicked.connect(self.removeCellsFromHwsd)
        grid.addWidget(w_filter_hwsd, irow, 2)

        # reads and set values from last run
        # ==================================
        read_config_file(self)

    # ===================================================
    def removeCellsFromHwsd(self):
        """

        """
        remove_cells_from_hwsd(self)

        return

    def genOsgbMeteogrid(self):
        """
        test coordinates only
        """
        write_osgb_meteogrid(self)

        return

    def genLookupFromHwsd(self):
        """

        """
        write_lookup_from_hwsd_xlsx(self)

        return

    # ===================================================
    def genCsvFromChess(self):
        """

        """
        read_chess_nc_write_meteo_csv(self)

        return

    # ===================================================

    def checkCvrtCoord(self):
        """

        """
        write_check_csv(self)

        return

    # ===================================================
    def genLookupFromChess(self):
        """

        """
        write_lookup_csv(self)

        return

    def checkUkCp18Clicked(self):

        pass

    def checkSoilFileClicked(self):

        check_soil_file(self)

    def writeMergePt2(self):

        write_merge_pt2(self)

    def genCsvFilesClicked(self):
        '''

        '''
        conversion_type = self.w_combo18.currentText()
        if conversion_type == 'CRU':
            generate_cru_csv_files(self)
        elif conversion_type == 'Hadley':
            generate_hadley_csv_files(self)
        else:
            generate_ukcp18_csv_files(self)

    def checkCsvFileClicked(self):

        check_csv_file(self)

    def reformatCsvFilesClicked(self):

        num_created = reformat_csv_files(self)
        print('Number of files created: {}'.format(num_created))

    def fetchUkCp18Dir(self):
        """
        """
        fname = self.w_lbl_ukcp18.text()
        fname = QFileDialog.getExistingDirectory(self, 'Select directory', fname)
        if fname != '':
            fname = normpath(fname)
            self.w_lbl_ukcp18.setText(fname)
            self.w_perturbs.setText(identify_ukcp18_dirs(self))

    def fetchNCsDir(self):
        """
        """
        fname = self.w_lbl_ncdir.text()
        fname = QFileDialog.getExistingDirectory(self, 'Select directory', fname)
        if fname != '':
            fname = normpath(fname)
            self.w_lbl_ncdir.setText(fname)
            identify_nc_files(self, fname)

    def fetchOutDir(self):

        fname = self.w_lbl_out_dir.text()
        fname = QFileDialog.getExistingDirectory(self, 'Select directory', fname)
        if fname != '':
            # TODO: need to check this is a directory with write permissions
            fname = normpath(fname)
            self.w_lbl_out_dir.setText(fname)

    def saveClicked(self):
        '''
        write last GUI selections
        '''
        write_config_file(self)

    def exitClicked(self):
        '''
        '''
        # write last GUI selections
        write_config_file(self)
        self.lgr.handlers[0].close()
        self.close()

def main():
    """
    """
    app = QApplication(sys.argv)  # create QApplication object
    form = Form() # instantiate form
    # display the GUI and start the event loop if we're not running batch mode
    form.show()             # paint form
    sys.exit(app.exec_())   # start event loop

if __name__ == '__main__':
    main()
