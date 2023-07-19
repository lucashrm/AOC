#!/usr/bin/env python3

''' --------------------------------------------------------------------------

Master GUI for the control of CIAO, the Calern Imaging Adaptive Observatory.

Author: Frantz Martinache (frantz.martinache@oca.eu)
-------------------------------------------------------------------------- '''

from pkgutil import get_data
from xaosim.QtMain import QtMain, QApplication
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import QThread, Qt, QDir
from PyQt5.QtWidgets import QLabel, QFileDialog, QWidget, QVBoxLayout, QInputDialog, QLineEdit, QGraphicsWidget
from PyQt5.QtGui import QImage, QPainter, QPen, QBrush, QGuiApplication
from PyQt5.QtChart import QChart, QChartView, QBarSet, QBarCategoryAxis, QBarSeries, QValueAxis

import ctypes
import _ctypes


import threading

import pyqtgraph as pg
# from pyqtgraph.GraphicsScene import MouseClickEvent
from pyqtgraph.GraphicsScene.mouseEvents import HoverEvent, MouseClickEvent, MouseDragEvent

import sys
import numpy as np
from numpy.linalg import solve
import matplotlib.cm as cm

import os
import shutil

import pdb
import time

import astropy.io.fits as pf
import argparse

import configparser
import io

import paho.mqtt.client as mqtt
import json

# =====================================================================
home = os.getenv('HOME')
conf_dir = home+'\\.config\\ciao\\'
if not os.path.exists(conf_dir):
    os.makedirs(conf_dir)

ciao_home = os.getenv('CIAO_HOME')


guidef = "AOC_gui.ui"
try:
    shutil.copy(guidef, conf_dir+guidef)
except:
    print("Problem with config directory?")

zer_gui = "zernike_gui.ui"
try:
    shutil.copy(zer_gui, conf_dir+zer_gui)
except:
    print("Problem with config directory?")

# pup_gui = "AOC_pupil.ui"
pup_gui = "AOC_pupil_v02.ui"  ## WITH NUMBERING TRANSPOSED
try:
    shutil.copy(pup_gui, conf_dir+pup_gui)
except:
    print("Problem with config directory?")

# from zoffsets_window import *

zoffsets_gui = "AOC_ZOffsets.ui"  ## WITH NUMBERING TRANSPOSED
try:
    shutil.copy(zoffsets_gui, conf_dir+zoffsets_gui)
except:
    print("Problem with config directory?")

from AOC_wfs import WFS
from AOC_wfc import *
from AOC_cam import Cam
from AOC_tel import Tel

# =====================================================================
# =====================================================================
myqt = 0 # to have myqt as a global variable
logfile = ".\\ciao_log.html"

def main(argv):
    global myqt

    parser = argparse.ArgumentParser(description="C2PU AO control GUI")
    parser.add_argument("--simu", dest="simu", action="store_true",
                        help="runs in simulated shared memory environment")
    parser.add_argument("--ncpu", dest="ncpu", action="store", type=int, default=1,
                        help="runs using the specified number of CPUs")

    args = parser.parse_args()
    parser.set_defaults(simu=False)
    parser.print_help()
    print("\n")

    ncpu = False
    if args.ncpu > 1:
        ncpu = args.ncpu
        print("\nRUNNING THE WFS USING %d CPUs!\n" % (ncpu,))

    QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_EnableHighDpiScaling, False)  # enable highdpi scaling
    QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, False)  # use highdpi icons
    myqt = QtMain()
    gui = MyWindow(simu=args.simu, ncpu=ncpu)
    myqt.mainloop()
    myqt.gui_quit()
    sys.exit()

# =====================================================================
#                               Tools
# =====================================================================
def arr2im(arr, vmin=False, vmax=False, pwr=1.0, cmap=None, gamma=1.0):
    ''' --------------------------------------------------------------
    Convert a numpy array into image for display:

    limits dynamic range, power coefficient, applies colormap and gamma
    --------------------------------------------------------------  '''
    arr2 = arr.astype('float')
    if vmin is False:
        mmin = arr2.min()
    else:
        mmin = vmin

    if vmax is False:
        mmax = arr2.max()
    else:
        mmax = vmax

    arr2 -= mmin
    if mmax != mmin:
        arr2 /= (mmax-mmin)

    arr2 = arr2**pwr

    #cm.jet.set_gamma(0.1)
    if cmap == None:
        mycmap = cm.jet
    else:
        mycmap = cmap

    res = mycmap(arr2)
    res[:,:,3] = gamma
    return(res)

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

# =====================================================================
#                         Thread specifics
# =====================================================================
class GenericThread(QtCore.QThread):
    ''' ---------------------------------------------------
    generic thread class used to externalize the execution
    of a function (calibration, closed-loop) to a separate
    thread.
    --------------------------------------------------- '''
    def __init__(self, function, *args, **kwargs):
        QtCore.QThread.__init__(self)
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def __del__(self):
        self.wait()

    def run(self):
        self.function(*self.args,**self.kwargs)
        return

# =====================================================================
#                     Chart Window
# =====================================================================
class ChartWindow(QWidget):
    def __init__(self):
        super(ChartWindow, self).__init__()
        uic.loadUi("chart_ui2.ui", self)

        self.parentWindow = None




# =====================================================================
#                     DM Pupil definition GUI
# =====================================================================
class PupilWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """

    class SubAperture:
        def __init__(self, _IDnumber = 0, _bEnabled = False, _threshold = 0.):
            self.IDnumber = _IDnumber
            self.bEnabled = _bEnabled           # bEnabled SHOULD BE AN INT32 VALUE
            self.threshold = _threshold

    def __init__(self, pView):
        # super().__init__()
        # layout = QVBoxLayout()
        # self.label = QLabel("DM Pupil Definition Tool")
        # layout.addWidget(self.label)
        # self.setLayout(layout)

        # PupilWindow(self,pupmask=self.pupil_mask.flatten(),vCurCellThreshold=self.vCurCellThreshold,vCurEnabled=self.vCurEnabled)
        super(PupilWindow, self).__init__()
        if not os.path.exists(conf_dir + pup_gui):
            uic.loadUi(pup_gui, self)
        else:
            uic.loadUi(conf_dir + pup_gui, self)

        self.pView = pView
        print('self.pView.pupil_mask >> ', self.pView.pupil_mask)

        self.pView.bShowCellSelectionFeedback = True

        self.sel_type = 0

        # self.full_pupil = pupmask
        # self.pupmask = pupmask
        self.vMapToSquare = np.abs(pView.pupil_mask > 0).nonzero()#.astype(np.int32)

        print('self.vMapToSquare >> ', self.vMapToSquare)
        
        self.nSubAps = len(self.pView.pupil_mask.flatten())
        print('self.nSubAps >> ', self.nSubAps)

        self.pB_Select.clicked.connect(self.select)
        self.pB_SelectEnabled.clicked.connect(self.select_enabled)
        self.pB_DeselectAll.clicked.connect(self.deselect_all)
        self.pB_Invert.clicked.connect(self.invert)
        self.pB_Enable.clicked.connect(self.enable)
        self.pB_Disable.clicked.connect(self.disable)
        self.pB_SetThreshold.clicked.connect(self.set_threshold)
        self.pB_Validate.clicked.connect(self.apply_threshold_map)
        self.pB_SetAsDefault.clicked.connect(self.save_as_default)

        # CONNECT ALL SUBAPERTURES CHECKS TO THE SAME FUNCTION
        for k in range(self.nSubAps):
            exec("self.check_cell_%d.stateChanged[int].connect(self.update_selected_subaps)" % (k+1,))
        
        self.vSelectedSupAps = []

        print('self.pView.CurSubApSel >> ',self.pView.vCurSubApSel.size)

        # CREATE A LIST WITH ALL SUBAPERTURE CHECKBOXES
        
        self.vSubApsChecks = [self.check_cell_1]
        for k in range(2, self.nSubAps+1):
            exec("self.vSubApsChecks.append(self.check_cell_%d)" % (k,))

        self.vSubAp = []

        if (self.pView.bSubApsThresholdsFileExits == False) and (self.pView.bSubApsEnabledFileExits == False):
            tmpArr = self.pView.pupil_mask.flatten()
            for k in range(self.nSubAps):
                mySA = self.SubAperture(_IDnumber = k, _threshold = self.pView.wfs.threshold)
                item = self.vSubApsChecks[k]
                if tmpArr[k] == 1:
                    self.pView.vCurSubApSel[k] = True
                    self.pView.vSubApEnable[k] = 1
                    self.pView.vSubApThresh[k] = mySA.threshold
                    mySA.bEnabled = True
                    item.setChecked(True)
                    item.setStyleSheet('color: green')
                    # item.setStyleSheet('color: white; background-color: green; selection-color: black; selection-background-color: red;')
                    # item.setToolTip(mySA.threshold)
                else:
                    self.pView.vCurSubApSel[k] = False
                    self.pView.vSubApEnable[k] = 0
                    mySA.bEnabled = False
                    item.setChecked(False)
                    item.setStyleSheet('color: red')
                    # item.setStyleSheet('QCheckBox::indicator {color: black; background-color:red; selection-color: yellow; selection-background-color: blue;}')
                    item.setStyleSheet('color: white; background-color: red; selection-color: black; selection-background-color: red;')

                self.vSubAp.append( mySA )

        elif (self.pView.bSubApsThresholdsFileExits == True) and (self.pView.bSubApsEnabledFileExits == True):
            for k in range(self.nSubAps):
                mySA = self.SubAperture(_IDnumber = k, _bEnabled = self.pView.vSubApEnable[k], _threshold = self.pView.vSubApThresh[k])
                self.vSubAp.append( mySA )

                item = self.vSubApsChecks[k]

                if mySA.bEnabled == True:
                    item.setChecked(True)
                    item.setStyleSheet('color: green')
                else:
                    item.setChecked(False)
                    item.setStyleSheet('color: red')
                    # item.setStyleSheet('QCheckBox::indicator {color: black; background-color:red; selection-color: yellow; selection-background-color: blue;}')
                    item.setStyleSheet('color: white; background-color: red; selection-color: black; selection-background-color: red;')

        self.combo_CellGroup.currentIndexChanged.connect(self.comboChange)
        self.combo_CellGroup.addItem("All")
        self.combo_CellGroup.addItem("None")
        self.combo_CellGroup.addItem("Full pupil")
        self.combo_CellGroup.addItem("Current pupil")
        self.combo_CellGroup.addItem("All enabled")
        self.combo_CellGroup.addItem("All disabled")

    def closeEvent(self, *args, **kwargs):
        # super(QWidget, self).closeEvent(*args, **kwargs)
        # self.pView.bShowCellSelectionFeedback = False
        self.pView.check_ShowPupil.setChecked(False)
        print("you just closed the pyqt window!!! you are awesome!!!")    # def __init__(self, parentWindow, pupmask, vCurCellThreshold, vCurEnabled):

    def save_as_default(self,):
        self.apply_threshold_map()
        fits.writeto(self.pView.subaps_threshold_file,self.pView.vSubApThresh,overwrite=True)
        fits.writeto(self.pView.subaps_enabled_file,self.pView.vSubApEnable,overwrite=True)


    # SETS THE SELECTION MODE
    def comboChange(self,):
        self.sel_type = self.combo_CellGroup.currentIndex()
        print('sel >> ', self.sel_type)
        # if sel == 0:
        #     self.vSelectedSupAps = self.vSubApsChecks
        # elif sel == 1:
        #     self.vSelectedSupAps = []
        # elif sel == 2:
        #     self.vSelectedSupAps = []
        #     for k in self.vMapToSquare:
        #         exec("self.vSelectedSupAps.append(self.check_cell_%d)" % (k+1,))
        # elif sel == 3:
        #     print('Not implemented yet!')

    def update_selected_subaps(self,):
        count = 1
        self.vSelectedSupAps = []
        self.pView.vCurSubApSel[:] = False
        # for item in self.vSubApsChecks:
        #     if item.isChecked() == True:
        #         self.vSelectedSupAps.append(item)
        #         print('check_cell_%d is checked!' % (count))
        #     count = count + 1
        for k in range(self.nSubAps):
            item = self.vSubApsChecks[k]
            if item.isChecked() == True:
                self.vSelectedSupAps.append(item)
                self.pView.vCurSubApSel[k] = True

    def showSelectionFeedback(self):
        self.pView.bShowCellSelectionFeedback = True
        # if self.bShowSelectionFeedback == True:
            
            # for k in range(self.nSubAps):
            #     if self.pView.vCurSubApSel[k] == True:
            #         self.pView.xx
                
    def select(self,):
        if self.sel_type == 0:
            # All
            self.pView.vCurSubApSel[:] = True
            self.vSelectedSupAps[:] = self.vSubApsChecks[:]
            for item in self.vSubApsChecks:
                item.setChecked(True)
        elif self.sel_type == 1:
            # None
            self.pView.vCurSubApSel[:] = False
            self.vSelectedSupAps = []
            for item in self.vSubApsChecks:
                item.setChecked(False)
        elif self.sel_type == 2:
            # Full pupil
            for k in range(self.nSubAps):
                item = self.vSubApsChecks[k]
                tmpArr = self.pView.pupil_mask.flatten()
                if tmpArr[k] == 1:
                    self.pView.vCurSubApSel[k] = True
                    item.setChecked(True)
                else:
                    self.pView.vCurSubApSel[k] = False
                    item.setChecked(False)
        elif self.sel_type == 3:
            # Current pupil
            print('(Not implemented yet.)')
        elif self.sel_type == 4:
            # All enabled
            for k in range(self.nSubAps):
                item = self.vSubApsChecks[k]
                if self.vSubAp[k].bEnabled == True:
                    self.pView.vCurSubApSel[k] = True
                    item.setChecked(True)
                else:
                    self.pView.vCurSubApSel[k] = False
                    item.setChecked(False)
        elif self.sel_type == 5:
            # All disabled
            for k in range(self.nSubAps):
                item = self.vSubApsChecks[k]
                if self.vSubAp[k].bEnabled == False:
                    self.pView.vCurSubApSel[k] = True
                    item.setChecked(True)
                else:
                    self.pView.vCurSubApSel[k] = False
                    item.setChecked(False)
                
        self.update_selected_subaps()

        print('Number of selected sub-apertures >> ',len(self.vSelectedSupAps))

    def select_enabled(self,):
        self.sel_type = 4
        self.select()
        # for item in self.vSubApsChecks:
        #     item.setChecked(False)
        # self.update_selected_subaps()

    def deselect_all(self,):
        self.sel_type = 1
        self.select()
        # for item in self.vSubApsChecks:
        #     item.setChecked(False)
        # self.update_selected_subaps()

    def invert(self,):
        self.vSelectedSupAps = []
        for item in self.vSubApsChecks:
            if item.isChecked() == True:
                item.setChecked(False)
            else:
                item.setChecked(True)
        self.update_selected_subaps()


    def enable(self,):
        for k in range(self.nSubAps):
            item = self.vSubApsChecks[k]
            self.pView.vSubApEnable[:] = 0
            if item.isChecked() == True:
                self.vSubAp[k].bEnabled = True
                item.setStyleSheet( 'color: green' )
                item.setChecked(False)
                self.pView.vSubApEnable[k] = 1

    def disable(self,):
        for k in range(self.nSubAps):
            item = self.vSubApsChecks[k]
            self.pView.vSubApEnable[:] = 1
            if item.isChecked() == True:
                self.vSubAp[k].bEnabled = False
                # item.setStyleSheet( 'color: red' )
                item.setStyleSheet('color: white; background-color: red; selection-color: black; selection-background-color: red;')
                item.setChecked(False)
                self.pView.vSubApEnable[k] = 0
        

    def set_threshold(self,):
        thresh_val = self.spinB_threshold.value()
        for k in range(self.nSubAps):
            # exec('if self.vSelectedSupAps[k].isChecked() == True:')
            if self.pView.vCurSubApSel[k] == True:
                self.vSubAp[k].threshold = thresh_val
                self.pView.vSubApThresh[k] = thresh_val
        print('self.pView.vSubApThresh >> ', self.pView.vSubApThresh)
        
        # print('Not implemented yet!')

    def apply_threshold_map(self,):
        
        for k in range(self.nSubAps):
            if self.vSubAp[k].bEnabled == True:
                self.pView.vSubApEnable[k] = 1
            else:
                self.pView.vSubApEnable[k] = 0

            self.pView.vSubApThresh[k] = self.vSubAp[k].threshold

        self.pView.update_pupil_thresholds()

# =====================================================================
#                     Zernike GUI object
# =====================================================================

class ZernikeWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self):
        super(ZernikeWindow, self).__init__()
        if not os.path.exists(conf_dir + zer_gui):
            uic.loadUi(zer_gui, self)
        else:
            uic.loadUi(conf_dir + zer_gui, self)

        # layout = QVBoxLayout()
        # self.label = QLabel("Zernike Window")
        # layout.addWidget(self.label)
        # self.setLayout(layout)

        self.parentWindow = None

        self.max_modes = 100
        self.vZerSelection = np.zeros(self.max_modes,dtype='int16')

        # self.chB_zer_ord1_all.setChecked(Qt.Unchecked)

        self.chB_zer_ord1_all.stateChanged[int].connect(self.update_ZerSelection)
        self.chB_zer_ord2_all.stateChanged[int].connect(self.update_ZerSelection)
        self.chB_zer_ord3_all.stateChanged[int].connect(self.update_ZerSelection)
        self.chB_zer_ord4_all.stateChanged[int].connect(self.update_ZerSelection)
        self.chB_zer_ord5_all.stateChanged[int].connect(self.update_ZerSelection)
        self.chB_zer_ord6_all.stateChanged[int].connect(self.update_ZerSelection)

    def update_ZerSelection(self,):
        if self.chB_zer_ord1_all.isChecked() == True:
            for k in range(1,3):
                exec('self.chB_zer_ord1_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,3):
                exec('self.chB_zer_ord1_%d.setChecked(Qt.Unchecked)' % k)

        if self.chB_zer_ord2_all.isChecked() == True:
            for k in range(1,4):
                exec('self.chB_zer_ord2_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,4):
                exec('self.chB_zer_ord2_%d.setChecked(Qt.Unchecked)' % k)

        if self.chB_zer_ord3_all.isChecked() == True:
            for k in range(1,5):
                exec('self.chB_zer_ord3_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,5):
                exec('self.chB_zer_ord3_%d.setChecked(Qt.Unchecked)' % k)

        if self.chB_zer_ord4_all.isChecked() == True:
            for k in range(1,6):
                exec('self.chB_zer_ord4_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,6):
                exec('self.chB_zer_ord4_%d.setChecked(Qt.Unchecked)' % k)

        if self.chB_zer_ord5_all.isChecked() == True:
            for k in range(1,7):
                exec('self.chB_zer_ord5_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,7):
                exec('self.chB_zer_ord5_%d.setChecked(Qt.Unchecked)' % k)

        if self.chB_zer_ord6_all.isChecked() == True:
            for k in range(1,8):
                exec('self.chB_zer_ord6_%d.setChecked(Qt.Checked)' % k)
        else:
            for k in range(1,8):
                exec('self.chB_zer_ord6_%d.setChecked(Qt.Unchecked)' % k)

        count = 0
        for k in range(1,7):
            for i in range(1,k+2):
                exec('self.vZerSelection[count] = 1 if self.chB_zer_ord%d_%d.isChecked() == True else 0' % (k,i))
                count += 1

        self.parentWindow.ZER_WFC.mode_select = self.vZerSelection[0:self.parentWindow.ZER_WFC.nmodes]

        print('self.vZerSelection: ',self.vZerSelection)

# =====================================================================
#                          Main GUI object
# =====================================================================

class MyWindow(QtWidgets.QMainWindow):
    ''' ------------------------------------------------------
    This is the meat of the program: the class that drives
    the GUI.

    Parameters:
    ----------
    - simu: boolean (set to True to run in simulated env)
    - ncpu: integer (set value to use the required number of CPUs)
    ------------------------------------------------------ '''
    def __init__(self, simu=True, ncpu=False):
        global index
        self.bInitDone = False

        self.simu    = simu # simulation or real thing?
        self.shm_cam = None
        self.shm_wfs = None
        self.vmin    = False
        self.vmax    = False
        self.mycmap  = cm.magma
        self.pwr     = 1.0
        self.thresh  = 0
        self.cam_counter = -1
        self.tel_ttx_cor = 0.
        self.tel_tty_cor = 0.
        self.tel_ttx_t = 0.
        self.tel_tty_cor = 0.
        self.tel_tty_cor = 0.
        self.tel_loop_keepgoing = False
        self.tel_mon_period_sec = 10.
        self.OffloadPixScale = 0.5

        self.DMstaticUpdate_period_sec = 0.5
        self.bShowCellSelectionFeedback = False

        # =========================================
        #               offsets
        # =========================================
        self.ostep = 0.5 # size of the offset step
        self.xoff  = 0.0 # keep track of total offset commands
        self.yoff  = 0.0 # sent to the WFS

        self.wfc_enabled        = False
        self.wfc_status_updated = True  # control GUI updates
        self.wfc_calibrating    = ""    # control variable

        self.lib_EventManager = ctypes.CDLL(r'C:\Users\lucas\Documents\STAGE\Sources\x64\Release\EventManager.dll')

        self.lib_EventManager.register_named_event.restype = ctypes.c_int
        self.lib_EventManager.register_named_event.argtypes = [ctypes.c_char_p]

        self.lib_EventManager.unregister_named_event.restype = ctypes.c_int
        self.lib_EventManager.unregister_named_event.argtypes = [ctypes.c_char_p]

        self.lib_EventManager.unregister_all.restype = ctypes.c_int
        self.lib_EventManager.unregister_all.argtypes = []

        self.lib_EventManager.trigger_named_event.restype = ctypes.c_int
        self.lib_EventManager.trigger_named_event.argtypes = [ctypes.c_char_p]

        self.lib_EventManager.chk_wait_status_named_event.restype = ctypes.c_int
        self.lib_EventManager.chk_wait_status_named_event.argtypes = [ctypes.c_char_p]

        self.lib_EventManager.wait_for_single_event.restype = ctypes.c_int
        self.lib_EventManager.wait_for_single_event.argtypes = [ctypes.c_char_p, ctypes.c_ulong]

        self.lib_EventManager.wait_for_multiple_events.restype = ctypes.c_int
        self.lib_EventManager.wait_for_multiple_events.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.POINTER(ctypes.c_char_p), ctypes.c_ulong]

        retval = self.lib_EventManager.register_named_event(('_AOC_GENERAL_TRIGGER').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_PUPMASK_AVAILABLE').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_ZER_CAL_READY').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_ZON_CAL_READY').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_DM_STATIC_MAP_UPDATED').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_REF_CENTROIDS_SET').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_OK').encode('utf-8'),1)
        retval = self.lib_EventManager.register_named_event(('_AOC_CALIB_DONE').encode('utf-8'),1)

        super(MyWindow, self).__init__()
        if not os.path.exists(conf_dir + guidef):
            uic.loadUi(guidef, self)
        else:
            uic.loadUi(conf_dir + guidef, self)

        # ==============================================
        #         prepare the image displays
        # ==============================================
        # -- first, the camera --
        self.gView_SH_cam.hideAxis('left')
        self.gView_SH_cam.hideAxis('bottom')

        self.imv_raw = pg.ImageItem()
        self.overlayGrid = pg.GraphItem()
        self.overlayRefs = pg.GraphItem()
        self.overlayCentroids = pg.GraphItem()
        self.overlaySlopes = pg.GraphItem()
        self.overlayPupil = pg.GraphItem()
        # self.scatter = pg.ScatterPlotItem()


        self.gView_SH_cam.addItem(self.imv_raw)
        self.gView_SH_cam.addItem(self.overlayGrid)
        self.gView_SH_cam.addItem(self.overlayRefs)
        self.gView_SH_cam.addItem(self.overlayCentroids)
        self.gView_SH_cam.addItem(self.overlayPupil)
        # self.gView_SH_cam.GraphItem.sigMouseClicked.connect(self.OnPupilSelectorClick)

        # -- second, the WFS data --
        self.gView_SH_info.hideAxis('left')
        self.gView_SH_info.hideAxis('bottom')

        self.imv_WFS = pg.ImageItem()
        self.gView_SH_info.addItem(self.imv_WFS)

        self.gView_Slopes.hideAxis('left')
        self.gView_Slopes.hideAxis('bottom')
        self.imv_Slopes = pg.ImageItem()
        self.gView_Slopes.addItem(self.imv_Slopes)


        # -- the tip-tilt log plot --
        self.logplotx = self.gView_TT_log.plot([0,200],[0,0],
                                               pen=(0,255,0), name="ttx")
        self.logploty = self.gView_TT_log.plot([0,200],[0,0],
                                               pen=(0,0,255), name="tty")

        self.nmodes = 20
        self.bg_vx = 0.5 + np.arange(20)
        self.yA = np.linspace(0.0, 3.0, num=20)
        self.yB = np.linspace(0.0, 3.0, num=20)
        self.mode_graphA = pg.BarGraphItem(x=self.bg_vx, height=self.yA, width=0.3, brush='b')
        self.mode_graphB = pg.BarGraphItem(x=self.bg_vx+0.3, height=self.yB, width=0.3, brush='r')
        self.gView_ModePlot.addItem(self.mode_graphA)
        self.gView_ModePlot.addItem(self.mode_graphB)


        # ==============================================
        #             GUI widget actions
        # ==============================================
        self.dspB_disp_min.valueChanged[float].connect(self.update_vmin)
        self.dspB_disp_max.valueChanged[float].connect(self.update_vmax)

        self.chB_min.stateChanged[int].connect(self.update_vmin)
        self.chB_min.stateChanged[int].connect(self.update_vmax)
        self.chB_AutoSetWFSThresh.stateChanged[int].connect(self.update_AutoSetWFSThresh)
        self.chB_min.setChecked(True)
        self.chB_max.setChecked(True)
        self.chB_AutoSetWFSThresh.setChecked(True)
        self.update_vmin()

        self.chB_nonlinear.stateChanged[int].connect(self.update_nonlinear)

        self.cmB_cbar.addItems(['gray',    'hot',   'jet',
                                'viridis', 'magma', 'inferno',
                                'plasma'])
        self.cmB_cbar.activated[str].connect(self.update_cbar)
        self.cmB_cbar.setCurrentIndex(3)
        self.update_cbar()

        self.chB_show_grid.stateChanged[int].connect(self.redraw_SH_grid)
        self.check_ShowCentroids.stateChanged[int].connect(self.redraw_SH_grid)
        self.check_ShowRefPos.stateChanged[int].connect(self.redraw_SH_grid)
        self.pB_ValidateLeakageGain.clicked.connect(self.validate_leakage_gain)

        self.pB_reset_TCS_offload.clicked.connect(self.reset_TCS_offload)

        self.spB_grid_x0.valueChanged[int].connect(self.redraw_SH_grid)
        self.spB_grid_y0.valueChanged[int].connect(self.redraw_SH_grid)
        self.spB_grid_dx.valueChanged[float].connect(self.redraw_SH_grid)
        self.spB_grid_dy.valueChanged[float].connect(self.redraw_SH_grid)
        #self.spB_grid_imin.valueChanged[int].connect(self.updt_SH_threshold)

        self.chB_LogSlopesIllum.stateChanged[int].connect(self.set_slope_logging)
        self.chB_LogModes.stateChanged[int].connect(self.set_mode_logging)
        self.spB_SlopeLogNumSamples.valueChanged[int].connect(self.update_slope_log_numsamp)

        # ==============================================
        #             top-menu actions
        # ==============================================
        self.actionLoad_configuration.triggered.connect(self.load_config)
        self.actionSave_config_to.triggered.connect(self.save_config_to_file)
        self.actionSave_config.triggered.connect(self.save_config)
        self.actionQuit.triggered.connect(self.exit)
        self.actionQuit.setShortcut('Ctrl+Q')

        if self.simu:
            self.shmf = "R:\\ciao_shcam.im.shm"
            self.dm_map = 'R:\\dmdisp.wf.shm'
            self.spB_nav.setValue(10)
        else:
            self.shmf = "R:\\ixon.im.shm"
            self.spB_nav.setValue(50)

        self.shm_cam = shm(self.shmf,packed=False,verbose=False)
        self.wfs = WFS(shmf=self.shmf, ncpu=ncpu)

        self.shm_slopes = shm("R:\\slope_history.im.shm",packed=False,verbose=False)
        self.shm_refcent = shm("R:\\reference_centroids.im.shm",packed=False,verbose=False)
        self.shm_slope_mask = shm("R:\\pupil_slope_mask.im.shm",packed=False,verbose=False)
        self.shm_modes = shm("R:\\modes_info.im.shm",packed=False,verbose=False)
        self.shm_LOmodes = shm('R:\\DM_Buffer_02.im.shm', verbose=False)
        self.shm_HOmodes = shm('R:\\DM_Buffer_03.im.shm', verbose=False)
        self.shm_timestamps = shm('R:\\timestamps.im.shm', verbose=False)


        self.pB_WFS_off.clicked.connect(self.wfs_stop)
        self.pB_WFS_on.clicked.connect(self.wfs_start)
        self.pB_WFS_ref.clicked.connect(self.wfs_set_ref)
        self.pB_WFS_clear_ref.clicked.connect(self.wfs_clear_ref)
        self.pB_grid_set.clicked.connect(self.wfs_regrid)

        self.pB_DevButtonA.clicked.connect(self.DevButtonA)
        self.pB_DevButtonB.clicked.connect(self.DevButtonB)
        self.check_LoopUpdateStaticMap.stateChanged[int].connect(self.DM_LoopStaticMapUpdate)
        self.check_EnableDMcmd.stateChanged[int].connect(self.DM_EnableDMcmd)

        # --- camera / DM setup control ---
        self.pB_DM_shutdown.clicked.connect(self.alpao_shutdown)
        self.pB_DM_start.clicked.connect(self.alpao_start)

        #self.mycam = Cam(fifo_dir="/home/ciaodev/bin/")
        self.mycam = Cam(fifo_dir="\\\\.\\pipe\\")
        if self.mycam.connected:
            self.enable_cam_gui(state=True)
            self.dspB_CamExpTime.setValue(self.mycam.cam_tint*1000)

        # self.pB_cam_stream.clicked.connect(self.mycam.stream)
        # self.pB_cam_stop.clicked.connect(self.mycam.pause)
        self.pB_cam_stream.clicked.connect(self.streamCam)
        self.pB_cam_stop.clicked.connect(self.pauseCam)

        self.pB_SetTempSetPoint.clicked.connect(self.SetCamSetPoint)
        self.tpB_Cooling.clicked.connect(self.ToggleCamCooling)
        self.pB_QueryCamTemperature.clicked.connect(self.QueryCamTemperature)
        self.tED_CurrentTemperature.setPlainText("N/A")
        self.combo_CameraMode.currentIndexChanged.connect(self.comboChangeCameraMode)
        self.combo_CameraMode.addItem("NORMAL")
        self.combo_CameraMode.addItem("FT")
        self.combo_CameraMode.setCurrentIndex(0)

        # self.pB_cam_tint_inc.clicked.connect(self.mycam.tint_inc)
        # self.pB_cam_tint_dec.clicked.connect(self.mycam.tint_dec)
        self.dspB_CamExpTime.valueChanged[float].connect(self.update_exp_time)
        # self.shm_cam.read_keyword(0)
        # exptime_musec = self.shm_cam.kwds[0]['value']
        # print('exptime_musec >> ',self.shm_cam.kwds[0])

        self.pB_cam_shutdown.clicked.connect(self.ixon_shutdown)

        # self.pB_Simulation.clicked.connect(self.set_simulation)
        # self.pB_StopSimulation.clicked.connect(self.stop_simulation)

        self.chB_Simulation.stateChanged[int].connect(self.set_simulation)
        self.pB_selectDirectory.clicked.connect(self.select_directory)

        self.pB_changePrefix.clicked.connect(self.change_prefix)

        self.chB_TCSOffload.stateChanged[int].connect(self.TCSOffload)
        self.check_InvertX.stateChanged[int].connect(self.OffloadInvertX)
        self.check_InvertY.stateChanged[int].connect(self.OffloadInvertY)

        self.chB_ShowZerWindow.stateChanged[int].connect(self.show_zer_window)

        self.bTCSOffloadInvertX = True
        self.bTCSOffloadInvertY = False
        self.check_InvertX.setChecked(self.bTCSOffloadInvertX)
        self.check_InvertY.setChecked(self.bTCSOffloadInvertY)

        # ==============================================
        self.show()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.refresh_all)
        self.timer.start(50)

        # --- wavefront control callbacks ---
        self.bDoMeasure = False
        self.bCloseLoop = False


        self.TT_WFC  = TT_WFC(simu=self.simu)               # tip-tilt wavefront control
        self.ZON_WFC = ZON_WFC(simu=self.simu)              # Zonal wavefront contrl
        # self.ZER_WFC = ZER_WFC(iz0=4, iz1=10) # Zernike wavefront control
        self.ZER_WFC = ZER_WFC(iz0=2, iz1=10, simu=self.simu) # Zernike wavefront control
        self.ZER_WFC.gain = 0.25
        self.dspB_zer_gain.setValue(self.ZER_WFC.gain)
        self.ZER_norm_calib_modes = np.ones(self.ZER_WFC.modes.shape[0])

        print('Number of ZER_WFC instances: %d' % sys.getrefcount(self.ZER_WFC))

        self.TT_WFC.reload_cal(fname=conf_dir+"TT_CAL.fits")
        self.ZON_WFC.reload_cal(fname=conf_dir+"ZON_CAL.fits")
        self.ZER_WFC.reload_cal(fname=conf_dir+"ZER_CAL.fits")

        self.ZON_calAmp = self.dspB_zo_a0.value()
        self.ZON_WFC.a0 = self.ZON_calAmp
        # self.dspB_zo_a0.setValue(self.ZON_WFC.a0)

        self.ZER_calAmp = self.dspB_zer_a0.value()
        self.ZER_WFC.a0 = self.ZER_calAmp

        self.dspB_tt_a0.setValue(self.TT_WFC.a0)


        self.dspB_zer_nz.valueChanged[int].connect(self.zer_update_nz)

        self.ZER_WFC.verbose = True
        self.ZER_WFC.nmodes = self.dspB_zer_nz.value() - 1
        self.ZER_calib_done = False

        self.ZON_WFC.nmodes = 97
        self.dms = 11       # NUMBER OF ACTUATORS ACROSS (SQUARE)

        self.pB_tt_abort.clicked.connect(self.tt_abort)
        self.pB_tt_cloop.clicked.connect(self.tt_cloop)
        self.pB_tt_recal.clicked.connect(self.tt_calib)
        self.pB_tt_reset.clicked.connect(self.tt_reset_correc)
        self.dspB_tt_gain.valueChanged[float].connect(self.tt_update_gain)
        self.dspB_OffloadPixScale.valueChanged[float].connect(self.update_OffloadPixScale)

        self.pB_zo_abort.clicked.connect(self.zo_abort)
        self.pB_zo_measure.clicked.connect(self.zo_measure)
        self.pB_zo_recal.clicked.connect(self.zo_calib)
        self.pB_zo_reset.clicked.connect(self.zo_reset_correc)
        self.dspB_zo_gain.valueChanged[float].connect(self.zo_update_gain)
        self.chB_zo_cloop.stateChanged[int].connect(self.zo_cloop)

        self.pB_zer_abort.clicked.connect(self.zer_abort)
        self.pB_zer_measure.clicked.connect(self.zer_measure)
        self.pB_zer_recal.clicked.connect(self.zer_calib)
        self.pB_zer_reset.clicked.connect(self.zer_reset_correc)
        self.dspB_zer_gain.valueChanged[float].connect(self.zer_update_gain)
        self.chB_zer_cloop.stateChanged[int].connect(self.zer_cloop)
        self.log_dir = ""

        # ====================================================
        # chart
        # ====================================================

        self.new_qc = QChart()
        self.new_qcv = QChartView(self.new_qc)
        self.new_qbs = QBarSeries()
        self.bar_set = QBarSet("p_value")
        self.new_qc.setMargins(QtCore.QMargins(-5,-5,-5,-10))
        self.new_qc.setContentsMargins(QtCore.QMarginsF(0, 0, 0, 0))
        self.new_qc.setBackgroundRoundness(0)
        self.new_qc.legend().hide()
        self.histo = []
        for i in range(1, 164):
            self.histo.append(i)
            self.bar_set.append(self.histo[i - 1])
        self.new_qbs.append(self.bar_set)
        self.new_qbs.setBarWidth(1)
        self.new_qc.addSeries(self.new_qbs)
        self.axis = QValueAxis()
        self.axis.setTickCount(10)
        self.axis.setLabelFormat("%d")
        self.new_qc.createDefaultAxes()
        self.new_qc.setAxisX(self.axis, self.new_qbs)
        self.axisY = QValueAxis()
        self.axisY.setRange(0, 2500)
        self.new_qc.setAxisY(self.axisY, self.new_qbs)
        self.vL_Chart.addWidget(self.new_qcv)

        self.qc = QChart()
        self.qcv = QChartView(self.qc)
        self.qbs = QBarSeries()
        self.set0 = QBarSet("Value")
        self.set0 << 0 << 100 << 200 << 300 << 400 << 500 << 600 << 700 << 800 << 900
        self.qbs.append(self.set0)
        self.qbs.setBarWidth(1)
        #self.qbs.setBarSpacing(0)
        self.qc.addSeries(self.qbs)
        #self.vL_chart.addWidget(self.qcv)

        self.axis = QValueAxis()
        self.axis.setTickCount(10)
        self.axis.setLabelFormat("%d")
        self.qc.createDefaultAxes()
        self.qc.setAxisX(self.axis, self.qbs)
        self.axisY = QValueAxis()
        self.axisY.setRange(0, 2500)
        self.qc.setAxisY(self.axisY, self.qbs)
        #self.vL_chart.setStyleSheet('background-color:grey')

        self.chart_win = None
        self.chB_showChart.stateChanged[int].connect(self.show_chart)

        # ====================================================
        # live tip-tilt offsets
        # ====================================================
        self.pB_tt_dx_p.clicked.connect(self.ttx_offset_p)
        self.pB_tt_dx_m.clicked.connect(self.ttx_offset_m)
        self.pB_tt_dy_p.clicked.connect(self.tty_offset_p)
        self.pB_tt_dy_m.clicked.connect(self.tty_offset_m)
        # ====================================================

        aoc_params = str(self.QueryAOCParams())
        aoc_params = aoc_params.replace("b'", "")
        aoc_params = aoc_params.replace("'", "")
        list_params = aoc_params.split(';')
        print(list_params)

        isSim = False

        if list_params[3] == "1" and list_params[4] == "1":
            print("Simulation mode")
            isSim = True
            self.chB_Simulation.setChecked(True)


        if isSim is False:
            print('Starting telescope offload thread.')
            self.tel = Tel()
            self.tel_mon_start()
        # self.tel_derot_offset = self.dspB_DerotOffset.value()
        # self.tel.offset_angle_deg = 75.
        self.derotator_angle_deg = 0.
        self.derotator_on = 0

        self.zer_win = None
        self.pup_win = None

        self.zer_update_nz()

        self.ttx_drift_avg = 0.
        self.tty_drift_avg = 0.
        self.vHO_global_amp = np.zeros(self.ZER_WFC.nmodes - 2)

        # ------------------
        #   Event monitoring thread
        # ------------------
        self.thread_event_monitor = threading.Thread(target=self.thread_eventmon_start, args=())
        self.thread_event_monitor.start()
        # self.thread_eventmon_start()# start the event monitoring thread

        # INITIALIZE DEFORMABLE MIRROR
        self.mycam.send_fifo_cmd(b'alpao init_dm')

        # MQTT CLIENT
        if isSim is False:
            print('Registering MQQT server...')
            username = 'neo'
            pwd = 'lagrange'
            client = mqtt.Client("JOVIAL")
            client.on_connect = self.on_connect
            client.on_message = self.on_message
            client.username_pw_set(username, password=pwd)
            # client.connect('10.150.1.100',port=1883)
            client.loop_start()


        self.check_ShowPupil.stateChanged[int].connect(self.show_pupil)
        self.win_pupil = 0

        # self.mycam.send_fifo_cmd(b'alpao load_flat \"C:\\Users\\AOC\\bin\\alpao_flat.fits\"')
        self.mycam.send_fifo_cmd(b'alpao load_flat \"C:\\Users\\AOC\\bin\\zygo_flat_97.fits\"')

        self.mycam.send_fifo_cmd(b'alpao load_map \"C:\\Users\\AOC\\bin\\flat_map_correction_97.fits\" 5')

        self.nSubAps = (self.dms-1)*(self.dms-1)
        print('self.nSubAps >> ', self.nSubAps)

        self.vCurSubApSel = np.zeros(self.nSubAps).astype(np.bool8)
        self.vSubApThresh = np.full(self.nSubAps, self.dspB_disp_min.value()).astype(np.float32)
        self.vSubApEnable = np.ones(self.nSubAps).astype(np.int32)

        self.ref_centroids_data = np.zeros((self.dms-1)*(self.dms-1)*2).astype(np.float32)
        self.ref_centroids_offsets = np.zeros((self.dms-1)*(self.dms-1)*2).astype(np.float32)

        # TRY TO LOAD PREVIOUS PUPIL SUBAPERTURES THRESHOLD/ENABLED FILE
        self.subaps_threshold_file = 'C:\\Users\\lucas\\Documents\\STAGE\\CIAO_Win\\data\\vSubApsThresholds_last.fits'
        self.subaps_enabled_file = 'C:\\Users\\lucas\\Documents\\STAGE\\CIAO_Win\\data\\vSubApsEnabled_last.fits'
        self.bSubApsThresholdsFileExits = os.path.exists(self.subaps_threshold_file)
        self.bSubApsEnabledFileExits = os.path.exists(self.subaps_enabled_file)

        if self.bSubApsThresholdsFileExits == True:
            print('=> Found previous subaperture threshold file and loading it.')
            self.vSubApThresh = np.array(fits.getdata(self.subaps_threshold_file)).astype(np.float32)
        else:
            print('=> No previous subaperture threshold file found (using default values).')

        if self.bSubApsEnabledFileExits == True:
            print('=> Found previous enabled subaperture file and loading it.')
            self.vSubApEnable = np.array(fits.getdata(self.subaps_enabled_file)).astype(np.int32)
        else:
            print('=> No previous enabled subaperture file found (using default values).')

        # WRITE THE FILE FOR THE IXON SERVER TO KNOW ABOUT THE THRESHOLDS
        fits.writeto('R:\\ixon_subaperture_thresholds.fits',self.vSubApThresh,overwrite=True)
        fits.writeto('R:\\ixon_subaperture_enabled.fits',self.vSubApEnable,overwrite=True)


        self.mycam.send_fifo_cmd(b'ixon set_thresh \"R:\\ixon_subaperture_thresholds.fits\" \"R:\\ixon_subaperture_enabled.fits\"')

        self.ref_centroids_data = self.shm_refcent.get_data(check=False, reform=True, sleepT=0.).flatten()
        # self.ref_centroids = self.shm_refcent.get_data(check=False, reform=True, sleepT=0.).flatten()
        self.ref_centroids = self.ref_centroids_data

        time.sleep(1.)
        print('COMMAND')
        self.lib_EventManager.trigger_named_event(('_AOC_DM_STATIC_MAP_UPDATED').encode('utf-8'))

        self.redraw_SH_grid()
        self.define_default_SH_pupil()

        self.DM_EnableDMcmd()

        # self.zoffsets_win = zoffsets_window(self)
        # self.zoffsets_win.show()

        self.bInitDone = True

        self.read_imagettes()

    def streamCam(self):
        self.combo_CameraMode.setEnabled(False)
        self.mycam.stream()

    def pauseCam(self):
        self.combo_CameraMode.setEnabled(True)
        self.mycam.pause()

    def on_connect(self, client, userdata, flags, rc):
        print('Connected with result code {0}'.format(str(rc)))
        client.subscribe("JOVIAL/derotator")

    def on_message(self, client, userdata, msg):
        self.derotator_angle_deg = float(json.loads(msg.payload)['theta_c'])
        self.derotator_on = float(json.loads(msg.payload)['on'])
        # self.tel.offset_angle_deg = 75. + self.derotator_angle_deg*self.derotator_on
        # self.tel.offset_angle_deg = -self.dspB_DerotOffset.value() + self.derotator_angle_deg*self.derotator_on
        # print('MQTT com happened...', self.dspB_DerotOffset.value(), self.derotator_angle_deg, self.derotator_on)


        # print('Message received >> ' + msg.topic + ' ' + str(msg.payload), str(self.derotator_angle_deg))
        # print('Angle sent to tel >> ' ,str(self.tel.offset_angle_deg))

    def show_chart(self):
        if self.chart_win == None:
            self.chart_win = ChartWindow()

            self.chart_win.vl_Chart.addWidget(self.qcv)
        if self.chB_showChart.isChecked() == True:
            self.chart_win.show()
            self.chart_win.parentWindow = self

        else:
            self.chart_win.hide()


    def show_pupil(self, checked):
        # print('"show_pupil" function unimplemented')
        # self.lib_EventManager.trigger_named_event(('_AOC_OK').encode('utf-8'))
        self.vCurCellThreshold = np.array(100).astype(np.float32)
        self.vCurEnabled = np.array(100).astype(np.bool8)
        print(self.pupil_mask.flatten())
        if self.pup_win == None:
            # self.pup_win = PupilWindow(self,pupmask=self.pupil_mask.flatten(),vCurCellThreshold=self.vCurCellThreshold,vCurEnabled=self.vCurEnabled)
            # self.pup_win = PupilWindow(self,pupmask=self.pupil_mask.flatten(),vCurCellThreshold=self.vCurCellThreshold,vCurEnabled=self.vCurEnabled)
            self.pup_win = PupilWindow(self)
            self.pup_win.pupmask=self.pupil_mask.flatten()
            # self.vCurCellThreshold=self.vCurCellThreshold,vCurEnabled=self.vCurEnabled)
            self.pup_win.parentWindow = self

        if self.check_ShowPupil.isChecked() == True:
            self.pup_win.show()
            self.bShowCellSelectionFeedback = True
        else:
            self.pup_win.hide()
            self.bShowCellSelectionFeedback = False

    def define_default_SH_pupil(self):
        # print('"define_default_SH_pupil" function unimplemented')
        nx = 10
        ny = 10
        self.pupil_mask = np.ones((nx, ny), dtype=np.int32) # corner mask

        # self.dmmask[np.abs(self.xdm) + np.abs(self.ydm) > 7] = 0.0

        self.pupil_mask[0,0:3] = 0
        self.pupil_mask[0,ny-3:] = 0
        self.pupil_mask[nx-1,0:3] = 0
        self.pupil_mask[nx-1,ny-3:] = 0

        self.pupil_mask[1,0:2] = 0
        self.pupil_mask[1,ny-2:] = 0
        self.pupil_mask[nx-2,ny-2:] = 0
        self.pupil_mask[nx-2,0:2] = 0

        self.pupil_mask[2,0] = 0
        self.pupil_mask[2,ny-1] = 0
        self.pupil_mask[nx-3,0] = 0
        self.pupil_mask[nx-3,ny-1] = 0

        # print('self.pupil_mask.size >> ',self.pupil_mask.size)

        # pen = QtGui.QPen()
        # pen.setWidth(0)
        # pen.setColor(QtGui.QColor(0, 255, 0, 127))

        # hoverPen = QtGui.QPen()
        # hoverPen.setWidth(0)
        # hoverPen.setColor(QtGui.QColor(255, 127, 0, 127))

        # brush = QtGui.QBrush()
        # brush.setColor(QtGui.QColor(0, 255, 0, 127))
        # brush.setStyle(Qt.SolidPattern)

        # pSize = 20

        # count = 0
        # for j in range(10):
        #     for i in range(10):
        #         if self.pupil_mask[i,j] == 1:
        #             # painter.drawRect(self.pos_xx[2*i, 0], self.pos_yy[2*j, 1], self.pos_xx[2*i+2, 0], self.pos_yy[2*j+2, 1])
        #             x1 = self.pos_xx[2*i, 0]
        #             y1 = self.pos_yy[2*j, 1]
        #             x2 = self.pos_xx[2*i+2, 0]
        #             y2 = self.pos_yy[2*j+2, 1]

        #             vpos = [0.5*(x1+x2),0.5*(y1+y2)]

        #             if count > 0:
        #                 self.vItemValidCells = self.vItemValidCells + [pg.ScatterPlotItem([vpos[0]], [vpos[1]], size=pSize, pen=pen, brush=brush,
        #                     symbol='s', hoverable=False, hoverPen=hoverPen, tip=False, pxMode=1, Name='subap'+str(count))]
        #             else:
        #                 self.vItemValidCells = [pg.ScatterPlotItem([vpos[0]], [vpos[1]], size=pSize, pen=pen, brush=brush,
        #                     symbol='s', hoverable=False, hoverPen=hoverPen, tip=False, pxMode=1, Name='subap' + str(count))]

        #             self.gView_SH_cam.addItem(self.vItemValidCells[count])

        #             # self.scatter.setData(vpos[0],vpos[1],pen=pen,brush=brush,size=10,pxMode=False)
        #             # print('Valid cell [i,j]:   ',i,j, self.pos_xx[2*i, 0], self.pos_yy[2*j, 1], self.pos_xx[2*i+2, 0], self.pos_yy[2*j+2, 1])
        #             count = count + 1


    def OnPupilSelectorClick(self,item,point):
        print('Clicked item name >> ', item.Name)

    def show_zer_window(self, checked):
        if self.zer_win == None:
            self.zer_win = ZernikeWindow()
            self.zer_win.parentWindow = self

        if self.chB_ShowZerWindow.isChecked() == True:
            self.zer_win.show()
        else:
            self.zer_win.hide()

    # =========================================================
    def load_config(self):
        filename = QFileDialog.getOpenFileName(self, 'Open File for loading', 'C:\\RSI\\SourcesC\\CIAO_win\\control\\','AOC cfg file (*.aoc)')
        print('Selected file:')
        print(filename[0])
        if filename != "":
            self.curConfigFile = filename[0]

        config = configparser.ConfigParser()
        content = config.read(self.curConfigFile)
        print(content)

        dx = float(config['WFSGrid']['DX'])
        dy = float(config['WFSGrid']['DY'])
        x0 = int(config['WFSGrid']['X0'])
        y0 = int(config['WFSGrid']['Y0'])
        vmin = int(config['WFSGrid']['DispMin'])
        vmax = int(config['WFSGrid']['DispMax'])
        nav = int(config['Calibration']['NumSampPerCalib'])
        self.OffloadPixScale = float(config['TCS_Offload']['OffloadPixScale'])
        self.tel_mon_period_sec = float(config['TCS_Offload']['tel_mon_period_sec'])

        print('dx: %.2f  dy: %.2f  x0 :%d  y0: %d' % (dx,dy,x0,y0))

        self.spB_grid_dx.setValue(dx)
        self.spB_grid_dy.setValue(dy)
        self.spB_grid_x0.setValue(x0)
        self.spB_grid_y0.setValue(y0)
        self.self.dspB_disp_min.setValue(vmin)
        self.self.dspB_disp_max.setValue(vmax)
        self.spB_nav.setValue(nav)
        self.dspB_OffloadPixScale.setValue(self.OffloadPixScale)

    def CfgSave(self):
        print('Saving configuration to "%s".' % ( self.curConfigFile ))

        # Create the configuration file as it doesn't exist yet
        cfgfile = open(self.curConfigFile, "w")

        # Add content to the file
        Config = configparser.ConfigParser()
        Config.add_section("WFSGrid")
        Config.set("WFSGrid", "DX", str(self.spB_grid_dx.value()))
        Config.set("WFSGrid", "DY", str(self.spB_grid_dy.value()))
        Config.set("WFSGrid", "X0", str(self.spB_grid_x0.value()))
        Config.set("WFSGrid", "Y0", str(self.spB_grid_y0.value()))
        # MIN/MAX THRESHOLDS FOR WFS IMAGE DISPLAY
        Config.set("WFSGrid", "DispMin", str(self.dspB_disp_min.value()))
        Config.set("WFSGrid", "DispMax", str(self.dspB_disp_max.value()))

        Config.add_section("Calibration")
        Config.set("Calibration", "NumSampPerCalib", str(self.spB_nav.value()))
        Config.add_section("TCS_Offload")
        Config.set("TCS_Offload", "OffloadPixScale", str(self.dspB_OffloadPixScale.value()))
        Config.set("TCS_Offload", "tel_mon_period_sec", str(self.tel_mon_period_sec))

        Config.write(cfgfile)
        cfgfile.close()


    # =========================================================
    def save_config_to_file(self):
        filename = QFileDialog.getSaveFileName(self, 'Open File for saving', 'C:\\RSI\\SourcesC\\CIAO_win\\control\\','AOC cfg file (*.aoc)')
        # print('Selected file:')
        # print(filename[0])
        if filename[0] != "":
            self.curConfigFile = filename[0]
            self.bSaveCfgFileDefined = True
        else:
            print('No file specified. Nothing done.')
            return

        self.CfgSave()

    # =========================================================
    def save_config(self):
        if self.bSaveCfgFileDefined == False:
            filename = QFileDialog.getSaveFileName(self, 'Open File for saving', 'C:\\RSI\\SourcesC\\CIAO_win\\control\\','AOC cfg file (*.aoc)')
            # print('Selected file:')
            # print(filename[0])
            if filename[0] != "":
                self.curConfigFile = filename[0]
                self.bSaveCfgFileDefined = True
            else:
                print('No file specified. Nothing done.')
                return

        print('Saving configuration to "%s".' % ( self.curConfigFile ))

        self.CfgSave()

        # # Add content to the file
        # Config = configparser.ConfigParser()
        # Config.add_section("WFSGrid")
        # Config.set("WFSGrid", "DX", str(self.spB_grid_dx.value()))
        # Config.set("WFSGrid", "DY", str(self.spB_grid_dy.value()))
        # Config.set("WFSGrid", "X0", str(self.spB_grid_x0.value()))
        # Config.set("WFSGrid", "Y0", str(self.spB_grid_y0.value()))
        # Config.write(cfgfile)
        # cfgfile.close()

    # =========================================================
    def set_slope_logging(self,):
        self.TT_WFC.bLogSlopes = self.chB_LogSlopesIllum.isChecked()

        if self.TT_WFC.bLogSlopes == True:
            self.log("Slope+Illum (TT) data logging to file enabled", "green")
            self.mycam.send_fifo_cmd(b'system log INIT')
            time.sleep(1.)
            self.mycam.send_fifo_cmd(b'system log ON')
            if self.log_dir == "":
                os.write(self.mycam.cmd_fifo, b'system query_log_dir')
                os.fsync(self.mycam.cmd_fifo)
                tmp = os.read(self.mycam.cmd_fifo, self.mycam.MAX_FIFO_LEN)
                os.fsync(self.mycam.cmd_fifo)
                print(tmp)
                self.log_dir = str(tmp) + 'calib\\'
                self.log_dir = self.log_dir.replace("b'", "")
                self.log_dir = self.log_dir.replace("'", "")
                fits.writeto(self.log_dir + 'ixon_subaperture_thresholds.fits', self.vSubApThresh, overwrite=True)
                fits.writeto(self.log_dir + 'ixon_subaperture_enabled.fits', self.vSubApEnable, overwrite=True)
                fits.writeto(self.log_dir + 'ActuatorsMappingIndices.fits', self.vSubApEnable, overwrite=True)
                shutil.copy('R:\\default_ref_pos_x.fits', self.log_dir)
                shutil.copy('R:\\default_ref_pos_y.fits', self.log_dir)
        else:
            self.log("Slope+Illum (TT) data logging to file disabled", "red")
            self.mycam.send_fifo_cmd(b'system log OFF')

    # =========================================================

    def set_simulation(self):
        if self.chB_Simulation.isChecked():
            self.mycam.send_fifo_cmd(b'system simulation ON')
            print("system simulation on")
        else:
            self.mycam.send_fifo_cmd(b'system simulation OFF')
            print("system simulation off")

    def select_directory(self):
        directory = QFileDialog()
        directory.setFileMode(QFileDialog.DirectoryOnly)
        directory_path = directory.getExistingDirectoryUrl(self, 'C:\\Users\\lucas\\Documents\\STAGE\\FITS\\')
        if directory_path:
            print(directory_path.path())
            path = directory_path.path()[1:]
            b = bytearray()
            b.extend(map(ord, path))
            self.mycam.send_fifo_cmd(b'system change_directory %b' % b)

    def change_prefix(self):
        text, ok = QInputDialog().getText(self, "QInputDialog().getText()",
                                          "Object Name:", QLineEdit.Normal,
                                          "obs")
        if ok and text:
            print(text)
            b = bytearray()
            b.extend(map(ord, text))
            self.mycam.send_fifo_cmd(b'system change_prefix %b' % b)

    def stop_simulation(self):
        self.mycam.send_fifo_cmd(b'system simulation OFF')
        print("system simulation off")

    # =========================================================
    def set_mode_logging(self,):
        self.ZER_WFC.bLogCmds = self.chB_LogModes.isChecked()
        self.ZER_WFC.bLogErrors = self.ZER_WFC.bLogCmds
        self.ZER_WFC.bLogSlopes = self.ZER_WFC.bLogCmds

        if self.ZER_WFC.bLogCmds == True:
            self.log("Modes (ZER) data logging to file enabled", "green")
            self.mycam.send_fifo_cmd(b'system log ON')
            if self.log_dir == "":
                os.write(self.mycam.cmd_fifo, b'system query_log_dir')
                os.fsync(self.mycam.cmd_fifo)
                tmp = os.read(self.mycam.cmd_fifo, self.mycam.MAX_FIFO_LEN)
                os.fsync(self.mycam.cmd_fifo)
                print(tmp)
                self.log_dir = str(tmp) + 'calib\\'
                self.log_dir = self.log_dir.replace("b'", "")
                self.log_dir = self.log_dir.replace("'", "")
                shutil.copy('R:\\default_ref_pos_x.fits', self.log_dir)
                shutil.copy('R:\\default_ref_pos_y.fits', self.log_dir)
        else:
            self.log("Modes (ZER) data logging to file disabled", "red")


    # =========================================================
    def update_slope_log_numsamp(self,):
        self.TT_WFC.log_length = self.spB_SlopeLogNumSamples.value()
        self.ZER_WFC.log_length = self.TT_WFC.log_length

    def SetCamSetPoint(self,):
        temp_goal = int(self.dspB_CamTemperature.value())
        print('Setting iXon temperature goal to ',temp_goal)
        self.mycam.set_tempgoal(temp_goal)

    def ToggleCamCooling(self,):
        if self.tpB_Cooling.isChecked():
            self.mycam.set_cooling(1)
        else:
            self.mycam.set_cooling(0)

    def QueryCamTemperature(self,):
        fTemp = self.mycam.get_fTemp()
        if fTemp == -999:
            self.tED_CurrentTemperature.setPlainText('N/A')
        else:
            self.tED_CurrentTemperature.setPlainText("%.3f" % (fTemp))

    def QueryAOCParams(self,):
        if self.mycam.connected:
            os.write(self.mycam.cmd_fifo,b'system query_params')
            os.fsync(self.mycam.cmd_fifo)
            tmp = os.read(self.mycam.cmd_fifo, self.mycam.MAX_FIFO_LEN)
            os.fsync(self.mycam.cmd_fifo)
            print("test")
            print(tmp)
            return tmp

    def comboChangeCameraMode(self):
        self.CameraMode = self.combo_CameraMode.currentIndex()
        if self.CameraMode == 0:
            self.mycam.send_fifo_cmd(b'ixon set_transfer_mode 0')
            self.log("Camera transfer mode set to NORMAL", "green")
        elif self.CameraMode == 1:
            self.mycam.send_fifo_cmd(b'ixon set_transfer_mode 1')
            self.log("Camera transfer mode set to FRAME TRANSFER (FT)", "green")

    def update_exp_time(self,):
        tint = self.dspB_CamExpTime.value()
        self.mycam.set_tint(tint*1e-3)
        self.refresh_stats()
        self.log("Exposure time set to %.4e ms" % (tint,))

    # =========================================================
    def disable_wfc(self, subsystem=""):
        ''' ----------------------------------------------------------
        Disable *all* buttons associated to the WFC loops (default)

        unless a valid subsystem is specified in which case, only that
        one is disabled.
        ---------------------------------------------------------- '''
        wfcs = ["tt", "zo"]

        if subsystem in wfcs:
            wfcs = [subsystem]
        else:
            self.wfc_enabled = False # WFC disabled globally
            print("global wfc disable")
        # for wfc in wfcs:
        #     exec('self.pB_%s_recal.setEnabled(False)' % (wfc,))
        #     exec('self.chB_%s_cloop.setEnabled(False)' % (wfc,))
        # self.pB_tt_recal.setEnabled(False)
        # self.pB_tt_cloop.setEnabled(False)

        # self.pB_zo_recal.setEnabled(False)
        # self.pB_zo_measure.setEnabled(False)

        self.chB_zer_cloop.setEnabled(True)
        self.chB_zo_cloop.setEnabled(True)

    # =========================================================
    def enable_wfc(self, subsystem=""):
        ''' ----------------------------------------------------------
        Enable *all* buttons associated to the WFC loops

        unless a valid subsystem is specified in which case, only that
        one is enabled.
        ---------------------------------------------------------- '''
        wfcs = ["tt", "zo"] # list of valid wfc loops

        if subsystem in wfcs:
            wfcs = [subsystem]
            print("global wfc ENABLE!")
        else:
            self.wfc_enabled = True # WFC enabled globally

        # for wfc in wfcs:
            # exec('self.pB_%s_recal.setEnabled(True)' % (wfc,))
            # exec('self.chB_%s_cloop.setEnabled(True)' % (wfc,))
            # print("%s enabled" % (wfc,))
        self.pB_tt_recal.setEnabled(True)
        self.pB_tt_cloop.setEnabled(True)

        self.pB_zo_recal.setEnabled(True)
        self.pB_zo_measure.setEnabled(True)

        self.pB_zer_recal.setEnabled(True)
        self.pB_zer_measure.setEnabled(True)

        self.chB_zer_cloop.setEnabled(True)
        self.chB_zo_cloop.setEnabled(True)

    # =========================================================
    def wfs_regrid(self,):
        ''' ----------------------------------------------------------
        Updates the grid used by the WFS to output signal
        ---------------------------------------------------------- '''
        self.wfs.update_grid(
            x0=self.spB_grid_x0.value(),
            y0=self.spB_grid_y0.value(),
            dx=self.spB_grid_dx.value(),
            dy=self.spB_grid_dy.value(),
            i0=self.spB_grid_imin.value())
        self.wfs.update_cells()
        self.wfs.define_SH_data()
        self.log("Updated grid set-up", "blue")

    # =========================================================
    def enable_wfs_gui(self, state=True):
        ''' Enables/disables buttons of the WFS GUI'''
        self.pB_WFS_on.setEnabled(state)
        self.pB_grid_set.setEnabled(state)

    def DevButtonA(self):
        # self.mycam.send_fifo_cmd(b'alpao update_combined')
        print('self.last_slope_x >> ', self.last_slope_x.shape)
        print('self.last_slope_y >> ', self.last_slope_y.shape)
        slope_signal = np.concatenate((self.last_slope_x, self.last_slope_y), axis=None)
        print('slope_signal >> ', slope_signal.shape)
        ee  = self.RRinv.dot(slope_signal)    # error signal
        print('ee >>', self.RRinv.shape, ee.shape)
        print(ee)

    def DevButtonB(self):
        self.mycam.send_fifo_cmd(b'alpao update_static_maps')

    def DM_LoopStaticMapUpdate(self):
        state = self.check_LoopUpdateStaticMap.isChecked()
        if state == True:
            self.startThread_loop_update_DMstatic()
        else:
            self.DMStaticUpdate_loop_keepgoing = False

    def DM_EnableDMcmd(self):
        if self.check_EnableDMcmd.isChecked():
            self.mycam.send_fifo_cmd(b'alpao enable_dm 1')
        else:
            self.mycam.send_fifo_cmd(b'alpao enable_dm 0')

    # =========================================================
    def wfs_start(self,):
        ''' ----------------------------------------------------------
        Starts the process measuring WFS data (ttxy + phot)
        ---------------------------------------------------------- '''
        self.enable_wfs_gui(state=False)
        self.wfsThread = GenericThread(self.wfs.loop)
        self.wfsThread.start()
        self.log("Wavefront sensing resumed", "green")
        self.wfc_status_updated = True
        self.enable_wfc()

    # =========================================================
    def wfs_stop(self,):
        ''' ----------------------------------------------------------
        Interrupts the process measuring WFS data (ttxy + phot)
        ---------------------------------------------------------- '''
        self.wfs.keepgoing = False
        self.enable_wfs_gui(state=True)
        self.log("Wavefront sensor paused", "orange")
        self.disable_wfc()

    # =========================================================
    def wfs_set_ref(self,):
        ''' ----------------------------------------------------------
        Sets the current tip-tilt position as the reference state
        ---------------------------------------------------------- '''
        # self.wfs.set_thread_params()
        # self.wfs.keepgoing = False
        self.enable_wfs_gui(state=True)
        # self.ZER_WFC.DM.dmmask = self.wfs.SH_phot
        # fits.writeto(r'C:\Home\Downloads\self.ZER_WFC.DM.dmmask.fits',self.ZER_WFC.DM.dmmask,overwrite=True)
        # self.wfs.calc_SH_data(ref=True)

        # X AND Y CENTROIDS ARE STORED IN THE SAME ARRAY (ALL X FOLLOWED BY ALL Y VALUES)
        self.ref_centroids_data = np.zeros((self.dms-1)*(self.dms-1)*2).astype(np.float32)
        for k in range(50):
            self.ref_centroids_data[:] = self.ref_centroids_data[:] + self.data_slopes[k,0:200]
            # print('k, data',k,ref_centroids)

        self.ref_centroids_data[:] = self.ref_centroids_data[:] / 50.

        # ADD THE MODES OFFSETS
        self.ref_centroids[:] = self.ref_centroids_data[:] + self.ref_centroids_offsets[:]

        self.shm_refcent.set_data(self.ref_centroids)
        self.shm_refcent.save_as_fits('R:\\ref_cent_last.fits')
        if self.log_dir != "":
            shutil.copy('R:\\ref_cent_last.fits', self.log_dir)
        self.lib_EventManager.trigger_named_event(('_AOC_REF_CENTROIDS_SET').encode('utf-8'))

        self.log("Current pos. set as reference", "blue")

    def wfs_set_ref_offsets(self):

        # ADD THE MODES OFFSETS
        self.ref_centroids[:] = self.ref_centroids_data[:] + self.ref_centroids_offsets[:]

        self.shm_refcent.set_data(self.ref_centroids)
        self.shm_refcent.save_as_fits('R:\\ref_cent_last.fits')

        self.lib_EventManager.trigger_named_event(('_AOC_REF_CENTROIDS_SET').encode('utf-8'))

        self.log("Current pos. set as reference", "blue")

    def read_imagettes(self):
        for i in range(100):
            image = fits.getdata('C:\\Users\\lucas\\Documents\\imagettes_fixed\\imagettes' + str(i) + '.fits')
            xs = len(image[0])
            ys = len(image)
            xc = xs / 2
            yc = ys / 2
            sum_x_col = 0
            sum_y_row = 0
            sum_x = 0
            sum_y = 0
            sum_glo = 0

            for x in range(xs):
                sum_x_col = 0
                for y in range(ys):
                    sum_x_col += image[y][x]
                    sum_glo += image[y][x]
                sum_x += sum_x_col * x
           # print("x: ", sum_x / sum_glo, " ", (sum_x / sum_glo) - xc)
            sum_glo = 0
            for y in range(ys):
                sum_y_row = 0
                for x in range(xs):
                    sum_y_row += image[y][x]
                    sum_glo += image[y][x]
                sum_y += sum_y_row * y
            slope_x = (sum_x / sum_glo) - xc
            slope_y = (sum_y / sum_glo) - yc
            #print("y: ", sum_y / sum_glo, " ", (sum_y / sum_glo) - yc)
            #print("image ", i, " with sum_y ", sum_y_row)
            # fichier = open("C:\\Users\\lucas\\Documents\\STAGE\\Misc\\data_python_slopes.txt", "a")
            # to_write = str(slope_x) + " " + str(slope_y)
            # fichier.write("\n" + to_write)
            # fichier.close()



    def wfs_clear_ref(self,):
        ''' ----------------------------------------------------------
        RESET THE REFERENCE CENTROIDS TO ZERO
        ---------------------------------------------------------- '''
        self.ref_centroids = np.zeros((self.dms-1)*(self.dms-1)*2).astype(np.float32)
        self.ref_centroids[0:100] = np.array(fits.getdata("R:\\default_ref_pos_x.fits"))
        self.ref_centroids[100:200] = np.array(fits.getdata("R:\\default_ref_pos_y.fits"))

        self.shm_refcent.set_data(self.ref_centroids)

        self.lib_EventManager.trigger_named_event(('_AOC_REF_CENTROIDS_SET').encode('utf-8'))

        self.log("Reference centroids set to zero!", "blue")

    # =========================================================
    def exit(self):
        try:
            pf.writeto(conf_dir+"ZER_CAL.fits",
                       self.ZER_WFC.RR,
                       overwrite=True)
        except:
            print("Zernike calibration matrices not saved")

        try:
            pf.writeto(conf_dir+"TT_CAL.fits",
                       self.TT_WFC.RR,
                       overwrite=True)
        except:
            print("Tip-tilt calibration matrices not saved")

        try:
            pf.writeto(conf_dir+"ZON_CAL.fits",
                       self.ZON_WFC.RR,
                       overwrite=True)
        except:
            print("Zonal calibration matrices not saved")

        # UNREGISTER ALL EVENTS AND FREE DLL
        retval = self.lib_EventManager.unregister_all()
        _ctypes.FreeLibrary(self.lib_EventManager._handle)

        sys.exit()

    # =========================================================
    def closeEvent(self, event):
        ''' ----------------------------------------------------------
        GUI call: activated when corner of window is clicked.
        Attempts to save the calibration matrices!
        ---------------------------------------------------------- '''
        self.TT_WFC.save_cal(conf_dir+"TT_CAL.fits")
        self.ZER_WFC.save_cal(conf_dir+"ZER_CAL.fits")
        self.ZON_WFC.save_cal(conf_dir+"ZON_CAL.fits")

        self.wfs.keepgoing     = False
        self.wfs.stop() # explicit for release of CPUs
        self.TT_WFC.keepgoing  = False
        self.ZER_WFC.keeggoing = False
        self.ZON_WFC.keeggoing = False

        time.sleep(1) # pause to shutdown external processes
        sys.exit()

    # =========================================================
    def refresh_stats(self):
        ''' ----------------------------------------------------------
        Displays the properties of the live image + some other info,
        such as: exposure time & tip-tilt.
        ---------------------------------------------------------- '''

        pt_levels = [0, 5, 10, 20, 50, 75, 90, 95, 99, 100]
        pt_values = np.percentile(self.data_cam, pt_levels)

        msg = "<pre>\n"
        msg += "exp.time = %7.3f ms\n" % (self.mycam.cam_tint * 1e3)

        print(self.data_cam.size)
        for i in range(128):
            res = np.sum(self.data_cam[i])
            #print(res)

        for i, ptile in enumerate(pt_levels):
            # self.set0 << pt_values[i] << pt_values[i] << pt_values[i] << pt_values[i] << pt_values[i]\
            #     << pt_values[i] << pt_values[i] << pt_values[i] \
            #     << pt_values[i] << pt_values[i]
            self.set0.replace(i, pt_values[i])
            msg += "p-tile %3d = %8.2f\n" % (ptile, pt_values[i])

        msg += "\ntt_xy = (%+.2f,%+.2f)\n" % (self.wfs.ttx_mean, self.wfs.tty_mean)
        msg += "</pre>"
        self.lbl_stats.setText(msg)

    # =========================================================
    def redraw_SH_grid(self):
        ''' ----------------------------------------------------------
        Displays a grid overlay on the camera image.
        ---------------------------------------------------------- '''

        self.SH_x0 = self.spB_grid_x0.value()
        self.SH_y0 = self.spB_grid_y0.value()
        self.SH_dx = self.spB_grid_dx.value()
        self.SH_dy = self.spB_grid_dy.value()

        isz = 128

        # nx = int((isz - self.SH_x0) / self.SH_dx + 2)
        nx = self.dms + 1
        self.xx = self.SH_x0 + np.arange(nx-1) * self.SH_dx
        print('test lucas np arange', np.arange(nx-1))
        if (isz - self.xx[-1] > self.SH_dx/2+2):
            self.xx = np.append(self.xx, isz)

        # ny = int((isz - self.SH_y0) / self.SH_dy + 2)
        ny = self.dms + 1
        self.yy = self.SH_y0 + np.arange(ny-1) * self.SH_dy
        print('test lucas np arange y', np.arange(ny - 1))
        if (isz - self.yy[-1] > self.SH_dy/2+2):
            self.yy = np.append(self.yy, isz)

        self.new_xx = np.zeros(10, dtype=np.int_)
        self.new_yy = np.zeros(10, dtype=np.int_)
        # DOUBLE THE NUMBER OF ELEMENTS TO DEFINE START AND END CORNER POINTS
        self.pos_xx = np.zeros((2 * nx, 2), dtype=np.int_)
        for i, myx in enumerate(self.xx):
            self.pos_xx[2*i, 0] = myx
            self.pos_xx[2*i, 1] = self.yy[0]
            self.pos_xx[2*i+1, 0] = myx
            self.pos_xx[2*i+1, 1] = min(self.yy[-1], isz) # -1 for display

        self.pos_yy = np.zeros((2 * ny, 2), dtype=np.int_)
        for i, myy in enumerate(self.yy):
            self.pos_yy[2*i,   1] = myy
            self.pos_yy[2*i,   0] = self.xx[0]
            self.pos_yy[2*i+1, 1] = myy
            self.pos_yy[2*i+1, 0] = min(self.xx[-1], isz)
        print('nx,ny : ',nx,ny)


        pos = np.append(self.pos_yy, self.pos_xx, axis=0)
        adj = np.arange(2*(nx+ny)).reshape(nx+ny, 2)

        # pos_test_x = [[]]
        # pos_test_y = [[]]
        # for jj in range(ny - 2):
        #     y0, y1 = self.yy[jj], self.yy[jj + 1]
        #     for ii in range(nx - 2):
        #         x0, x1 = self.xx[ii], self.xx[ii + 1]
        #         # pos_test_x[jj].append(x1 - x0)
        #         # pos_test_y[jj].append(y1 - y0)
        #         # print("c'est le test x y", pos_test_x, pos_test_y)

        self.SH_xx = self.xx
        self.SH_yy = self.yy

        nx = int((isz - self.SH_x0) / self.SH_dx + 2)
        xx = np.round(self.SH_x0+np.arange(nx-1)*self.SH_dx).astype(int)
        if (isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, isz)

        ny = int((isz - self.SH_y0) / self.SH_dy + 2)
        yy = np.round(self.SH_y0 + np.arange(ny - 1) * self.SH_dy).astype(int)
        if (isz - yy[-1] > self.SH_dy / 2 + 2):
            yy = np.append(yy, isz)

        print(self.SH_xx)

        size_array_x = [[0 for a in range(len(yy) - 1)] for a in range(len(yy) - 1)]
        size_array_y = [[0 for a in range(len(yy) - 1)] for a in range(len(yy) - 1)]

        for ii in range(0, len(self.SH_yy) - 1):
            y0, y1 = yy[ii], yy[ii + 1]
            for jj in range(0, len(self.SH_yy) - 1):
                x0, x1 = xx[jj], xx[jj + 1]
                size_xx = x1 - x0
                size_yy = y1 - y0
                print(ii)
                size_array_x[ii][jj] = size_xx
                size_array_y[ii][jj] = size_yy

        print("size_array_x", size_array_x)

        # pos_xx_new = np.zeros((2 * nx, 2), dtype=np.int_)
        # for i, myx in enumerate(size_array_x):
        #     self.pos_xx[2 * i, 0] = myx
        #     self.pos_xx[2 * i, 1] = self.yy[0]
        #     self.pos_xx[2 * i + 1, 0] = myx
        #     self.pos_xx[2 * i + 1, 1] = min(self.yy[-1], isz)  # -1 for display
        #
        # pos_yy_new = np.zeros((2 * ny, 2), dtype=np.int_)
        # for i, myy in enumerate(size_array_y):
        #     self.pos_yy[2 * i, 1] = myy
        #     self.pos_yy[2 * i, 0] = self.xx[0]
        #     self.pos_yy[2 * i + 1, 1] = myy
        #     self.pos_yy[2 * i + 1, 0] = min(self.xx[-1], isz)
        #print(self.yy[-1])

        pos_xx_new = np.zeros((2 * nx, 2), dtype=np.int_)
        print(pos_xx_new)
        act_pos_x = 0
        for i in range(len(size_array_x[0])):
            pos_xx_new[2 * i, 0] = act_pos_x
            pos_xx_new[2 * i, 1] = self.yy[0]
            pos_xx_new[2 * i + 1, 0] = act_pos_x
            pos_xx_new[2 * i + 1, 1] = min(self.yy[-1], isz)
            self.new_xx[i] = act_pos_x
            act_pos_x += size_array_x[0][i]


        pos_yy_new = np.zeros((2 * ny, 2), dtype=np.int_)
        print(pos_yy_new)
        act_pos_y = 128
        ind = 0
        for i in reversed(range(len(size_array_y[0]))):
            self.new_yy[ind] = act_pos_y
            act_pos_y -= size_array_y[i][0]
            ind += 1

        act_pos_y = 0
        for i in range(len(size_array_y[0])):
            pos_yy_new[2 * i, 1] = act_pos_y
            pos_yy_new[2 * i, 0] = self.xx[0]
            pos_yy_new[2 * i + 1, 1] = act_pos_y
            pos_yy_new[2 * i + 1, 0] = min(self.xx[-1], isz)
            act_pos_y += size_array_y[i][0]

        print(size_array_y)
        print(pos_xx_new)
        # print("test", self.pos_xx)

        new_pos = np.append(pos_yy_new, pos_xx_new, axis=0)

        print("new pos", new_pos)

        print('nx,ny : ', nx, ny)


        self.vItemValidCells = np.array([])

        if self.chB_show_grid.isChecked():
            print('pos >>', pos)
            print('adj >>', adj)
            self.overlayGrid.setData(pos=pos, adj=adj,
                                     pen=(0, 255, 0), size=0,
                                     # symbolPen=None, symbolBrush=None,
                                     pxMode=False)
        else:
            self.overlayGrid.setData(pos=new_pos, adj=adj,
                                  pen=(255, 0, 0), size=0,
                                    symbolPen=None, symbolBrush=None,
                                    pxMode=False)

        # if self.chB_show_grid.isChecked():
        #     print('pos >>',pos)
        #     print('adj >>',adj)
        #     self.overlayGrid.setData(pos=new_pos, adj=adj,
        #                          pen=(255, 0, 0), size=0,
        #                          # symbolPen=None, symbolBrush=None,
        #                          pxMode=False)
        # else:
        #     self.overlayGrid.setData(pos=new_pos, adj=adj, pen=None, size=0,
        #                          #symbolPen=None, symbolBrush=None,
        #                          pxMode=False)

            if len(self.vItemValidCells) > 0:
                for i in range(len(self.vItemValidCells)):
                    self.gView_SH_cam.removeItem(self.vItemValidCells[i])

        # if self.bOverlaySelectedCells == True:
        #     for k in range (self.pup_win.nSubAps):
        #         if self.vCurSubApSel[k] == True:


    # =========================================================
    def validate_leakage_gain(self):
        leakage_gain = self.dspB_leakage_gain.value()
        self.mycam.send_fifo_cmd(('system set_leakage_gain ' + str(leakage_gain)).encode())
        self.log("Leakage gain updated = %.4f" % (leakage_gain,), "blue")

    # =========================================================
    def refresh_all(self):
        ''' ----------------------------------------------------------
        Refresh the display

        Independently from the "under the hood" engine, the display
        gets refreshed every now and then, to give visual feedback
        to the user.
        ---------------------------------------------------------- '''
        if self.cam_counter < self.shm_cam.get_counter():
            #print(self.shm_cam.get_counter())
            self.data_cam = np.flipud(self.shm_cam.get_data(check=False, reform=True))
            self.cam_counter = self.shm_cam.get_counter()
            self.imv_raw.setImage(arr2im(self.data_cam.T,
                                         vmin=self.vmin, vmax=self.vmax,
                                         pwr=self.pwr,
                                         cmap=self.mycmap), border=1)

        self.refresh_wfs_display()
        # -------- GET THE MODES HISTORY FROM SHARED MEMORY
        self.data_modes = (np.asarray(self.shm_modes.get_data(check=False, reform=True, sleepT=0.)))[0:50,0:self.ZER_WFC.nmodes]
        # print('data_modes.shape >> ', self.data_modes.shape)
        self.proj_ttx = self.data_modes[:,0] * self.ZER_calAmp
        self.proj_tty = self.data_modes[:,1] * self.ZER_calAmp
        self.logplotx.setData(self.proj_ttx)
        self.logploty.setData(self.proj_tty)

        # GET THE MOST RECENT MEASURED MODE PROJECTIONS VALUES
        if False:
            mode_amp = self.data_modes[49,:] * self.ZER_calAmp
            self.mode_graphA.setOpts(x=self.bg_vx, height=mode_amp)
            if 'mode_amp_res' in locals():
                self.mode_graphB.setOpts(x=self.bg_vx+0.3, height=self.mode_amp_res)
            else:
                self.mode_graphB.setOpts(x=self.bg_vx+0.3, height=self.yB)

            # for k in range(2, self.ZER_WFC.modes.shape[0]):
            #      = np.sum(self.shm_HOmodes * DMact_modes[k*nbAct:(k+1)*nbAct] = tmpModes[k,self.ZER_WFC.DM.vActMapping]



        # GET THE LOW ORDER (TT) COMMAND MAP AND PROJECT IT ON THE TT CALIBRATION MODES TO OBTAIN THE DM TT ABSOLUTE VALUES
        if self.bTCSOffloadInvertX == True:
            signX = -1
        else:
            signX = +1

        if self.bTCSOffloadInvertY == True:
            signY = -1
        else:
            signY = +1


        self.shm_LOmodes.read_keyword(0)
        self.shm_LOmodes.read_keyword(1)
        LO_cmd = self.shm_LOmodes.get_data(check=False).flatten()
        # kw_TTx = self.shm_LOmodes.kwds[0]['value']
        kw_TTx = np.sum(LO_cmd * (self.ZER_WFC.modes[0,:,:].flatten())[self.ZER_WFC.DM.vActMapping]) * self.ZER_norm_calib_modes[0]
        # kw_TTy = self.shm_LOmodes.kwds[1]['value']
        kw_TTy = np.sum(LO_cmd * (self.ZER_WFC.modes[1,:,:].flatten())[self.ZER_WFC.DM.vActMapping]) * self.ZER_norm_calib_modes[1]
        self.ttx_drift_avg = signX * self.OffloadPixScale * kw_TTx * self.ZER_calAmp
        self.tty_drift_avg = signY * self.OffloadPixScale * kw_TTy * self.ZER_calAmp
        # print('kw_TTx, kw_TTy >> ', kw_TTx, kw_TTy)

        self.yB[0] = kw_TTx
        self.yB[1] = kw_TTy

        HO_cmd = self.shm_HOmodes.get_data(check=False).flatten()
        for k in range(2, self.ZER_WFC.nmodes):
            self.shm_HOmodes.read_keyword(k-2)
            # kw_HO = self.shm_HOmodes.kwds[k]['value']
            kw_HO = np.sum(HO_cmd * (self.ZER_WFC.modes[k,:,:].flatten())[self.ZER_WFC.DM.vActMapping]) * self.ZER_norm_calib_modes[k]
            # print('k, kw_HO >> ', k, kw_HO)
            self.vHO_global_amp[k-2] = kw_HO * self.ZER_calAmp
            self.yB[k] = kw_HO

        timestamps = (np.asarray(self.shm_timestamps.get_data(check=False, reform=True, sleepT=0.)))[0:50,:]
        #print((timestamps[0][100] - timestamps[0][99]) / 10000000)
        self.refresh_stats()


    # =========================================================
    def refresh_wfs_display(self):
        ''' ----------------------------------------------------------
        Specific display refresh for the WFS data (ttxy + phot)
        ---------------------------------------------------------- '''
        # self.imv_WFS.setImage(
        #     arr2im(np.concatenate(self.wfs.SH_comb),
        #            vmin=-5, vmax=5, pwr=1, cmap=self.mycmap),
        #     border=1)

        self.data_slopes = (np.asarray(self.shm_slopes.get_data(check=False, reform=True, sleepT=0.)))[0:50,:]
        # print('data_slopes >> ',self.data_slopes.shape)
        self.imv_Slopes.setImage(
            arr2im(self.data_slopes[0:50,200:400].T,
                    vmin=-5., vmax=+5.,
                    pwr=1, cmap=self.mycmap),
                    border=1)

        last_slope_map = np.zeros([10,30])
        self.last_slope_x = self.data_slopes[49,200:300]
        self.last_slope_y = self.data_slopes[49,300:400]
        self.last_illum = self.data_slopes[49,400:500]
        last_slope_map[0:10,0:10] = np.flipud(self.last_slope_x.reshape([10,10]))
        last_slope_map[0:10,10:20] = np.flipud(self.last_slope_y.reshape([10,10]))
        last_slope_map[0:10,20:30] = np.flipud(self.last_illum.reshape([10,10]))
        # NORMALIZE ILLUMINATION TO MAX
        if np.sum(last_slope_map[0:10,20:30]) > 0:
            last_slope_map[0:10,20:30] = last_slope_map[0:10,20:30] / np.max(last_slope_map[0:10,20:30])
        else:
            last_slope_map[0:10,20:30] = 0
        self.imv_WFS.setImage(arr2im(last_slope_map.T,vmin=-5, vmax=5, pwr=1, cmap=self.mycmap),border=1)


        # UPDATE CENTROID AND REFERENCE POSITIONS
        reference_pos = np.zeros(((self.dms-1)*(self.dms-1), 2), dtype=np.float32)
        centroid_pos = np.zeros(((self.dms-1)*(self.dms-1), 2), dtype=np.float32)
        slopes = np.zeros(((self.dms-1)*(self.dms-1), 2), dtype=np.float32)

        try:
            _ = self.data_slopes.shape
        except:
            print('"self.data_slopes" has not been initialized.')
            return

        # centroid_pos IS ONLY VALID IN DISPLAY SPACE (VIZUALIZATION ONLY!!)
        centroid_pos[:,0] = self.data_slopes[49,0:100]
        centroid_pos[:,1] = self.data_slopes[49,100:200]
        slopes[:,0] = self.data_slopes[49,200:300]
        slopes[:,1] = np.flipud(self.data_slopes[49,300:400])
        # print('reference_pos/ref_centroids shapes >> ',reference_pos.shape,self.ref_centroids.shape)
        reference_pos[:,0] = self.ref_centroids[0:100]
        reference_pos[:,1] = self.ref_centroids[100:200]
        # print(reference_pos)
        flipped_yy = np.flipud(self.new_yy)
        flipped_xx = np.flipud(self.xx)
        x = 0
        y = 0
        # print("slopes: ", slopes)
        for j in reversed(range(self.dms-1)):
            for i in reversed(range(self.dms-1)):
                crd = i + j*(self.dms-1)
                # centroid_pos[crd,0] = centroid_pos[crd,0] + self.xx[i]
                # centroid_pos[crd,1] = centroid_pos[crd,1] + self.yy[j]
                # reference_pos[crd,0] = reference_pos[crd,0] + self.xx[i]
                # reference_pos[crd,1] = reference_pos[crd,1] + self.yy[j]
                #print(crd, self.new_yy[y])
                #print(crd,centroid_pos[crd,0] + self.new_xx[i], centroid_pos[crd,1] + self.new_yy[j], self.new_xx[i], self.new_yy[j])
                centroid_pos[crd,0] = centroid_pos[crd,0] + self.new_xx[i] + 0.5
                centroid_pos[crd,1] = self.new_yy[j] - centroid_pos[crd, 1] - 0.5
                reference_pos[crd,0] = reference_pos[crd,0] + self.new_xx[i] + 0.5
                reference_pos[crd,1] = self.new_yy[j] - reference_pos[crd,1] - 0.5
                # print(crd, centroid_pos[crd,0], centroid_pos[crd,1])
                x += 1
            y += 1
            x = 0
        # centroid_y = self.data_slopes[49,100:200]
        # centroid_pos = np.append(centroid_y, centroid_x, axis=0)
        # print('reference_pos >> ',centroid_pos[15,0],centroid_pos[15,1])
        dpi = QGuiApplication.primaryScreen().logicalDotsPerInch()
        screenWidth = QGuiApplication.primaryScreen().size().width()
        penWidth = max(1, int(1 / (dpi / screenWidth)))
        if self.check_ShowRefPos.isChecked():
            symBrush = pg.mkBrush([255,255,255,127])
            crossSize = 1
            customSymbol = QtGui.QPainterPath()
            customSymbol.moveTo(-crossSize, 0)
            customSymbol.lineTo(crossSize, 0)
            customSymbol.moveTo(0, -crossSize)
            customSymbol.lineTo(0, crossSize)
            pen = QtGui.QPen(QtGui.QColor(255, 0, 0))
            pen.setWidth(penWidth)
            self.overlayRefs.setData(pos=reference_pos,
                                pen=pen, symbol=customSymbol, symbolBrush=None, size=1,
                                pxMode=False, antialias=True, compositionMode=QtGui.QPainter.CompositionMode.CompositionMode_SourceOver)
        else:
            self.overlayRefs.setData(pos=reference_pos,
                                pen=None, size=0,
                                pxMode=False, compositionMode=QtGui.QPainter.CompositionMode.CompositionMode_SourceOver)

        if self.check_ShowCentroids.isChecked():
            crossSize = 1
            customSymbol = QtGui.QPainterPath()
            customSymbol.moveTo(-crossSize, -crossSize)
            customSymbol.lineTo(crossSize, crossSize)
            customSymbol.moveTo(-crossSize, crossSize)
            customSymbol.lineTo(crossSize, -crossSize)
            pen = QtGui.QPen(QtGui.QColor(255, 0, 0))
            pen.setWidth(penWidth)
            self.overlayCentroids.setData(pos=centroid_pos,
                                pen=pen, symbol=customSymbol, symbolBrush=None, size=1,
                                pxMode=False, ntialias=True, compositionMode=QtGui.QPainter.CompositionMode.CompositionMode_SourceOver)
        else:
            self.overlayCentroids.setData(pos=centroid_pos,
                                pen=None, size=0,
                                pxMode=False, compositionMode=QtGui.QPainter.CompositionMode.CompositionMode_SourceOver)

        # self.bShowCellSelectionFeedback = True
        if self.bShowCellSelectionFeedback == True:
            myListPos = []
            myListAdj = []
            ptCounter = 0
            selCount = 0
            isz = 128
            for k in range(0, self.nSubAps):
                if self.vCurSubApSel[k] == True:
                    ix = k % 10
                    iy = k // 10
                    myListPos.append([[self.xx[ix]+1,isz-(self.yy[iy]+1)],[self.xx[ix]+1,isz-(self.yy[iy+1]-2)],[self.xx[ix+1]-2,isz-(self.yy[iy+1]-2)],[self.xx[ix+1]-2,isz-(self.yy[iy]+1)]])
                    myListAdj.append([[ptCounter,ptCounter+1],[ptCounter+1,ptCounter+2],[ptCounter+2,ptCounter+3],[ptCounter+3,ptCounter]])
                    ptCounter += 4
                    selCount += 1

            if selCount > 0:
                vCellPos = np.array(myListPos[0])
                vCellAdj = np.array(myListAdj[0]).astype(np.int32)
                for k in range(1, selCount):
                    vCellPos = np.concatenate((vCellPos, np.array(myListPos[k])), axis=0)
                    vCellAdj = np.concatenate((vCellAdj, np.array(myListAdj[k]).astype(np.int32)), axis=0)
                # print('vCellPos')
                # print(vCellPos)
                # print('vCellAdj')
                # print(vCellAdj)

                # vCellPos = np.array([[self.xx[0]+1,self.yy[0]+1],[self.xx[0]+1,self.yy[1]-2],[self.xx[1]-2,self.yy[1]-2],[self.xx[1]-2,self.yy[0]+1]])
                # vAdj = np.array([[0,1],[1,2],[2,3],[3,0]])

                self.overlayPupil.setData(pos=vCellPos, adj=vCellAdj,
                                    pen=(0,255,0,255), size=0,
                                    # symbolPen=None, symbolBrush=None,
                                    pxMode=False)
            else:
                self.overlayPupil.setData(pen=None, pxMode=False)
        else:
            self.overlayPupil.setData(pen=None, pxMode=False)

    # =========================================================
    def update_cbar(self):
        ''' ----------------------------------------------------------
        Updates the colorbar used to display images in the GUI
        ---------------------------------------------------------- '''
        cbar = str(self.cmB_cbar.currentText()).lower()
        try:
            exec('self.mycmap = cm.%s' % (cbar,))
        except:
            print("colormap %s is not available" % (cbar,))
            self.mycmap = cm.jet

    # =========================================================
    def update_vmin(self):
        ''' ----------------------------------------------------------
        Sets a new minimum value for the SH image display
        ---------------------------------------------------------- '''

        # if self.bInitDone == True:
        #     # self.vSubApThresh = self.pup_win.vCellThreshold
        #     self.vSubApThresh = np.full((self.dms-1)*(self.dms-1), self.dspB_disp_min.value()).astype(np.float32)
        #     print('self.wfs.nCells >> ', (self.dms-1)*(self.dms-1))
        #     self.vSubApThresh[0] = 10000;
        #     self.vSubApThresh[1] = 10000;
        #     self.vSubApThresh[2] = 10000;
        #     self.vSubApThresh[7] = 10000;
        #     self.vSubApThresh[8] = 10000;
        #     self.vSubApThresh[9] = 10000;
        #     self.vSubApThresh[10] = 10000;
        #     self.vSubApThresh[11] = 10000;
        #     self.vSubApThresh[18] = 10000;
        #     self.vSubApThresh[19] = 10000;
        #     self.vSubApThresh[20] = 10000;
        #     self.vSubApThresh[29] = 10000;
        #     self.vSubApThresh[70] = 10000;
        #     self.vSubApThresh[79] = 10000;
        #     self.vSubApThresh[80] = 10000;
        #     self.vSubApThresh[81] = 10000;
        #     self.vSubApThresh[88] = 10000;
        #     self.vSubApThresh[89] = 10000;
        #     self.vSubApThresh[90] = 10000;
        #     self.vSubApThresh[91] = 10000;
        #     self.vSubApThresh[92] = 10000;
        #     self.vSubApThresh[97] = 10000;
        #     self.vSubApThresh[98] = 10000;
        #     self.vSubApThresh[99] = 10000;
        #     fits.writeto('R:\\ixon_subaperture_thresholds.fits',self.vSubApThresh,overwrite=True)
        #     self.mycam.send_fifo_cmd(b'ixon set_thresh \"R:\\ixon_subaperture_thresholds.fits\"')

        #self.cam_counter -= 1
        self.vmin = False
        if self.chB_min.isChecked():
            self.vmin = self.dspB_disp_min.value()

            if self.chB_AutoSetWFSThresh.isChecked():
                self.spB_grid_imin.setValue(self.vmin)

    # =========================================================
    def update_vmax(self):
        ''' ----------------------------------------------------------
        Sets a new maximum value for the SH image display
        ---------------------------------------------------------- '''
        self.cam_counter -= 1
        self.vmax = False
        if self.chB_max.isChecked():
            self.vmax = self.dspB_disp_max.value()

    # =========================================================
    def update_pupil_thresholds(self):
        ''' ----------------------------------------------------------
        Updates the per-subapertures threshold values by updating the
        "R:\\ixon_subaperture_thresholds.fits" file.
        ---------------------------------------------------------- '''
        fits.writeto('R:\\ixon_subaperture_thresholds.fits',self.vSubApThresh,overwrite=True)
        fits.writeto('R:\\ixon_subaperture_enabled.fits',self.vSubApEnable,overwrite=True)
        self.mycam.send_fifo_cmd(b'ixon set_thresh \"R:\\ixon_subaperture_thresholds.fits\" \"R:\\ixon_subaperture_enabled.fits\"')

        # fits.writeto('R:\\ixon_subaperture_enabled.fits',self.vSubApEnable,overwrite=True)
        # self.mycam.send_fifo_cmd(b'ixon set_thresh \"R:\\ixon_subaperture_enabled.fits\"')

        print('Threshold map updated.')


    # =========================================================
    def update_AutoSetWFSThresh(self):
        if self.chB_AutoSetWFSThresh.isChecked():
            self.spB_grid_imin.setValue(int(self.vmin))

    # =========================================================
    def update_nonlinear(self):
        ''' ----------------------------------------------------------
        Switch between linear and square root display for SH image
        ---------------------------------------------------------- '''
        self.cam_counter -= 1
        self.pwr = 1.0
        if self.chB_nonlinear.isChecked():
            self.pwr = 0.5

    # =========================================================
    def log(self, message="", color="black", crt=True):
        ''' --------------------------------------------------
        Convenient logging function.

        Parameter:
        - msg   : the string to log (html text)
        - color : a valid html color name (default "black"
        - crt   : a flag for carriage return (default: True)

        Add a time stamp !
        Write to a log file !
        -------------------------------------------------- '''
        global myqt
        myline = "<b><font color='%s'>%s %s</font></b>" % (
            color, time.strftime('%D %H:%M:%S'), message)
        if crt:
            myline += "<br>"
        try:
            myqt.gui_do(self.text_log.insertHtml, myline)
            with open(logfile, "a") as mylogfile:
                mylogfile.write(myline+"\n")

        except:
            print("'%s' did not log?" % (myline,))
            pass

        # these two lines scroll the text display down
        temp = self.text_log.textCursor()
        myqt.gui_do(self.text_log.setTextCursor, temp)

    # =========================================================
    #        wavefront control callbacks !!
    # =========================================================

    # =========================================================
    # ---                      TIP-TILT                     ---
    # =========================================================

    def tt_abort(self):
#        msg = self.TT_WFC.abort()
        if self.ZER_WFC.keepgoing == False:
            self.bDoMeasure = False
        self.TT_WFC.keepgoing = False
#        self.log(msg, "red")
        self.enable_wfc(subsystem="tt")

    # =========================================================
    def tt_reset_correc(self):
        self.TT_WFC.reset()
        self.log("TT correction reset!", "green")

    # =========================================================
    def tt_calib(self):
        self.TT_WFC.keepgoing = False
        nav = self.spB_nav.value()
        a0 = self.dspB_tt_a0.value()
        self.TT_WFC.nav = nav
        self.disable_wfc(subsystem="tt")
        self.log("TT loop calibration starts", "blue")
        self.log("cal: nav=%d, a0=%.2f" % (nav, a0), "blue")
        self.ttcal_Thread = GenericThread(self.TT_WFC.calibrate,
                                          a0=a0, reform=False)
        self.wfc_calibrating = "TT"
        self.ttcal_Thread.start()

    # =========================================================
    def tt_update_gain(self):
        gain = self.dspB_tt_gain.value()
        self.TT_WFC.gain = gain
        self.log("TT loop gain updated =%.2f" % (gain,), "blue")

    # =========================================================
    def tt_cloop(self):
        self.bCloseLoop = True
        self.TT_WFC.bApplyCorrection = True
        self.bDoMeasure = True
        self.corMode = 0
        # self.nmodes = self.TT_WFC.nmodes
        self.ttloop_Thread = GenericThread(self.TT_WFC.cloop)
        self.log("TT closed loop!", "green")
        self.ttloop_Thread.start()
        self.disable_wfc(subsystem="tt")

    # =========================================================
    # ---                     Zernike                       ---
    # =========================================================
    def zer_abort(self):
#        msg = self.ZER_WFC.abort()
#        self.log(msg, "red")
        if self.TT_WFC.keepgoing == False:
            self.bDoMeasure = False
        self.ZER_WFC.keepgoing = False
        self.enable_wfc(subsystem="zer")
        self.log("Zernike loop stopped!", "green")

    # =========================================================
    def zer_reset_correc(self):
        self.ZER_WFC.reset()
        self.log("Zernike correction reset!", "green")

    # =========================================================
    def zer_calib(self):
        self.calib_mode = 'ZER'
        nmodes = self.ZER_WFC.modes.shape[0]
        nCalFrames = self.spB_nav.value()
        szMode = self.ZER_WFC.modes.shape[1]*self.ZER_WFC.modes.shape[2]
        nbAct = self.ZER_WFC.DM.nba
        self.ZER_calAmp = self.dspB_zer_a0.value()
        self.ZER_WFC.gain = self.ZER_calAmp

        tmpModes = self.ZER_calAmp * self.ZER_WFC.modes.reshape((nmodes, szMode)).astype(np.float32)
        DMact_modes = np.zeros(nmodes*nbAct).astype(np.float32)
        print('self.ZER_WFC.modes.shape >> ', tmpModes.shape)

        for k in range(nmodes):
            DMact_modes[k*nbAct:(k+1)*nbAct] = tmpModes[k,self.ZER_WFC.DM.vActMapping]
            self.ZER_norm_calib_modes[k] = 1./np.sum(tmpModes[k,self.ZER_WFC.DM.vActMapping]*tmpModes[k,self.ZER_WFC.DM.vActMapping])

        calFile = 'R:\\ZER_CAL.fits'
        fits.writeto(calFile, DMact_modes.reshape((nmodes, nbAct)), overwrite=True)
        fits.setval(calFile, 'NMODES', value=nmodes)
        fits.setval(calFile, 'NCALFRM', value=nCalFrames)
        self.lib_EventManager.trigger_named_event(('_AOC_ZER_CAL_READY').encode('utf-8'))

        if self.lib_EventManager.wait_for_single_event(('_AOC_OK').encode('utf-8'), 5000) != 0:
            self.log("[zer_calib] TIMEOUT while waiting for '_AOC_OK' event. Calibration aborted.", "red")
        else:
            self.mycam.send_fifo_cmd(b"system calib_zer")
            os.write(self.mycam.cmd_fifo, b'system query_log_dir')
            os.fsync(self.mycam.cmd_fifo)
            tmp = os.read(self.mycam.cmd_fifo, self.mycam.MAX_FIFO_LEN)
            os.fsync(self.mycam.cmd_fifo)
            print(tmp)
            calFileCpy = str(tmp) + 'calib\\' + 'ZER_CAL.fits'
            calFileCpy = calFileCpy.replace("b'", "")
            calFileCpy = calFileCpy.replace("'", "")
            self.log_dir = str(tmp) + 'calib\\'
            self.log_dir = self.log_dir.replace("b'", "")
            self.log_dir = self.log_dir.replace("'", "")
            fits.writeto(calFileCpy, DMact_modes.reshape((nmodes, nbAct)), overwrite=True)
            fits.setval(calFileCpy, 'NMODES', value=nmodes)
            fits.setval(calFileCpy, 'NCALFRM', value=nCalFrames)
            fits.writeto(self.log_dir + 'ixon_subaperture_thresholds.fits', self.vSubApThresh, overwrite=True)
            fits.writeto(self.log_dir + 'ixon_subaperture_enabled.fits', self.vSubApEnable, overwrite=True)
            fits.writeto(self.log_dir + 'ActuatorsMappingIndices.fits', self.vSubApEnable, overwrite=True)
            shutil.copy('R:\\default_ref_pos_x.fits', self.log_dir)
            shutil.copy('R:\\default_ref_pos_y.fits', self.log_dir)


    def process_calibration(self, calib_mode):
        if calib_mode == 'ZER':
            print('Processing ZER calibration')
            self.RR = fits.getdata(r'R:\ZER_RecModes.fits')
            self.ZER_calib_done = True
        elif calib_mode == 'ZON':
            print('Processing ZON calibration')
            self.RR = fits.getdata(r'R:\ZON_RecModes.fits')
            self.ZON_calib_done = True

        self.RTR = self.RR.dot(self.RR.T)
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        u, s, vh = np.linalg.svd(self.RR.T, full_matrices=False, compute_uv=True, hermitian=False)
        print('u.shape, s.shape, vh.shape >> ', u.shape, s.shape, vh.shape)
        self.RRsvd = vh.T.dot(np.diag(1./s).dot(u.T))
        print('Singular values >> ', s)
        fits.writeto(r'R:\RRinv.fits',self.RRinv,overwrite=True)
        fits.writeto(r'R:\RRsvd.fits',self.RRsvd,overwrite=True)

        # SANITY CHECK
        fits.writeto(r'R:\RRinvRRT.fits',np.dot(self.RRinv,self.RR.T),overwrite=True)
        fits.writeto(r'R:\RRsvdRRT.fits',np.dot(self.RRsvd,self.RR.T),overwrite=True)

        # COPY IN LOG FILES
        fits.writeto(self.log_dir + 'RRinv.fits', self.RRinv, overwrite=True)
        fits.writeto(self.log_dir + 'RRsvd.fits', self.RRsvd, overwrite=True)

        # SANITY CHECK
        fits.writeto(self.log_dir + 'RRinvRRT.fits', np.dot(self.RRinv, self.RR.T), overwrite=True)
        fits.writeto(self.log_dir + 'RRsvdRRT.fits', np.dot(self.RRsvd, self.RR.T), overwrite=True)

        self.bDoMeasure = True

    # =========================================================
    def zer_update_gain(self):
        gain = self.dspB_zer_gain.value()
        self.ZER_WFC.gain = gain
        self.mycam.send_fifo_cmd(('system set_ZER_gain ' + str(gain)).encode())
        self.log("ZERNIKE loop gain updated =%.2f" % (gain,), "blue")

    # =========================================================
    def zer_measure(self):
        self.bDoMeasure = True
        self.corMode = 1
        self.nmodes = self.ZER_WFC.nmodes
        self.zerloop_Thread = GenericThread(self.ZER_WFC.cloop)
        self.log("Zernike closed loop!", "green")
        self.zerloop_Thread.start()
        self.disable_wfc(subsystem="zer")

    # =========================================================
    def zer_cloop(self):
        if self.chB_zer_cloop.isChecked():
            # self.LogDataSlopes = np.zeros([self.log_length,3,10,10],dtype='float32')
            # # self.LogDataErrors = np.zeros([self.log_length,self.nmodes])
            # # self.LogDataCmds = np.zeros([self.log_length,self.nmodes])
            # self.LogDataDMProjModes = np.zeros([self.log_length,self.nmodes])
            # self.LogTimeStamp = np.zeros(self.log_length,dtype='float64')

            self.bCloseLoop = True
            # self.ZER_WFC.bApplyCorrection = True
            self.mycam.send_fifo_cmd(b'alpao direct_dm_update ON')
            print('Direct DM update set to ON')
            self.mycam.send_fifo_cmd(b'system close_loop_ZER ON')
            print('Close loop enabled')
        else:
            self.bCloseLoop = False
            # self.ZER_WFC.bApplyCorrection = False
            self.mycam.send_fifo_cmd(b'alpao direct_dm_update OFF')
            print('Direct DM update set to OFF')
            self.mycam.send_fifo_cmd(b'system close_loop_ZER OFF')
            print('Close loop disabled')
            self.zer_reset_correc()

        self.bDoMeasure = True

        # self.corMode = 1
        # self.nmodes = self.ZER_WFC.nmodes
        # self.zerloop_Thread = GenericThread(self.ZER_WFC.cloop)
        # self.log("Zernike closed loop!", "green")
        # self.zerloop_Thread.start()
        # self.disable_wfc(subsystem="zer")

    # =========================================================
    def zer_update_nz(self):
        print('Number of ZER_WFC instances: %d' % sys.getrefcount(self.ZER_WFC))
        self.nz = self.dspB_zer_nz.value()
        # self.ZER_WFC = ZER_WFC(iz0=4, iz1=self.nz)
        self.ZER_WFC = ZER_WFC(iz0=2, iz1=self.nz, simu=self.simu)
        self.ZER_WFC.nmodes = self.nz - 1
        self.vHO_global_amp = np.zeros(self.ZER_WFC.nmodes - 2)
        print('Number of ZER_WFC modes: %d' % self.ZER_WFC.nmodes)
        self.ZER_WFC.zer_mode_amp = np.zeros(self.ZER_WFC.nmodes)
        self.ZER_WFC.zer_mode_amp_resid = np.zeros(self.ZER_WFC.nmodes)
        self.bg_vx = 0.5 + np.arange(self.ZER_WFC.nmodes)
        yA = np.linspace(0.0, 3.0, num=self.ZER_WFC.nmodes)
        yB = np.linspace(0.0, 3.0, num=self.ZER_WFC.nmodes)
        self.mode_graphA.setOpts(x=self.bg_vx,height=yA)
        self.mode_graphB.setOpts(x=self.bg_vx,height=yB)
        self.log("Zernike basis of modes updated. Recalibrate", "orange")

    # =========================================================
    # ---                     ZONAL                         ---
    # =========================================================
    def zo_abort(self):
#        msg = self.ZON_WFC.abort()
#        self.log(msg, "red")
        self.bDoMeasure = False
        self.ZON_WFC.keepgoing = False
        self.enable_wfc(subsystem="zo")

    # =========================================================
    def zo_reset_correc(self):
        self.ZON_WFC.reset()
        self.log("ZON correction reset!", "green")

    # =========================================================
    def zo_calib(self):
        # self.ZON_WFC.keepgoing = False
        # nav = self.spB_nav.value()
        # a0 = self.dspB_zo_a0.value()
        # self.ZON_WFC.nav = nav
        # self.disable_wfc(subsystem="zo")
        # self.log("ZON loop calibration starts", "blue")
        # self.log("cal: nav=%d, a0=%.2f" % (nav, a0), "blue")
        # self.zocal_Thread = GenericThread(self.ZON_WFC.calibrate,
        #                                   a0=a0, reform=False)
        # self.wfc_calibrating = "ZON"
        # self.zocal_Thread.start()
        #
        # self.bg_vx = 0.5 + np.arange(self.ZON_WFC.nmodes)
        # yA = np.linspace(0.0, 3.0, num=self.ZON_WFC.nmodes)
        # yB = np.linspace(0.0, 3.0, num=self.ZON_WFC.nmodes)
        # self.mode_graphA.setOpts(x=self.bg_vx,height=yA)
        # self.mode_graphB.setOpts(x=self.bg_vx,height=yB)
        nmodes = self.ZON_WFC.modes.shape[0]
        nCalFrames = self.spB_nav.value()
        szMode = self.ZON_WFC.modes.shape[1]*self.ZON_WFC.modes.shape[2]
        nbAct = self.ZON_WFC.DM.nba
        self.ZON_calAmp = self.dspB_zo_a0.value()

        tmpModes = self.ZON_calAmp * self.ZON_WFC.modes.reshape((nmodes, szMode))
        DMact_modes = np.zeros(nmodes*nbAct)
        print('self.ZON_WFC.modes.shape >> ', tmpModes.shape)

        for k in range(nmodes):
            # print('tmpModes >> ', k, tmpModes[k,:])
            DMact_modes[k*nbAct:(k+1)*nbAct] = tmpModes[k,self.ZON_WFC.DM.vActMapping]

        calFile = 'R:\\ZON_CAL.fits'
        fits.writeto(calFile, DMact_modes.reshape((nmodes, nbAct)), overwrite=True)
        fits.setval(calFile, 'NMODES', value=nmodes)
        fits.setval(calFile, 'NCALFRM', value=nCalFrames)
        self.lib_EventManager.trigger_named_event(('_AOC_ZON_CAL_READY').encode('utf-8'))

    # =========================================================
    def zo_update_gain(self):
        gain = self.dspB_zo_gain.value()
        self.ZON_WFC.gain = gain
        self.log("Zonal loop gain updated =%.2f" % (gain,), "blue")

    # =========================================================
    def zo_measure(self):
        self.bDoMeasure = True
        self.corMode = 2
        # self.ZON_WFC.bApplyCorrection = False
        self.nmodes = self.ZON_WFC.nmodes
        self.zonloop_Thread = GenericThread(self.ZON_WFC.cloop)
        self.log("ZON closed loop!", "green")
        self.zonloop_Thread.start()
        self.disable_wfc(subsystem="zo")

    # =========================================================
    def zo_cloop(self):
        # self.bDoMeasure = True
        # self.corMode = 2
        # self.ZON_WFC.bApplyCorrection = True
        # self.nmodes = self.ZON_WFC.nmodes
        # self.zonloop_Thread = GenericThread(self.ZON_WFC.cloop)
        # self.log("ZON closed loop!", "green")
        # self.zonloop_Thread.start()
        # self.disable_wfc(subsystem="zo")
        if self.chB_zo_cloop.isChecked():
            self.bCloseLoop = True
            self.ZON_WFC.bApplyCorrection = True
            print('Close loop enabled')
        else:
            self.bCloseLoop = False
            self.ZON_WFC.bApplyCorrection = False
            print('Close loop disabled')

        self.bDoMeasure = True

    # =========================================================
    #         call-backs for tip-tilt offsets
    # =========================================================

    def ttx_offset_p(self):
        self.tt_offset(axis="x", sign="-")
    def ttx_offset_m(self):
        self.tt_offset(axis="x", sign="+")
    def tty_offset_p(self):
        self.tt_offset(axis="y", sign="-")
    def tty_offset_m(self):
        self.tt_offset(axis="y", sign="+")
    # =========================================================

    def tt_offset(self, axis="x", sign="+"):
        ''' ----------------------------------------
        processes the offset requests from user
        ---------------------------------------- '''

        exec("self.wfs.SH_%sref %s= self.ostep" % (axis, sign))
        exec("self.%soff %s= self.ostep" % (axis, sign))

        if axis == "x": # this is not pretty but works for now
            exec("self.log('%s-axis pointing offset = %.1f pixels', 'orange')" % (axis.upper(), self.xoff))
        else:
            exec("self.log('%s-axis pointing offset = %.1f pixels', 'orange')" % (axis.upper(), self.yoff))

    # =========================================================
    def ixon_shutdown(self):
        print('Shutting down iXon!')
        self.mycam.quit()
        self.enable_cam_gui(state=False)
        self.log("iXon camera was shut down", "blue")
        self.bCamStreamOn = False

    # =========================================================
    def enable_cam_gui(self, state=True):
        ''' Enables/disables *all* buttons of the camera GUI'''

        self.pB_cam_stream.setEnabled(state)
        self.pB_cam_stop.setEnabled(state)
        self.pB_cam_tint_inc.setEnabled(state)
        self.pB_cam_tint_dec.setEnabled(state)
        self.pB_cam_shutdown.setEnabled(state)


    # =========================================================
    def alpao_shutdown(self):
        #os.system("touch "+ciao_home+"dm_stop")
        os.system("copy NUL " + home + "\\bin\\dm_stop")
        self.log("alpao DM was shut down", "blue")

    # =========================================================
    def alpao_start(self):
        os.system("cd " + ciao_home + "\\dm_server" + "& start /B .\\alpao_server")
        self.log("alpao DM was turned on!", "blue")


    def update_OffloadPixScale(self,):
        self.OffloadPixScale = self.dspB_OffloadPixScale.value()

    # =========================================================
    def startThread_loop_update_DMstatic(self,):
        ''' ----------------------------------------------------------
        Starts the process of monitoring the telescope pointing data
        ---------------------------------------------------------- '''
        # self.telThread = GenericThread(self.tel.loop)
        self.Thread_DMStaticUpdate = GenericThread(self.Thread_DMStaticUpdate_loop)
        self.Thread_DMStaticUpdate.start()
        self.log("Static DM map update thread created and started", "green")

    def Thread_DMStaticUpdate_loop(self,):
        self.DMStaticUpdate_loop_keepgoing = True
        e = threading.Event()
        while self.DMStaticUpdate_loop_keepgoing:
            event_is_set = e.wait(self.DMstaticUpdate_period_sec)
            if event_is_set:
                print('processing event')
            else:
                self.mycam.send_fifo_cmd(b'alpao update_static_maps')
                print('>> DM static map updated <<')

    # =========================================================
    def tel_mon_start(self,):
        ''' ----------------------------------------------------------
        Starts the process of monitoring the telescope pointing data
        ---------------------------------------------------------- '''
        # self.telThread = GenericThread(self.tel.loop)
        self.telThread = GenericThread(self.tel_loop)
        self.telThread.start()
        self.log("Telescope monitoring resumed", "green")

    def tel_loop(self,):
        self.tel_loop_keepgoing = True
        e = threading.Event()
        # sign = 1.
        while self.tel_loop_keepgoing:
            event_is_set = e.wait(self.tel_mon_period_sec)

            # print('event set: %s' % event_is_set)
            if event_is_set:
                print('processing event')
            else:
                print('TTX/Y last avg: %.3e %.3e' % (self.ttx_drift_avg,self.tty_drift_avg))
                if (np.abs(self.ttx_drift_avg) <= 10.) and (np.abs(self.tty_drift_avg) <= 10):
                    self.tel.offload_tiptilt(self.ttx_drift_avg, self.tty_drift_avg, self.dspB_DerotOffset.value() + self.derotator_angle_deg)
                else:
                    print('Offload corrections too large (>10arcsec). Nothing done.')
                    self.ttx_drift_avg = 0.
                    self.tty_drift_avg = 0.
                    self.chB_TCSOffload.setChecked(False)



    def TCSOffload(self,):
        self.ttx_drift_avg = 0.
        self.tty_drift_avg = 0.

        if self.chB_TCSOffload.isChecked():
            self.tel.bEnableExtTTOffload = True
            print('TCS offload enabled and all values zeroed.')
        else:
            self.tel.bEnableExtTTOffload = False
            print('TCS offload disabled and all values zeroed.')

    def OffloadInvertX(self):
        if self.check_InvertX.isChecked() == True:
            self.bTCSOffloadInvertX = True
        else:
            self.bTCSOffloadInvertX = False

    def OffloadInvertY(self):
        if self.check_InvertY.isChecked() == True:
            self.bTCSOffloadInvertY = True
        else:
            self.bTCSOffloadInvertY = False

    def reset_TCS_offload(self):
        self.ttx_drift_avg = 0.
        self.tty_drift_avg = 0.

    # =========================================================
    def thread_eventmon_start(self,):
        ''' ----------------------------------------------------------
        Starts the process of monitoring the telescope pointing data
        ---------------------------------------------------------- '''
        # self.telThread = GenericThread(self.tel.loop)
        self.event_mon_loop_end = False
        self.thread_eventmon = GenericThread(self.thread_eventmon_loop)
        self.thread_eventmon.start()
        self.log("Event monitoring thread started.", "green")

    def make_clist(self, lst):
        return (ctypes.c_char_p * len(lst))(*[x.encode('utf-8') for x in lst])

    # =========================================================
    def thread_eventmon_loop(self,):
        # loop_counter = 0
        print('Event monitor loop started...')
        while(self.event_mon_loop_end == False):
            # print('Event monitor loop running...', loop_counter)

            # wr = self.lib_EventManager.wait_for_single_event(('_AOC_GENERAL_TRIGGER').encode('utf-8'),0)
            # if wr == 258:
            #     self.log("Waiting for '_AOC_GENERAL_TRIGGER' event returned TIMEOUT (258).", "orange")

            strEventList = ['_AOC_OK','_AOC_CALIB_DONE']
            wr = self.lib_EventManager.wait_for_multiple_events(len(strEventList),False,self.make_clist(strEventList),0)
            if wr == 258:
                # self.log('Waiting for multiple event list (' + str(len(strEventList)) + ') returned TIMEOUT (258).', 'orange')
                do_nothing = 1
            elif wr == 0:
                self.log(('Event named ' + strEventList[0] + ' has been signaled.'), 'blue')
            elif wr == 1:
                self.log(('Event named ' + strEventList[1] + ' has been signaled.'), 'blue')
                self.process_calibration(self.calib_mode)
                self.mycam.send_fifo_cmd(b'system set_cmd_matrix')

            # loop_counter = loop_counter + 1
            time.sleep(0.1)



# ==========================================================
# ==========================================================
if __name__ == "__main__":
    main(sys.argv)
