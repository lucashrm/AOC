#!/usr/bin/env python3

''' --------------------------------------------------------------------------

Master GUI for the control of CIAO, the Calern Imaging Adaptive Observatory.

Author: Frantz Martinache (frantz.martinache@oca.eu)
-------------------------------------------------------------------------- '''

from xaosim.QtMain import QtMain, QApplication
from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtCore import QThread, Qt
from PyQt5.QtWidgets import QLabel, QFileDialog, QWidget
from PyQt5.QtGui import QImage

import threading

import pyqtgraph as pg
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


guidef = "ciao_gui_LAb.ui"
try:
    shutil.copy(guidef, conf_dir+guidef)
except:
    print("Problem with config directory?")

zer_gui = "zernike_gui.ui"
try:
    shutil.copy(zer_gui, conf_dir+zer_gui)
except:
    print("Problem with config directory?")

from ciao_wfs import WFS
from ciao_wfc import *
from ciao_cam import Cam
from ciao_tel import Tel

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

        self.max_modes = 90
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

class MyWindow(QtGui.QMainWindow):
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

        # =========================================
        #               offsets
        # =========================================
        self.ostep = 0.5 # size of the offset step
        self.xoff  = 0.0 # keep track of total offset commands
        self.yoff  = 0.0 # sent to the WFS

        self.wfc_enabled        = False
        self.wfc_status_updated = True  # control GUI updates
        self.wfc_calibrating    = ""    # control variable

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
        self.overlay = pg.GraphItem()

        self.gView_SH_cam.addItem(self.imv_raw)
        self.gView_SH_cam.addItem(self.overlay)

        # -- second, the WFS data --
        self.gView_SH_info.hideAxis('left')
        self.gView_SH_info.hideAxis('bottom')

        self.imv_WFS = pg.ImageItem()
        self.gView_SH_info.addItem(self.imv_WFS)

        # -- the tip-tilt log plot --
        self.logplotx = self.gView_TT_log.plot([0,200],[0,0],
                                               pen=(0,255,0), name="ttx")
        self.logploty = self.gView_TT_log.plot([0,200],[0,0],
                                               pen=(0,0,255), name="tty")

        self.nmodes = 20
        self.bg_vx = 0.5 + np.arange(20)
        yA = np.linspace(0.0, 3.0, num=20)
        yB = np.linspace(0.0, 3.0, num=20)
        self.mode_graphA = pg.BarGraphItem(x=self.bg_vx, height=yA, width=0.3, brush='b')
        self.mode_graphB = pg.BarGraphItem(x=self.bg_vx+0.3, height=yB, width=0.3, brush='r')
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

        self.pB_WFS_off.clicked.connect(self.wfs_stop)
        self.pB_WFS_on.clicked.connect(self.wfs_start)
        self.pB_WFS_ref.clicked.connect(self.wfs_set_ref)
        self.pB_grid_set.clicked.connect(self.wfs_regrid)

        # --- camera / DM setup control ---
        self.pB_DM_shutdown.clicked.connect(self.alpao_shutdown)
        self.pB_DM_start.clicked.connect(self.alpao_start)

        #self.mycam = Cam(fifo_dir="/home/ciaodev/bin/")
        self.mycam = Cam(fifo_dir="\\\\.\\pipe\\")
        if self.mycam.connected:
            self.enable_cam_gui(state=True)

        self.pB_cam_stream.clicked.connect(self.mycam.stream)
        self.pB_cam_stop.clicked.connect(self.mycam.pause)

        self.pB_SetTempSetPoint.clicked.connect(self.SetCamSetPoint)
        self.tpB_Cooling.clicked.connect(self.ToggleCamCooling)
        self.pB_QueryCamTemperature.clicked.connect(self.QueryCamTemperature)
        self.tED_CurrentTemperature.setPlainText("N/A")
        # self.dspB_CamTemperature.valueChanged[float].connect(self.camSetPoint)

        self.pB_cam_tint_inc.clicked.connect(self.mycam.tint_inc)
        self.pB_cam_tint_dec.clicked.connect(self.mycam.tint_dec)
        self.pB_cam_shutdown.clicked.connect(self.ixon_shutdown)

        self.chB_TCSOffload.stateChanged[int].connect(self.TCSOffload)

        self.chB_ShowZerWindow.stateChanged[int].connect(self.show_zer_window)

        # ==============================================
        self.show()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.refresh_all)
        self.timer.start(10)

        # --- wavefront control callbacks ---
        self.bDoMeasure = False
        self.bCloseLoop = False


        self.TT_WFC  = TT_WFC(simu=self.simu)               # tip-tilt wavefront control
        self.ZON_WFC = ZON_WFC(simu=self.simu)              # Zonal wavefront contrl
        # self.ZER_WFC = ZER_WFC(iz0=4, iz1=10) # Zernike wavefront control
        self.ZER_WFC = ZER_WFC(iz0=2, iz1=10, simu=self.simu) # Zernike wavefront control

        print('Number of ZER_WFC instances: %d' % sys.getrefcount(self.ZER_WFC))

        self.TT_WFC.reload_cal(fname=conf_dir+"TT_CAL.fits")
        self.ZON_WFC.reload_cal(fname=conf_dir+"ZON_CAL.fits")
        self.ZER_WFC.reload_cal(fname=conf_dir+"ZER_CAL.fits")

        self.dspB_zo_a0.setValue(self.ZON_WFC.a0)
        self.dspB_zer_a0.setValue(self.ZER_WFC.a0)
        self.dspB_tt_a0.setValue(self.TT_WFC.a0)

        self.dspB_zer_nz.valueChanged[int].connect(self.zer_update_nz)

        self.ZER_WFC.verbose = True
        self.ZER_WFC.nmodes = self.dspB_zer_nz.value() - 1

        self.ZON_WFC.nmodes = 97
        self.dms = 11       # NUMBER OF ACTUATORS ACROSS (SQUARE)

        # for sensor in ["tt", "zer", "zo"]:
        #     exec("self.pB_%s_abort.clicked.connect(self.%s_abort)" % (sensor, sensor))
        #     exec("self.pB_%s_measure.clicked.connect(self.%s_measure)" % (sensor, sensor))
        #     exec("self.pB_%s_recal.clicked.connect(self.%s_calib)" % (sensor, sensor))
        #     exec("self.pB_%s_reset.clicked.connect(self.%s_reset_correc)" % (sensor, sensor))
        #     exec("self.dspB_%s_gain.valueChanged[float].connect(self.%s_update_gain)" % (sensor, sensor))
        #     exec("self.%s_update_gain()" % (sensor,))
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

        # ====================================================
        # live tip-tilt offsets
        # ====================================================
        self.pB_tt_dx_p.clicked.connect(self.ttx_offset_p)
        self.pB_tt_dx_m.clicked.connect(self.ttx_offset_m)
        self.pB_tt_dy_p.clicked.connect(self.tty_offset_p)
        self.pB_tt_dy_m.clicked.connect(self.tty_offset_m)
        # ====================================================

        self.tel = Tel()
        self.tel_mon_start()
        self.tel.offset_angle_deg = 75.

        self.zer_win = None

        self.zer_update_nz()

        self.derotator_angle_deg = 0.
        self.derotator_on = 0

        username = 'neo'
        pwd = 'lagrange'
        client = mqtt.Client("JOVIAL")
        client.on_connect = self.on_connect
        client.on_message = self.on_message
        client.username_pw_set(username, password=pwd)
        client.connect_async('10.150.1.100',1883)
        client.loop_start()

    def on_connect(self, client, userdata, flags, rc):
        # print('Connected with result code {0}'.format(str(rc)))
        client.subscribe("JOVIAL/derotator")

    def on_message(self, client, userdata, msg):
        self.derotator_angle_deg = float(json.loads(msg.payload)['theta_c'])
        self.derotator_on = float(json.loads(msg.payload)['on'])
        self.tel.offset_angle_deg = 75.+self.derotator_angle_deg*self.derotator_on

        # print('Message received >> ' + msg.topic + ' ' + str(msg.payload), str(self.derotator_angle_deg))
        # print('Angle sent to tel >> ' ,str(self.tel.offset_angle_deg))

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
        else:
            self.log("Slope+Illum (TT) data logging to file disabled", "red")

    # =========================================================
    def set_mode_logging(self,):
        self.ZER_WFC.bLogCmds = self.chB_LogModes.isChecked()
        self.ZER_WFC.bLogErrors = self.ZER_WFC.bLogCmds
        self.ZER_WFC.bLogSlopes = self.ZER_WFC.bLogCmds

        if self.ZER_WFC.bLogCmds == True:
            self.log("Modes (ZER) data logging to file enabled", "green")
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
        x0 = self.spB_grid_x0.value()
        y0 = self.spB_grid_y0.value()
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
        self.wfs.set_thread_params()
        self.wfs.keepgoing = False
        self.enable_wfs_gui(state=True)
        # self.ZER_WFC.DM.dmmask = self.wfs.SH_phot
        # fits.writeto(r'C:\Home\Downloads\self.ZER_WFC.DM.dmmask.fits',self.ZER_WFC.DM.dmmask,overwrite=True)
        self.wfs.calc_SH_data(ref=True)
        self.log("Current pos. set as reference", "blue")

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

        for i, ptile in enumerate(pt_levels):
            msg += "p-tile %3d = %8.2f\n" % (ptile, pt_values[i])
        msg += "\ntt_xy = (%+.2f,%+.2f)\n" % (self.wfs.ttx_mean, self.wfs.tty_mean)
        msg += "</pre>"
        self.lbl_stats.setText(msg)

    # =========================================================
    def redraw_SH_grid(self):
        ''' ----------------------------------------------------------
        Displays a grid overlaid the camera image.
        ---------------------------------------------------------- '''
        self.SH_x0 = self.spB_grid_x0.value()
        self.SH_y0 = self.spB_grid_y0.value()
        self.SH_dx = self.spB_grid_dx.value()
        self.SH_dy = self.spB_grid_dy.value()

        isz = 128

        # nx = int((isz - self.SH_x0) / self.SH_dx + 2)
        nx = self.dms + 1
        xx = self.SH_x0 + np.arange(nx-1) * self.SH_dx
        if (isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, isz)

        # ny = int((isz - self.SH_y0) / self.SH_dy + 2)
        ny = self.dms + 1
        yy = self.SH_y0 + np.arange(ny-1) * self.SH_dy
        if (isz - yy[-1] > self.SH_dy/2+2):
            yy = np.append(yy, isz)

        print('nx,ny : ',nx,ny)

        posx = np.zeros((2 * nx, 2), dtype=np.int_)
        for i, myx in enumerate(xx):
            posx[2*i, 0] = myx
            posx[2*i, 1] = yy[0]
            posx[2*i+1, 0] = myx
            posx[2*i+1, 1] = min(yy[-1], isz-1) # -1 for display

        posy = np.zeros((2 * ny, 2), dtype=np.int_)
        for i, myy in enumerate(yy):
            posy[2*i,   1] = myy
            posy[2*i,   0] = xx[0]
            posy[2*i+1, 1] = myy
            posy[2*i+1, 0] = min(xx[-1], isz-1)

        pos = np.append(posy, posx, axis=0)
        adj = np.arange(2*(nx+ny)).reshape(nx+ny, 2)

        self.SH_xx = xx
        self.SH_yy = yy

        if self.chB_show_grid.isChecked():
            self.overlay.setData(pos=pos, adj=adj,
                                 pen=(255,0,0,255), size=0,
                                 #symbolPen=None, symbolBrush=None,
                                 pxMode=False)
        else:
            self.overlay.setData(pos=pos, adj=adj, pen=None, size=1,
                                 #symbolPen=None, symbolBrush=None,
                                 pxMode=False)
    # =========================================================
    def refresh_all(self):
        ''' ----------------------------------------------------------
        Refresh the display

        Independently from the "under the hood" engine, the display
        gets refreshed every now and then, to give visual feedback
        to the user.
        ---------------------------------------------------------- '''
        if self.cam_counter < self.shm_cam.get_counter():
            self.data_cam = self.shm_cam.get_data(check=False, reform=True)
            self.cam_counter = self.shm_cam.get_counter()
            self.imv_raw.setImage(arr2im(self.data_cam.T,
                                         vmin=self.vmin, vmax=self.vmax,
                                         pwr=self.pwr,
                                         cmap=self.mycmap), border=1)
        # --------
#        print(self.wfs.ttx_log)
        self.logplotx.setData(self.wfs.ttx_log)
        self.logploty.setData(self.wfs.tty_log)

        if self.bDoMeasure == True:

            # self.bTTLoop = False
            if self.ZON_WFC.keepgoing == True:
                mode_amp = np.zeros(self.ZON_WFC.nmodes)
                mode_amp_res = np.zeros(self.ZON_WFC.nmodes)
                mode_amp_res[0:self.ZON_WFC.zer_mode_amp.size] = self.ZON_WFC.zer_mode_amp_resid[0:]

                self.mode_graphA.setOpts(x=self.bg_vx,height=mode_amp)
                self.mode_graphB.setOpts(x=self.bg_vx+0.3,height=mode_amp_res)
            elif self.ZER_WFC.keepgoing == True:
            # print('CM: nmodes=%d' % self.ZER_WFC.nmodes)
                mode_amp = np.zeros(self.ZER_WFC.nmodes)
                mode_amp_res = np.zeros(self.ZER_WFC.nmodes)
                # mode_amp[2:self.ZER_WFC.zer_mode_amp.size+2] = self.ZER_WFC.zer_mode_amp
                if self.TT_WFC.keepgoing == True:
                    mode_amp[0:2] = self.TT_WFC.zer_mode_amp[0:2]
                    mode_amp[2:] = self.ZER_WFC.zer_mode_amp[2:]
                    mode_amp_res[0:2] = self.TT_WFC.zer_mode_amp_resid[:]
                    mode_amp_res[2:self.ZER_WFC.zer_mode_amp_resid.size] = self.ZER_WFC.zer_mode_amp_resid[2:]
                else:
                    mode_amp[0:self.ZER_WFC.zer_mode_amp.size] = self.ZER_WFC.zer_mode_amp[0:]
                    mode_amp_res[0:self.ZER_WFC.zer_mode_amp_resid.size] = self.ZER_WFC.zer_mode_amp_resid

                self.mode_graphA.setOpts(x=self.bg_vx,height=mode_amp)
                self.mode_graphB.setOpts(x=self.bg_vx+0.3,height=mode_amp_res)
            # self.mode_graphB.setOpts(height=self.ZER_WFC.zer_mode_amp_resid)

        self.refresh_wfs_display()
        self.refresh_stats()

        # print('Frame drop vector >> ',self.wfs.vFrameDrop)


        # if (self.wfs.tavg_count != 0):
        #     pixScale = 10.
        #     self.tel.ttx_tavg = self.wfs.ttx_tavg / self.wfs.tavg_count * pixScale
        #     self.tel.tty_tavg = self.wfs.tty_tavg / self.wfs.tavg_count * pixScale
        #     self.wfs.ttx_tavg = 0.
        #     self.wfs.tty_tavg = 0.
        #     self.wfs.tavg_count = 0.

        #     # self.tel.ttx_cor = self.TT_WFC.ttx_mode * 10.
        #     # self.tel.tty_cor = self.TT_WFC.tty_mode * 10.
        #     # self.TT_WFC.ttx_mode = 0.
        #     # self.TT_WFC.tty_mode = 0.

        # -------
        if self.wfc_status_updated: # only update GUI if needed
            self.wfc_status_updated = False

            if self.wfc_enabled:
                self.enable_wfc()
            else:
                self.disable_wfc()

            if self.TT_WFC.calib_on or self.TT_WFC.loop_on:
                self.disable_wfc(subsystem="tt")
            else:
                self.enable_wfc(subsystem="tt")

            if self.ZER_WFC.calib_on or self.ZER_WFC.loop_on:
                self.disable_wfc(subsystem="zer")
            else:
                self.enable_wfc(subsystem="zer")

            if self.ZON_WFC.calib_on or self.ZON_WFC.loop_on:
                self.disable_wfc(subsystem="zo")
            else:
                self.enable_wfc(subsystem="zo")


        # -------
        if self.wfc_calibrating == "TT": # query thread status
            if self.ttcal_Thread.isFinished():
                self.enable_wfc(subsystem="tt")
                self.log("TT calibration is over", "green")
                self.wfc_calibrating = ""

        # -------
        if self.wfc_calibrating == "ZER": # query thread status
            if self.zercal_Thread.isFinished():
                self.enable_wfc(subsystem="zer")
                self.log("ZER calibration is over", "green")
                self.wfc_calibrating = ""

        # -------
        if self.wfc_calibrating == "ZON": # query thread status
            if self.zocal_Thread.isFinished():
                self.enable_wfc(subsystem="zo")
                self.log("ZON calibration is over", "green")
                self.wfc_calibrating = ""


    # =========================================================
    def refresh_wfs_display(self):
        ''' ----------------------------------------------------------
        Specific display refresh for the WFS data (ttxy + phot)
        ---------------------------------------------------------- '''
        self.imv_WFS.setImage(
            arr2im(np.concatenate(self.wfs.SH_comb),
                   vmin=-5, vmax=5, pwr=1, cmap=self.mycmap),
            border=1)

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
        self.ZER_WFC.keepgoing = False
        nav = self.spB_nav.value()
        a0 = self.dspB_zer_a0.value()
        self.ZER_WFC.nav = nav
        self.disable_wfc(subsystem="zer")
        self.log("Zernike loop calibration starts", "blue")
        self.log("cal: nav=%d, a0=%.2f" % (nav, a0), "blue")
        self.zercal_Thread = GenericThread(self.ZER_WFC.calibrate,
                                           a0=a0, reform=False)
        self.wfc_calibrating = "ZER"
        self.zercal_Thread.start()

        print('self.ZER_WFC.nmodes: %d' % self.ZER_WFC.nmodes)

        self.bg_vx = 0.5 + np.arange(self.ZER_WFC.nmodes)
        yA = np.linspace(0.0, 3.0, num=self.ZER_WFC.nmodes)
        yB = np.linspace(0.0, 3.0, num=self.ZER_WFC.nmodes)
        self.mode_graphA.setOpts(x=self.bg_vx,height=yA)
        self.mode_graphB.setOpts(x=self.bg_vx,height=yB)

    # =========================================================
    def zer_update_gain(self):
        gain = self.dspB_zer_gain.value()
        self.ZER_WFC.gain = gain
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
            self.bCloseLoop = True
            self.ZER_WFC.bApplyCorrection = True
            print('Close loop enabled')
        else:
            self.bCloseLoop = False
            self.ZER_WFC.bApplyCorrection = False
            print('Close loop disabled')

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
        self.ZON_WFC.keepgoing = False
        nav = self.spB_nav.value()
        a0 = self.dspB_zo_a0.value()
        self.ZON_WFC.nav = nav
        self.disable_wfc(subsystem="zo")
        self.log("ZON loop calibration starts", "blue")
        self.log("cal: nav=%d, a0=%.2f" % (nav, a0), "blue")
        self.zocal_Thread = GenericThread(self.ZON_WFC.calibrate,
                                          a0=a0, reform=False)
        self.wfc_calibrating = "ZON"
        self.zocal_Thread.start()

        self.bg_vx = 0.5 + np.arange(self.ZON_WFC.nmodes)
        yA = np.linspace(0.0, 3.0, num=self.ZON_WFC.nmodes)
        yB = np.linspace(0.0, 3.0, num=self.ZON_WFC.nmodes)
        self.mode_graphA.setOpts(x=self.bg_vx,height=yA)
        self.mode_graphB.setOpts(x=self.bg_vx,height=yB)

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
                pixScale = 10.
                if self.wfs.tavg_count != 0.:
                    self.tel.ttx_tavg = self.wfs.ttx_tavg / self.wfs.tavg_count * self.OffloadPixScale
                    self.tel.tty_tavg = self.wfs.tty_tavg / self.wfs.tavg_count * self.OffloadPixScale
                    self.wfs.ttx_tavg = 0.
                    self.wfs.tty_tavg = 0.
                    self.wfs.tavg_count = 0.

                    print('TTX/Y last avg: %.3e %.3e' % (self.tel.ttx_tavg,self.tel.tty_tavg))
                    print('Frame drop vector: ',self.wfs.vFrameDrop)
                    self.wfs.vFrameDrop = np.zeros(self.wfs.vFrameDrop.shape[0])

                    self.tel.offload_tiptilt(self.tel.ttx_tavg, self.tel.tty_tavg, self.tel.offset_angle_deg)


    def TCSOffload(self,):
        if self.chB_TCSOffload.isChecked():
            self.tel.bEnableExtTTOffload = True
            print('TCS offload enabled')
        else:
            self.tel.bEnableExtTTOffload = False
            print('TCS offload disabled')

# ==========================================================
# ==========================================================
if __name__ == "__main__":
    main(sys.argv)
