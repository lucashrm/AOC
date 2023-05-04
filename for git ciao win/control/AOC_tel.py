#!/usr/bin/env python

import time
import socket
import threading
import numpy as np

class Tel():
    ''' -----------------------------------------------------------------------
    Simple telescope connection interface

    Primarily used to get pointing information from the telescope
    and to send tip-tilt offload commands
    ----------------------------------------------------------------------- '''

    # =========================================================================
    def __init__(self,):
        self.TCS_IP   = "10.150.10.18"
        self.TCS_PORT = 3737
        # self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        # self.sock.settimeout(1)
        self.connected = False
        self.keepgoing = False
        # self.mon_period_sec = 30.
        self.ttx_cor = 0.
        self.tty_cor = 0.
        self.ttx_tavg = 0.
        self.tty_tavg = 0.
        self.bEnableExtTTOffload = False

    # =========================================================================
    def connect(self,):
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self.sock.settimeout(1)
            self.sock.connect((self.TCS_IP, self.TCS_PORT))
            self.sock.getblocking()
            self.connected = True
            # print("Socket connected!")
        except:
            print("Socket connection refused")
        return self.connected

    # =========================================================================
    def close(self,):
        self.sock.close()
        self.connected = False

    # # =========================================================================
    # def get_hour_angle(self):
    #     if self.connect() != True:
    #         print('WARNING! Could not connect to TCS! Nothing done.')
    #         return False
    #     self.sock.send(b'#TCS_GET_POINTING_INFO:?\r')
    #     self.connected = True
    #     reply = self.sock.recv(1024)[1:-1].split(b' ')
    #     self.tcs_alpha = float(reply[0])
    #     self.tcs_delta = float(reply[1])
    #     self.tcs_azi = float(reply[2])
    #     self.tcs_alt = float(reply[3])
    #     self.tcs_TS = float(reply[4])
    #     self.tcs_HA = float(reply[5])
    #     self.tcs_focus = float(reply[6])
    #     self.sock.close()
    #     self.connected = False
    #     return True

    # =========================================================================
    def get_pointing_info(self):
        if self.connect() != True:
            print('WARNING! Could not connect to TCS! Nothing done.')
            self.tcs_HA = 0.
            return False
        self.sock.send(b'#TCS_GET_POINTING_INFO:?\r')
        self.connected = True
        try:
            reply = self.sock.recv(1024)[1:-1].split(b' ')
        except:
            print('Failed to receive data')
            return False
        self.tcs_alpha = float(reply[0])
        self.tcs_delta = float(reply[1])
        self.tcs_azi = float(reply[2])
        self.tcs_alt = float(reply[3])
        self.tcs_TS = float(reply[4])
        self.tcs_HA = float(reply[5])
        self.tcs_focus = float(reply[6])
        self.close()
        self.connected = False
        # print('ALPHA: %.6f  DELTA: %.6f  AZI: %.6f  ALT: %.6f  TS: %.6f  HA: %.6f  FOCUS: %.6f' %
        #     (float(reply[0]),float(reply[1]),float(reply[2]),float(reply[3]),float(reply[4]),float(reply[5]),float(reply[6])))
        return True

    # =========================================================================
    def set_guiding_offsets(self, ra=0.0, dec=0.0):
        if self.connect() != True:
            print('WARNING! Could not connect to TCS! Nothing done.')
            return False
        self.connected = True
        cmd = b'#TCS_SET_GUIDING_OFFSETS: %.2f %.2f?\r' % (ra, dec)
        self.sock.send(cmd)
        print(cmd)
        reply = self.sock.recv(1024)[1:-1].split(b' ')
        self.close()
        print(reply)
        return True

    # =========================================================================
    def move_focus_relative(self, inc=0.0):
        if self.connect() != True:
            print('WARNING! Could not connect to TCS! Nothing done.')
            return False
        self.connected = True
        cmd = b'#TCS_MOVE_FOCUS_RELATIVE: %.1f?\r' % (inc)
        self.sock.send(cmd)
        print(cmd)
        reply = self.sock.recv(1024)[1:-1].split(b' ')
        self.close()
        print(reply)
        return True

    # =========================================================================
    def offload_tiptilt(self, ttx, tty, h0_deg=0.0):
        # if not self.connected:
        #     self.connect()

        self.get_pointing_info()

        h0_rad = h0_deg*np.pi/180.0             # mystery angle
        houra = self.tcs_HA + h0_rad        # hour angle in rad
        # print('Correction angle in degrees >> ', houra*180/np.pi, h0_deg)

        # myrot = np.array([[-np.sin(houra), -np.cos(houra)],
        #                   [-np.cos(houra), np.sin(houra)]])
        myrot = np.array([[ np.cos(houra), np.sin(houra)],
                          [-np.sin(houra), np.cos(houra)]])

        mycmd = np.dot(myrot, [-ttx, tty])   # in arcseconds
        # print('HA: %.6f' % self.tcs_HA)
        print(mycmd)
        if self.bEnableExtTTOffload == True:
            # print('>> Sending offsets to TCS')
            self.set_guiding_offsets(ra=mycmd[0], dec=mycmd[1])


# ==========================================================
# ==========================================================
if __name__ == "__main__":
    test = Tel()
    test.connect()
    if test.connected:
        print("Connected to TCS!")
    else:
        print("Not connected to TCS")

    test.close()
