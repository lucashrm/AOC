#!/usr/bin/env python

import os
import time
import pdb

class Cam():
    ''' -----------------------------------------------------------------------
    Simple camera connection interface (via the fifo)

    ----------------------------------------------------------------------- '''

    # =========================================================================
    def __init__(self, fifo_dir="\\\\.\\pipe\\"):
        self.fifo_in  = fifo_dir + "ixon_fifo_in"
        self.fifo_out = fifo_dir + "ixon_fifo_out"
        self.connected = False
        self.streaming = False
        self.was_streaming = False

        self.MAX_FIFO_LEN = 256

        print(fifo_dir)

        # WINDOWS: USING os.path.exists(self.fifo_in) ON THE FIFO FILE BREAKS THE PIPE!!
#        if os.path.exists(self.fifo_in):
        try:
            self.cmd_fifo = os.open(self.fifo_in,os.O_RDWR)
            self.connected = True
            print("fifo in SUCCESSFULLY opened! (%s)" %(self.fifo_in))
            print("=====================================================")
        except:
            print("could not open the fifo in write mode (%s)" %(self.fifo_in))
            print("=====================================================")
#        else:
#            print("expected fifo does not exist (%s)" %(self.fifo_in))
#            print("=====================================================")

        # try:
        #     self.io_fifo = os.open(self.fifo_out,os.O_WRONLY)
        #     self.io_connected = True
        #     print("fifo in SUCCESSFULLY opened! (%s)" %(self.fifo_out))
        #     print("=====================================================")
        # except:
        #     print("could not open the fifo in write mode (%s)" %(self.fifo_out))
        #     print("=====================================================")

        # exposure time control: only a finite number of possibilities
        self.cam_tints = [0.00001, #0.00002, 0.00005,
                          0.00010, 0.00020, 0.00050,
                          0.00100, 0.00200, 0.00500,
                          0.01000, 0.02000, 0.05000,
                          0.10000, 0.20000, 0.50000]

        self.cam_tint = 0.0002 # fall back on 0.2 ms
        self.cam_itint = self.cam_tints.index(self.cam_tint)

        self.cam_fTemp = -999.

        if self.connected:
            self.cam_tint = 0.0002 # fall back on 0.2 ms
            #try:
            self.get_tint()
            #except:

            #pdb.set_trace()
            self.cam_itint = self.cam_tints.index(self.cam_tint)
            #self.set_tint(self.cam_tints[self.cam_itint])

    # =========================================================================
    def stream(self,):
        if self.connected:
            os.write(self.cmd_fifo,b"stream")
            os.fsync(self.cmd_fifo)
#            self.cmd_fifo.flush()
            self.streaming = True

    # =========================================================================
    def pause(self,):
        self.streaming = False
        if self.connected:
            os.write(self.cmd_fifo,b"abort")
            os.fsync(self.cmd_fifo)
#            self.cmd_fifo.flush()

    # =========================================================================
    def quit(self,):
        print('Cam: shutting down...')
        self.streaming = False
        if self.connected:
            os.write(self.cmd_fifo,b"quit")
            os.fsync(self.cmd_fifo)
#            self.cmd_fifo.flush()
        os.close(self.cmd_fifo)
        self.connected = False

    # =========================================================================
    def set_tint(self, tint):
        self.was_streaming = False
        if self.streaming is True:
            self.was_streaming = True
            self.pause()
        print('New exposure time command: %s' % (b"tint %.6f" % (tint,)))
        os.write(self.cmd_fifo,b"tint %.6f" % (tint,))
#        self.cmd_fifo.flush()
        if self.was_streaming:
            self.stream()
            self.streaming = True

    # =========================================================================
    def set_tempgoal(self, temp_goal):
        self.was_streaming = False
        if self.streaming is True:
            self.was_streaming = True
            self.pause()
        print('New exposure time command: %s' % (b"tempg %d" % (temp_goal,)))
        os.write(self.cmd_fifo,b"tempg %d" % (temp_goal,))
        os.fsync(self.cmd_fifo)
#        self.cmd_fifo.flush()

    def set_cooling(self, state):
        if state == 0:
            print('Setting camera cooling OFF')
            os.write(self.cmd_fifo,b"cooling 0")
            os.fsync(self.cmd_fifo)
        elif state == 1:
            print('Setting camera cooling ON')
            os.write(self.cmd_fifo,b"cooling 1")
            os.fsync(self.cmd_fifo)

    def get_fTemp(self,):
        if self.connected:
            os.write(self.cmd_fifo,b"fTemp?")
            os.fsync(self.cmd_fifo)
            tmp = float(os.read(self.cmd_fifo, self.MAX_FIFO_LEN))
            os.fsync(self.cmd_fifo)
            print("Current temperature:", tmp)
            self.cam_fTemp = tmp

            return self.cam_fTemp
        else:
            return -999

    # =========================================================================
    def get_tint(self):
        if self.connected:
            os.write(self.cmd_fifo,b"tint?")
            os.fsync(self.cmd_fifo)
            # self.cmd_fifo.flush()
            # with os.open(self.fifo_out, os.O_RDONLY) as fifo:
            #     tmp = float(fifo.read())
            #     print("integration time read:", tmp)

            tmp = float(os.read(self.cmd_fifo, self.MAX_FIFO_LEN))
            os.fsync(self.cmd_fifo)
            print("integration time read:", tmp)
            self.cam_tint = tmp # assign tint to current structure


    # =========================================================================
    def tint_dec(self,):
        self.cam_itint = max(self.cam_itint-1, 0)
        print("!!!!!!", self.cam_itint)
        self.cam_tint = self.cam_tints[self.cam_itint]
        if self.connected:
            self.set_tint(self.cam_tint)

    # =========================================================================
    def tint_inc(self,):
        self.cam_itint = min(self.cam_itint+1, len(self.cam_tints)-1)
        print("!!!!!!", self.cam_itint)
        self.cam_tint = self.cam_tints[self.cam_itint]
        if self.connected:
            self.set_tint(self.cam_tint)

# ==========================================================
# ==========================================================
if __name__ == "__main__":
    test = Cam(fifo_dir="\\\\.\\pipe\\")
    if test.connected:
        test.stream()
    else:
        print("Not connected to camera fifo")
