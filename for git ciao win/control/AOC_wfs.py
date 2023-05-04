#!/usr/bin/env python3

import sys
import numpy as np
from xaosim.shmlib import shm
import datetime
from datetime import datetime
import time
#import fitsio
#from fitsio import FITS,FITSHDR
import astropy.io.fits as fits

import os
import ctypes
import _ctypes

from scipy.ndimage import center_of_mass as comass
# would be worth giving this function a shot!

import multiprocessing as mp

# !!! multi-CPU not working yet !!!

ciao_home = os.getenv('CIAO_HOME')

def ext_centroid_position_cell(args, **kwargs):
    return WFS.centroid_position_cell(*args, **kwargs)

def ext_max_cell(args, **kwargs):
    return WFS.max_cell(*args, **kwargs)

def centroid_position(img):
    ''' --------------------------------------------------------------
    Wrapper function to make it easy to change different algorithms
    -------------------------------------------------------------- '''
    #return comass(img) # the scipy center_of_mass?
    return maximum_position(img) # another possibility
    #return centroid_position_0(img)

def centroid_position_0(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    xxc = np.arange(img.shape[1])
    yyc = np.arange(img.shape[0])
    mprofx = img.mean(0)
    mprofy = img.mean(1)
    mprofx -= mprofx.min()
    mprofy -= mprofy.min()
    denomx = np.sum(mprofx)
    denomy = np.sum(mprofy)

    if denomx > 1:
        xc = np.sum(mprofx * xxc) / denomx
    else:
        xc = 0.0
    if denomy > 1:
        yc = np.sum(mprofy * yyc) / denomy
    else:
        yc = 0.0

    print('xc,yc >> %f %f' % (xc,yc))
    return (yc, xc)

def maximum_position(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    return np.unravel_index(np.argmax(img, axis=None), img.shape)

# =============================================================================
# =============================================================================

class WFS():
    # =========================================================
    def __init__(self, shmf="R:\\SHcam.im.shm", ncpu=False):
        ''' -------------------------------------------------------
        Default Wavefront Sensor (WFS) class constructor.

        Parameters:
        ----------
        - shmf: shared memory data structure where to read images
        - ncpu: (integer) number of cpus to use for computation
        ------------------------------------------------------- '''
        # -------------------------
        # to be provided externally
        # -------------------------
        self.SH_x0, self.SH_y0 = 0, 0   # properties of the grid
        self.SH_dx, self.SH_dy = 12.8, 12.8 # properties of the grid
        # self.SH_dx, self.SH_dy = 11.6, 11.6 # properties of the grid
        self.threshold = 180.0
        self.nCells = 0

        # -------------------------------
        # external C++ multi-threaded DLL
        # -------------------------------
        self.bUseMT = False
        self.nThreads = 8
        self.bThreadInitialized = False

        self.lib = ctypes.CDLL(r'C:\Users\lucas\Documents\STAGE\CIAO_Win\DEV_CIAO\x64\Release\DLLfromPython.dll')
        self.lib.calc_SH_data.restype = ctypes.c_int
        self.lib.calc_SH_data.argtypes = [ctypes.c_int, #
            ctypes.c_int,                               #
            np.ctypeslib.ndpointer(dtype=np.int32),     # vImg
            ctypes.c_float,                             # bckgd
            ctypes.c_int,                               #
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.float32),
            np.ctypeslib.ndpointer(dtype=np.float32),
            np.ctypeslib.ndpointer(dtype=np.float32)]

        self.lib.StartThread_ThreadManager.restype = ctypes.c_int
        self.lib.StartThread_ThreadManager.argtypes = [
            ctypes.c_int,                               # nThreads
            ctypes.c_int,                               # sizex
            ctypes.c_int,                               # sizey
            ctypes.c_int,                               # nx
            ctypes.c_int,                               # ny
            ]

        self.lib.set_thread_params.restype = ctypes.c_int
        self.lib.set_thread_params.argtypes = [
            ctypes.c_int,                               # threadID
            ctypes.c_int,                               # nsub
            np.ctypeslib.ndpointer(dtype=np.int32),     # v_xx
            np.ctypeslib.ndpointer(dtype=np.int32),     # v_yy
            np.ctypeslib.ndpointer(dtype=np.int32),     # vIndex_i
            np.ctypeslib.ndpointer(dtype=np.int32),     # vIndex_j
            np.ctypeslib.ndpointer(dtype=np.float32),   # SH_xtmp
            np.ctypeslib.ndpointer(dtype=np.float32),   # SH_ytmp
            np.ctypeslib.ndpointer(dtype=np.float32)    # SH_phot
            ]

        self.lib.trigger_calc_SH_data_mt.restype = ctypes.c_int
        self.lib.trigger_calc_SH_data_mt.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32),
            ctypes.c_float
            ]

        # -----------------------------
        # shared memory data structures
        # -----------------------------
        self.shm_im = shm(shmf)                  # shared mem. image
        self.isz = self.shm_im.mtdata['size'][1] # should be 128

        verbose = False
        if verbose:
            print(shmf)
            print("self.shm_im.mtdata['imname'] = %s" % (self.shm_im.mtdata['imname']))
            print("self.shm_im.mtdata['bimname'] = %s" % (self.shm_im.mtdata['bimname']))
            print("self.shm_im.mtdata['version'] = %s" % (self.shm_im.mtdata['version']))
            print("self.shm_im.mtdata['bversion'] = %s" % (self.shm_im.mtdata['bversion']))
            print("self.shm_im.mtdata['naxis'] = %d" % (self.shm_im.mtdata['naxis']))
            print("self.shm_im.mtdata['size'] = %d %d" % (self.shm_im.mtdata['size'][1],self.shm_im.mtdata['size'][2]))
            print(self.shm_im.mtdata['size'])

        self.sizex = self.shm_im.mtdata['size'][1]
        self.sizey = self.shm_im.mtdata['size'][2]
        print('SH image size x/y: %d %d' % (self.sizex,self.sizey,))

        self.im_cnt = self.shm_im.get_counter()  # image counter
        print("self.im_cnt: %d" % self.im_cnt)
        self.update_cells()
        self.define_SH_data()

        self.ttx_mean = 0.0
        self.tty_mean = 0.0
        self.log_len = 200 # length of the tip-tilt log
        self.ttx_log = []
        self.tty_log = []

        self.ttx_tavg = 0.0
        self.tty_tavg = 0.0
        self.tavg_count = 0.0

        self.max_frame_drop = 100
        # KEEP self.vFrameDrop[0] FOR DROP INTERVALS GREATER THAN 20 FRAMES
        self.vFrameDrop = np.zeros(self.max_frame_drop,dtype=np.float32)

        self.prev_im_cnt = 0

        # ------------------
        # control structures
        # ------------------
        self.keepgoing = False # close-loop flag

        self.mproc = False # multi-processor option
        if ncpu is not False:
            ncpumax = mp.cpu_count()
            print("Multi processor option selected")
            print("%d CPUs available: %d requested" % (ncpumax, ncpu))
            self.mproc = True
            self.ncpu = ncpu
            self.pool = mp.Pool(ncpu) # split computation on CPU pool


    # =========================================================
    def update_cells(self):
        ''' -------------------------------------------------------
        Updates the SH cells positions for new SH grid parameters
        ------------------------------------------------------- '''
        print('self.SH_dx/dy: ',self.SH_dx,self.SH_dy)

        nx = int((self.isz - self.SH_x0) / self.SH_dx + 2)
        xx = np.round(self.SH_x0+np.arange(nx-1)*self.SH_dx).astype(int)
        if (self.isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, self.isz)
        self.SH_xx = xx

        print('nx, xx, self.SH_xx >> ',nx, xx, self.SH_xx)

        ny = int((self.isz - self.SH_y0) / self.SH_dy + 2)
        yy = np.round(self.SH_y0+np.arange(ny-1)*self.SH_dy).astype(int)
        if (self.isz - yy[-1] > self.SH_dy/2+2):
            yy = np.append(yy, self.isz)
        self.SH_yy = yy

        # -------------------------------------------------------------
        # with these coordinates computed, we can now define the slices
        # of the image that will correspond to the different cells
        # Hopefully, these slices will accelerate the computation !!
        # -------------------------------------------------------------
        self.SH_cells = []
        for jj in range(ny-2):
            for ii in range(nx-2):
                self.SH_cells.append([slice(yy[jj], yy[jj+1]), slice(xx[ii], xx[ii+1])])

    # =========================================================
    def define_SH_data(self):
        ''' -------------------------------------------------------
        Creates empty arrays that will be used to store the information
        extracted from the analysis of the SH images: photometry, slope
        in the x and y directions, position of the reference

        Creates the shared memory data structures that will be used
        to communicate this information with other processes:

        - shm_phot_inst: current illumination of the apertures
        - shm_comb     : xy slopes & photometry combo

        Currently not used: shm_SNR
        ------------------------------------------------------- '''
        print('Entering WFS::define_SH_data')
#        print(self.SH_xx)
#        print("xx,yy >> %d %d" %(self.SH_xx , self.SH_yy))
        xx, yy   = self.SH_xx , self.SH_yy  # position of the SH cell edges
#        print('xx,yy >> ')
#        print(xx)
#        print(yy)
        ncx, ncy = xx.size - 1, yy.size - 1 # number of cells
#        print(xx.size)
#        print(yy.size)
        # -------------
        # empty arrays
        # -------------
        # if (self.bThreadInitialized == False):
        self.SH_phot = np.zeros((ncy, ncx),dtype=np.float32) # photometry
        self.SH_xtmp = np.zeros((ncy, ncx),dtype=np.float32)
        self.SH_ytmp = np.zeros((ncy, ncx),dtype=np.float32)
        self.SH_xslp = np.zeros((ncy, ncx),dtype=np.float32) # x-slope
        self.SH_yslp = np.zeros((ncy, ncx),dtype=np.float32) # y-slope
        self.SH_xref = np.zeros((ncy, ncx),dtype=np.float32) # x-reference
        self.SH_yref = np.zeros((ncy, ncx),dtype=np.float32) # y-reference
        self.SH_comb = np.zeros((3, ncy, ncx),dtype=np.float32) # xy slopes + phot combined

        self.SH_phot[:] = 0.
        self.SH_xtmp[:] = 0.
        self.SH_ytmp[:] = 0.
        # # ----------------------------
        # # additional utility arrays
        # # ----------------------------
        # self.SH_xtmp = np.zeros((ncy, ncx),dtype=np.float32)
        # self.SH_ytmp = np.zeros((ncy, ncx),dtype=np.float32)
        print(type(self.SH_xtmp[0,0]))

        print('ncx,ncy : %d %d' % (ncx, ncy))

        self.shm_SNR = shm(
            r'R:\SNR.im.shm', data=self.SH_phot, verbose=False) # experiment!
        self.shm_phot_inst = shm(
            r'R:\phot_inst.im.shm', data=self.SH_phot, verbose=False)

        self.shm_comb = shm(r'R:\comb.im.shm',data=self.SH_comb,verbose=False)

        for jj in range(ncy):
            #y0, y1 = int(np.round(yy[jj])), int(np.round(yy[jj+1]))
            y0, y1 = yy[jj], yy[jj+1]
            for ii in range(ncx):
                #x0, x1 = int(np.round(xx[ii])), int(np.round(xx[ii+1]))
                x0, x1 = xx[ii], xx[ii+1]
                self.SH_xref[jj,ii] = np.round(0.5 * (x1 - x0)) + 1
                self.SH_yref[jj,ii] = np.round(0.5 * (y1 - y0)) + 1

        fits.writeto("SHrefx.fits", self.SH_xref, overwrite=True)
        fits.writeto("SHrefy.fits", self.SH_yref, overwrite=True)

        if (self.bThreadInitialized == False):
            # START THE THREADS NOW (IN IDLE STATE)
            retval = self.lib.StartThread_ThreadManager(self.nThreads,self.sizex,self.sizey,ncx,ncy)
            if (retval != 0):
                print('ERROR! Could not start the threads!')
            else:
                print('SUCCESS: %d threads started!' % (self.nThreads))

        self.set_thread_params()


    # =========================================================
    def update_grid(self, **kwargs):
        ''' -------------------------------------------------------
        update the parameters that define the SH cells
        ------------------------------------------------------- '''
        valid_keys = {"x0": "SH_x0", "y0": "SH_y0",
                      "dx": "SH_dx", "dy": "SH_dy",
                      "i0": "threshold"}
        if kwargs is not None:
            for key, value in kwargs.items():
                if key in valid_keys:
                    setattr(self, valid_keys[key], value)
                else:
                    print("%s: invalid key" % (key,))
                    print("valid keys are: " +
                          len(valid_keys.keys()) * "%s " % \
                          tuple(valid_keys.keys()))
        else:
            print("No parameter was updated")

    # =========================================================
    def centroid_position_cell(self, aslice):
        ''' -------------------------------------------------------
        Class specific wrapper function that returns the centroid
        for a given slice of the live image
        ------------------------------------------------------- '''
        return centroid_position(self.live_img[aslice])

    # =========================================================
    def max_cell(self, aslice):
        ''' -------------------------------------------------------
        Another such wrapper that should eventually be merged with
        the centroid since both info can easily be aggregated in
        one go.
        ------------------------------------------------------- '''
        return self.live_img[aslice].max()

    # =========================================================
    def calc_SH_data_new(self, data=None, ref=True):
        ''' -------------------------------------------------------
        Calculates the Shack-Hartman data for the current image.

        Parameters:
        - ref: if True -> current image taken as new reference (bool)
        - data: if not None, info computed from the provided image
        ------------------------------------------------------- '''
        ncx, ncy = self.SH_xx.size - 1, self.SH_yy.size - 1

        if data is None:
            self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True,sleepT=0)
        else:
            self.live_img = data

        self.im_cnt = self.shm_im.get_counter()  # update image counter
        self.live_img[self.live_img < self.threshold] = 0.0#self.threshold
        #self.live_img -= self.threshold
        #self.live_img[self.live_img <=0] = 0.0

        # =======================================================
        # trying to write this the right way now...
        # optional self.pool.map(...)
        # despite my efforts, the pool trick doesn't seem to work with python2.7

        if self.mproc is True:
            # -------------------------------
            # multi-processor scenario (TBC)
            # -------------------------------
            pool = mp.Pool(self.ncpu) # split computation on CPU pool
            cen_list = np.array(list(pool.map(ext_centroid_position_cell, zip([self]*len(self.SH_cells), self.SH_cells))))
            max_list = np.array(list(pool.map(ext_max_cell, zip([self]*len(self.SH_cells), self.SH_cells))))
            pool.close()
        else:
            # -------------------------------
            # this is the single CPU scenario
            # -------------------------------

            # get the list of centroid positions
            cen_list = np.array(list(map(self.centroid_position_cell, self.SH_cells)))
            # and get the list of maximum values
            max_list = np.array(list(map(self.max_cell, self.SH_cells)))

        self.SH_phot = max_list.reshape((ncy,ncx)).astype('float')
        self.SH_ytmp = cen_list[:,0].reshape((ncy, ncx))
        self.SH_xtmp = cen_list[:,1].reshape((ncy, ncx))# * 0.0

        # populate the different data structures with the new info
        # =========================================================
        if ref is True:
            print("NEW REFERENCE!")
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()
            fits.writeto("SHrefx.fits", self.SH_xref, clobber=True)
            fits.writeto("SHrefy.fits", self.SH_yref, clobber=True)

        self.SH_xslp = self.SH_xtmp - self.SH_xref
        self.SH_yslp = (self.SH_ytmp - self.SH_yref)

        #self.SH_phot -= self.threshold
        discard_them = self.SH_phot <= 1500#self.threshold * 2
        self.SH_phot[discard_them] = 0.0
        self.SH_xslp[discard_them] = 0.0
        self.SH_yslp[discard_them] = 0.0

        self.SH_comb[0] = self.SH_xslp
        self.SH_comb[1] = self.SH_yslp
        self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

        self.shm_comb.set_data(self.SH_comb)
        self.shm_phot_inst.set_data(self.SH_phot)

        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
        self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)

#        print(self.ttx_mean, self.tty_mean)
        self.update_log()

    def set_thread_params(self):
        print(self.SH_xx)
        print(self.SH_yy)

        print('self.bThreadInitialized: %r' % (self.bThreadInitialized))

        if (self.bUseMT == True):
            # if (self.bSetThreadParams == False):
            ncx, ncy = self.SH_xx.size - 1, self.SH_yy.size - 1
            self.nCells = ncx * ncy
            nsub = self.nCells // self.nThreads
            cellCounter = 0
            print('nCells, nThreads, nsub: %d %d %d' %(self.nCells,self.nThreads,nsub))
            self.vInd_i = np.zeros((self.nThreads-1,nsub),dtype=int)
            self.vInd_j = np.zeros((self.nThreads-1,nsub),dtype=int)
            for k in range(self.nThreads - 1):
                for u in range(nsub):
                    self.vInd_i[k,u] = cellCounter % ncx
                    self.vInd_j[k,u] = cellCounter // ncx
                    cellCounter += 1
                #     print('threadID, cellCounter, i, j: %d %d %d %d' % (k, cellCounter, self.vInd_i[k,u], self.vInd_j[k,u]))
                # print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_i[k,:])))
                # print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_j[k,:])))
                retval = self.lib.set_thread_params(k,nsub,self.SH_xx,self.SH_yy,self.vInd_i[k,:],self.vInd_j[k,:],self.SH_xtmp,self.SH_ytmp,self.SH_phot)

            rem_sub = self.nCells - nsub * (self.nThreads - 1)
            self.vInd_i2 = np.zeros(rem_sub,dtype=int)
            self.vInd_j2 = np.zeros(rem_sub,dtype=int)
            for u in range(rem_sub):
                self.vInd_i2[u] = cellCounter % ncx
                self.vInd_j2[u] = cellCounter // ncx
                cellCounter += 1
            #     print('threadID, cellCounter, i, j: %d %d %d %d' % (self.nThreads-1, cellCounter, self.vInd_i2[u], self.vInd_j2[u]))
            # print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_i2)))
            # print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_j2)))
            retval = self.lib.set_thread_params(self.nThreads-1,rem_sub,self.SH_xx,self.SH_yy,self.vInd_i2,self.vInd_j2,self.SH_xtmp,self.SH_ytmp,self.SH_phot)
            # self.bThreadInitialized = False
                # print('Thread initialized!')
            # else:
            #     print('Threads already initialized.')
            #     for k in range(self.nThreads - 1):
            #         print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_i[k,:])))
            #         print('vInd_j Pyt: [ %s]' % ' '.join(map(str, self.vInd_j[k,:])))
            #     print('vInd_i Pyt: [ %s]' % ' '.join(map(str, self.vInd_i2)))
            #     print('vInd_j Pyt: [ %s]' % ' '.join(map(str, self.vInd_j2)))



    # =========================================================
    def calc_SH_data(self, data=None, ref=True):
        ''' -------------------------------------------------------
        Calculates the Shack-Hartman data for the current image.

        Parameters:
        - ref: if True -> current image taken as new reference (bool)
        - data: if not None, info computed from the provided image

        This version of the function was idenfied as limiting the
        frequency of the loop. It is for now kept as a reference.
        ------------------------------------------------------- '''
        # clockA = time.perf_counter()
        xx, yy   = self.SH_xx , self.SH_yy
        ncx, ncy = xx.size - 1, yy.size - 1

        if data is None:
            self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True,sleepT=0)
        else:
            self.live_img = data

        # clockB = time.perf_counter()
        # print(clockB - clockA)

        #~ fitsio.write("./xx.fits",xx)
        #~ fitsio.write("./yy.fits",yy)
        #~ fitsio.write("./live_img.fits",self.live_img)

        self.im_cnt = self.shm_im.get_counter()  # image counter
        diff_cnt = self.im_cnt - self.prev_im_cnt
        if diff_cnt <= 20:
            self.vFrameDrop[diff_cnt] = self.vFrameDrop[diff_cnt] + 1
        else:
            self.vFrameDrop[0] = self.vFrameDrop[0] + 1
        self.prev_im_cnt = self.im_cnt

        bckgd = self.live_img.mean()
        # print('bkg >> %f' %(bckgd,))


        vImg = self.live_img.astype(int)

        # clockA = time.perf_counter()
        if (self.bUseMT == True):
            # retval = self.lib.trigger_calc_SH_data_mt(self.sizex,self.sizey,vImg,ncx,ncy,xx,yy,self.SH_xtmp,self.SH_ytmp,self.SH_phot)
            retval = self.lib.trigger_calc_SH_data_mt(vImg,bckgd)
        else:
            retval = self.lib.calc_SH_data(self.sizex,self.sizey,vImg,bckgd,ncx,ncy,xx,yy,self.SH_xtmp,self.SH_ytmp,self.SH_phot)
        # clockB = time.perf_counter()
        # print(clockB - clockA)
        # print("calc_SH_data timing: %.3e %.3e %d" % (time.time() - time_start, time_diff/mycount,mycount))

        # print("calc_SH_data timing: %.3e %.3e %d" % (time.time() - time_start, time_diff/mycount,mycount))

        if ref == True:
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()
            # fits.writeto("C:\\Home\\Downloads\\SH_xref.fits",self.SH_xref,overwrite=True)
            # fits.writeto("C:\\Home\\Downloads\\SH_yref.fits",self.SH_yref,overwrite=True)
            # fits.writeto("C:\\Home\\Downloads\\SH_phot.fits", self.SH_phot, overwrite=True)
            fits.writeto("C:\\Home\\Downloads\\vImg.fits",vImg,overwrite=True)

        self.SH_xslp = self.SH_xtmp - self.SH_xref
        self.SH_yslp = self.SH_ytmp - self.SH_yref

        self.SH_phot -= self.threshold
        self.SH_xslp[self.SH_phot < 0] = 0.0
        self.SH_yslp[self.SH_phot < 0] = 0.0

        self.SH_comb[0] = self.SH_xslp
        self.SH_comb[1] = self.SH_yslp
        self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

        self.shm_comb.set_data(self.SH_comb)
        self.shm_phot_inst.set_data(self.SH_phot)

        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
        self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)

        self.ttx_tavg += self.ttx_mean
        self.tty_tavg += self.tty_mean
        self.tavg_count += 1.

        self.update_log()


    # =========================================================
    def update_log(self,):
#        print(self.ttx_mean,self.tty_mean)
        self.ttx_log.append(self.ttx_mean)
        self.tty_log.append(self.tty_mean)
#        print(self.ttx_log[len(vLog)-1],self.tty_log[len(vLog)-1])
        if len(self.ttx_log) > self.log_len:
            self.ttx_log.pop(0)
            self.tty_log.pop(0)

    # =========================================================
    def loop(self,):
        self.keepgoing = True
        # counter = 0
        # acctime = 0
        # start_time = time.perf_counter()

        # self.bLogRTData = True
        # self.log_length = 10
        # self.LogDataSlopes = np.zeros([self.log_length,3,10,10],dtype='float32')
        # self.LogTimeStamp = np.zeros(self.log_length,dtype='float64')
        # self.LogDataCommands = np.zeros([self.log_length,self.nmodes])
        # self.LogDataErrors = np.zeros([self.log_length,self.nmodes])

        log_counter = 0

        clock_zero = time.perf_counter()

        while self.keepgoing:
            # clockA = time.perf_counter()
            self.calc_SH_data(ref=False)

            # if self.bLogRTData == True:
            #     # print('signal dims >> ',np.shape(signal))
            #     self.LogDataSlopes[log_counter] = self.SH_comb
            #     self.LogTimeStamp[log_counter] = clockA - clock_zero
            #     log_counter += 1
            #     if log_counter >= self.log_length:
            #         now = datetime.now()
            #         dt_string = now.strftime("%Y.%m.%d_%Hh%Mm%Ss")
            #         fits.writeto('C:\\Data\\' + dt_string + '_SlopeX.fits',self.LogDataSlopes[:,0,:,:].reshape([self.log_length,100]),overwrite=True)
            #         fits.writeto('C:\\Data\\' + dt_string + '_SlopeY.fits',self.LogDataSlopes[:,1,:,:].reshape([self.log_length,100]),overwrite=True)
            #         fits.writeto('C:\\Data\\' + dt_string + '_Illum.fits',self.LogDataSlopes[:,2,:,:].reshape([self.log_length,100]),overwrite=True)
            #         fits.writeto('C:\\Data\\' + dt_string + '_timestamps.fits',self.LogTimeStamp[:],overwrite=True)
            #         break
            # clockB = time.perf_counter()
            # dtime = (time.perf_counter() - clockA)
            # acctime += dtime
            # print(dtime)
            # counter += 1
            # if counter % 100 == 0:
            #     # freq = counter / (time.perf_counter() - start_time)
            #     # print('Loop frequency >> %.3f' % (freq))
            #     print('calc_SH_data avergae timing >> %.3e' % (acctime / counter))
            #     counter = 0
            #     # start_time = time.perf_counter()
        print("WFS monitoring is now off")

    # =========================================================
    def start(self,):
        print("start %s" % ("YEAH!", ))
        print(self.SH_xx)
        print(self.SH_yy)
        self.calc_SH_data(ref=True)

    # =========================================================
    def stop(self,):
        # release the CPUs used by multi-processing
        if self.mproc is True:
            print("OK")
            print("relieving the %d CPUs from duty" % (self.ncpu))
            self.pool.close()
# ==========================================================
# ==========================================================
if __name__ == "__main__":
    mon = WFS(shmf="R:\\ixon.im.shm")
    mon.start()
