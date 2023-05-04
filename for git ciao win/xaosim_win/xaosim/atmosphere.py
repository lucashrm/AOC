#!/usr/bin/env python3

''' ====================================================================
This is the atmospheric simulation module of XAOSIM.

It defines a generic atmospheric phase screen class that is to be used
in dynamic simulations.
==================================================================== '''

import numpy as np
import scipy as sci
import random as rnd
import threading
import time
from . import wavefront as wft
from .shmlib import shm
from scipy.ndimage.interpolation import shift
from scipy.ndimage.fourier import fourier_shift
from scipy import special as spe
import astropy.io.fits as fits

idl_sharp = lambda x, y: np.dot(x.T, y.T).T

# ===========================================================
# ===========================================================
class Phscreen(object):
    '''Atmospheric Kolmogorov-type phase screen.

    ====================================================================

    Class Attributes:
    ----------------
    - csz     : size (csz x csz) of the phase screen       (in pixels)
    - pdiam   : diameter of the aperture within this array (in pixels)
    - rndarr  : uniformly distributed random array         (csz x csz)
    - kolm    : the original phase screen                  (csz x csz)
    - kolm2   : the oversized phase screen    ((csz + pdiam) x (csz + pdiam))
    - qstatic : an optional quasi static aberration    (pdiam x pdiam)
    - rms     : total phase screen rms value           (in nanometers)
    - rms_i   : instant rms inside the pupil           (in nanometers)

    Comment:
    -------
    While the attributes are documented here for reference, the prefered
    way of interacting with them is via the functions defined within the
    class.
    ====================================================================

    '''
    # ==================================================
    def __init__(self, name="MaunaKea", csz = 512,
                 lsz=8.0, r0=0.2, L0=10.0,
                 fc=24.5, correc=1.0,
                 #shmf='phscreen.wf.shm', shdir='/dev/shm/'):
                 shmf='phscreen.wf.shm', shdir='R:\\'):

        ''' Kolmogorov type atmosphere + qstatic error

        -----------------------------------------------------
        Parameters:
        ----------
        - name  : a string describing the instrument
        - csz   : the size of the Fourier array
        - lsz   : the screen linear size (in meters)
        - r0    : the Fried parameter (in meters)
        - L0    : the outer scale parameter (in meters)
        - shmf  : file name to point to shared memory
        - shdir : location of the shm "files"
        -----------------------------------------------------
        '''
        self.shmf    = shmf
        self.shdir   = shdir
        self.csz     = csz
        self.lsz     = lsz
        self.r0      = r0
        self.L0      = L0

        self.rms_i   = 0.0
        self.correc  = correc
        self.fc      = fc
        self.kolm    = wft.atmo_screen(csz, lsz, r0, L0, fc, correc).real

        self.qstatic = np.zeros((self.csz, self.csz))
        self.shm_phs = shm(shdir + shmf,
                           data = self.qstatic, verbose=False)

        self.kolm2   = np.tile(self.kolm, (2,2))
        #self.kolm2   = self.kolm2[:self.sz+self.pdiam,:self.sz+self.pdiam]

        self.keepgoing = False

        self.offx = 0 # x-offset on the "large" phase screen array
        self.offy = 0 # y-offset on the "large" phase screen array

        self.ttc     = False   # Tip-tilt correction flag

        # auxilliary array (for tip-tilt correction)
        self.xx, self.yy  = np.meshgrid(np.arange(self.csz)-self.csz//2,
                                        np.arange(self.csz)-self.csz//2)
        self.xxnorm2 = np.sum(self.xx**2)
        self.yynorm2 = np.sum(self.yy**2)

        # self.bInitialize = True
        self.ScreenSize_pix = self.csz
        self.ScreenSize_m = 1.25
        self.r0_m = self.r0
        self.L0_m = self.L0
        self.LayerSpeed_mps = 4.
        self.TimeStep_s = 0.002

        self.CumulProgressPix = [0.]
        self.matrixA = np.empty([self.ScreenSize_pix,self.ScreenSize_pix*2],dtype='float64')
        self.matrixB = np.empty([self.ScreenSize_pix,self.ScreenSize_pix],dtype='float64')
        self.seed = None
        self.PhaseScreen_Buffer = np.empty([self.ScreenSize_pix+1,self.ScreenSize_pix],dtype='float64')
        self.PhaseScreen_rad = np.empty([self.ScreenSize_pix,self.ScreenSize_pix],dtype='float64')

        print('self.PhaseScreen_Buffer: ' , hex(id(self.PhaseScreen_Buffer)))
        print('self.PhaseScreen_rad: ' , hex(id(self.PhaseScreen_rad)))

        # INITIALIZE PHASE SCREEN ROUTINE
        self.bUseProgresivePhaseScreen = False
        if self.bUseProgresivePhaseScreen == True:
            self.progressive_phase_screen(True, self.L0_m, self.r0_m, self.LayerSpeed_mps, self.TimeStep_s)

        print('HereD:')
        print('self.PhaseScreen_Buffer: ' , hex(id(self.PhaseScreen_Buffer)))
        print('self.PhaseScreen_rad: ' , hex(id(self.PhaseScreen_rad)))
        print(np.shape(self.PhaseScreen_Buffer))
        print(np.shape(self.PhaseScreen_rad))


    # ==============================================================
    def start(self, delay=0.1, screen_speed=0.1, factor=1., wf_mode='atmo_fft'):
        ''' ----------------------------------------
        High-level accessor to start the thread of
        the phase screen server infinite loop
        ---------------------------------------- '''
        if not self.keepgoing:

            self.kolm2   = np.tile(self.kolm, (2,2))
            #self.kolm2   = self.kolm2[:self.sz+self.pdiam,:self.sz+self.pdiam]

            self.keepgoing = True
            t = threading.Thread(target=self.__loop__, args=(delay, screen_speed, factor, wf_mode))
            t.start()
            print("The *ATMO* phase screen server was started")
        else:
            print("The *ATMO* phase screen server is already running")

    # ==============================================================
    def freeze(self):
        ''' ----------------------------------------
        High-level accessor to interrupt the thread
        of the phase screen server infinite loop
        ---------------------------------------- '''
        if self.keepgoing:
            self.keepgoing = False
        else:
            print("The *ATMO* server was frozen")


    # ==============================================================
    def stop(self):
        ''' ----------------------------------------
        High-level accessor to interrupt the thread
        of the phase screen server infinite loop
        ---------------------------------------- '''
        if self.keepgoing:
            self.kolm2[:] = 0.0
            time.sleep(0.5)
            self.keepgoing = False
            print("The *ATMO* server was stopped")
        else:
            print("The *ATMO* server was not running")


    # ==============================================================
    def set_qstatic(self, qstatic=None):
        if qstatic is not None:
            if qstatic.shape == (self.csz, self.csz):
                self.qstatic = qstatic
                print("quasi-static phase screen updated!")
            else:
                print("could not update quasi-static phase screen")
                print("array should be %d x %d (dtype=%s)" % (
                    self.csz, self.csz, str(self.qstatic.dtype)))

            # if simulation is not already active, update phase screen!
            if self.keepgoing == False:
                self.shm_phs.set_data(self.qstatic)

        else:
            print("no new quasi-static screen was provided!")

    # ==============================================================
    def update_screen(self, correc=None, fc=None, r0=None, L0=None):
        ''' ------------------------------------------------
        Generic update of the properties of the phase-screen

        ------------------------------------------------ '''
        if r0 is not None:
            self.r0 = r0

        if L0 is not None:
            self.L0 = L0

        if correc is not None:
            self.correc = correc

        if fc is not None:
            self.fc = fc

        self.kolm    = wft.atmo_screen(
            self.csz, self.lsz, self.r0, self.L0, self.fc, self.correc).real

        self.kolm2   = np.tile(self.kolm, (2,2))

        if self.keepgoing is False:
            # case that must be adressed:
            # amplitude changed when atmo is frozen!
            if False:
                subk = self.kolm2[self.offx:self.offx+self.csz,
                                  self.offy:self.offy+self.csz].copy()

            if True:
                subk = self.xx * np.cos(self.time * 2. * np.pi / (TimeStep_s * 10.)) + self.yy * np.sin(self.time * 2. * np.pi / (TimeStep_s * 10.))

            if self.ttc is True:
                ttx = np.sum(subk*self.xx) / self.xxnorm2
                tty = np.sum(subk*self.yy) / self.yynorm2
                subk -= ttx * self.xx + tty * self.yy

            self.rms_i = subk.std()
            self.shm_phs.set_data(subk + self.qstatic)

    # ==============================================================
    def progressive_phase_screen(self, bInitialize=True, L0_m=10., r0_m=0.1, LayerSpeed_mps=4., TimeStep_s=0.001):

        # print('self.PhaseScreen_Buffer: ' , hex(id(self.PhaseScreen_Buffer)))
        # print('self.PhaseScreen_rad: ' , hex(id(self.PhaseScreen_rad)))
        # print('self.CumulProgressPix/CumulProgressPix: ' , hex(id(self.CumulProgressPix[0])), hex(id(self.CumulProgressPix[0])))

        size = np.int(self.ScreenSize_pix)
        Zsize = 2
        Ps = np.double(self.ScreenSize_m / np.double(size))
        Z = np.zeros(size*Zsize)

        # FOR INITIALIZATION, ASSUME r0_m/L0_m = 1.D, AND SCALE B LATER BY (L0_m/r0_m)**(5/6)
        if bInitialize == False:
            self.PhaseScreen_Copy = np.copy(self.PhaseScreen_Buffer[0:size,:])

            # fits.writeto(r'C:\Home\Downloads\PhaseScreen_Copy.fits',PhaseScreen_Copy,overwrite=True)

            LayerSpeed_PixPerSample = np.abs(size * self.LayerSpeed_mps / self.ScreenSize_m * self.TimeStep_s)

            ResProgressPix = LayerSpeed_PixPerSample - np.int(LayerSpeed_PixPerSample)

            DiffProgress = self.CumulProgressPix[0] + ResProgressPix
            AdditionalCol = np.int(LayerSpeed_PixPerSample) + np.int(DiffProgress)

            # ;pm, 'CumulProgressPix, LayerProgressPix, DiffProgress, AdditionalCol >> ',CumulProgressPix,ResProgressPix,DiffProgress,AdditionalCol

            DiffProgress = DiffProgress - np.int(DiffProgress)

            # print('AdditionalCol : ',AdditionalCol)

            if AdditionalCol >= 1:
                # ;pm, 'Shifting PhaseScreen_Buffer by one column and fill PhaseScreen_Copy.'
                # ;; AN ADDITIONAL COLUMN HAS ALREADY BEEN ADDED AT THE PREVIOUS ITERATION
                self.PhaseScreen_Copy[:,:] = self.PhaseScreen_Buffer[1:size+1,:]

                for u in range(1,AdditionalCol):
                    # print('Adding column(s): ',AdditionalCol)
                    for i in range(0,size*Zsize):
                        ii = i % Zsize
                        jj = i // Zsize
                        Z[i] = self.PhaseScreen_Copy[size-Zsize + ii,jj]

                    beta_vect = np.random.normal(loc=0.0,scale=1.,size=size)
                    beta_vect = beta_vect / np.sqrt(np.std(beta_vect))

                    # tmpA = matrixA#Z
                    # tmpB = (matrixB * (L0_m/r0_m)**(5.D/6.D))#beta_vect
                    tmpA = np.dot(self.matrixA, Z)
                    tmpB = np.dot(self.matrixB.T * (L0_m/r0_m)**(5./6.), beta_vect)
                    # tmpA = idl_sharp(matrixA, Z)
                    # tmpB = idl_sharp(matrixB * (L0_m/r0_m)**(5./6.), beta_vect)

                    X = tmpA + tmpB

                    # NewPhase = SHIFT(PhaseScreen_Copy,-1L,0)
                    # NewPhase[size-1L,*] = X
                    NewPhase = np.roll(self.PhaseScreen_Copy,-1,axis=0)
                    # NewPhase[size-1,*] = X
                    # NewPhase = fourier_shift(np.fft.fftn(PhaseScreen_Copy), np.array([-1, 0. ]))
                    # NewPhase[size-1,0:] = X
                    NewPhase[size,:] = X

                    self.PhaseScreen_Copy[:,:] = NewPhase

                    # ;pm, 'Added one more column to PhaseScreen_Copy ',u,'/',AdditionalCol-1L

                # ;; ADD A COLUMN TO THE SCREEN BUFFER ARRAY PhaseScreen_Buffer
                for i in range(0,size*Zsize):
                    ii = i % Zsize
                    jj = i // Zsize
                    Z[i] = self.PhaseScreen_Copy[size-Zsize + ii,jj]

                beta_vect = np.random.normal(loc=0.0,scale=1.,size=size)
                beta_vect = beta_vect / np.sqrt(np.std(beta_vect))

                # print('matrixA: ',np.shape(self.matrixA),np.shape(self.matrixB))
                # tmpA = matrixA#Z
                # tmpB = (matrixB * (L0_m/r0_m)**(5.D/6.D))#beta_vect
                tmpA = np.dot(self.matrixA, Z)
                tmpB = np.dot(self.matrixB.T * (self.L0_m/self.r0_m)**(5./6.), beta_vect)
                # tmpA = idl_sharp(matrixA, Z)
                # tmpB = idl_sharp(matrixB * (L0_m/r0_m)**(5./6.), beta_vect)

                X = tmpA + tmpB

                self.PhaseScreen_Buffer[0:size,:] = self.PhaseScreen_Copy[:,:]
                self.PhaseScreen_Buffer[size,:] = X

            # fits.writeto(r'C:\Home\Downloads\PhaseScreen_Buffer.fits',self.PhaseScreen_Buffer,overwrite=True)
            # print('DiffProgress: ',DiffProgress)
            # DiffProgress = -1.
            # ShiftedPhase = np.roll(PhaseScreen_Buffer,DiffProgress,axis=1)
            FourierShiftedPhase = fourier_shift(np.fft.fftn(self.PhaseScreen_Buffer), np.array([-DiffProgress, 0. ]))
            ShiftedPhase = np.fft.ifftn(FourierShiftedPhase).real
            # shifted_image_cubic = cv2.warpAffine(image, Mtrans, image.shape[::-1], flags=cv2.INTER_CUBIC)

            # fits.writeto(r'C:\Home\Downloads\ShiftedPhase.fits',ShiftedPhase,overwrite=True)
            self.PhaseScreen_rad[:,:] = ShiftedPhase[0:size,:]

            # print('1. CumulProgressPix: ',self.CumulProgressPix[0], hex(id(self.CumulProgressPix[0])))
            self.CumulProgressPix[0] += np.abs(ResProgressPix)
            # print('2. CumulProgressPix: ',self.CumulProgressPix[0], hex(id(self.CumulProgressPix[0])))

            self.CumulProgressPix[0] = self.CumulProgressPix[0] - np.int(self.CumulProgressPix[0])
            # print('3. CumulProgressPix: ',self.CumulProgressPix[0], hex(id(self.CumulProgressPix[0])))
        else:
            # print('Initializing progressive phase screen (Assemat)')
            self.CumulProgressPix[0] = 0.

            # PhaseScreen_rad = VK_aberration_screen(ScreenSize_pix,ScreenSize_m,r0_m,L0_m,seed,0,'')
            # np.copyto(PhaseScreen_rad,self.kolm.real)
            self.PhaseScreen_rad[:,:] = fits.getdata(r'c:\Home\Downloads\testVK.fits')
            # print('HereE:')
            # print('PhaseScreen_rad: ' , hex(id(self.PhaseScreen_rad)))
            # print(np.shape(self.PhaseScreen_rad))
            # fits.writeto(r'c:\home\downloads\test.fits',PhaseScreen_rad,overwrite=True)
            self.PhaseScreen_Copy = np.copy(self.PhaseScreen_rad)

            ZZt = np.zeros([Zsize*size,Zsize*size])
            XZt = np.zeros([size,Zsize*size])
            XXt = np.zeros([size,size])

            counti = 0
            countj = 0

            for i in range(0,size*Zsize):
                countj = 0
                yi = counti // Zsize
                xi = counti - yi*Zsize

                for j in range(0,size*Zsize):
                    if i == j:
                        countj += 1
                        continue
                    yj = countj // Zsize
                    xj = countj - yj*Zsize

                    r = Ps * np.sqrt(np.double(xj - xi)**2 + np.double(yj - yi)**2)

                    # ;ZZt[i,j] = (L0_m/r0_m)**(5.D/3.D) * GAMMA(11.D/6.D) / (2.D**(5.D/6.D) * np.pi**(8.D/3.D)) * $
                    ZZt[i,j] = spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * \
                        (24./5. * spe.gamma(6./5.))**(5./6.) * (2.*np.pi*r / L0_m)**(5./6.) * \
                        spe.kv( 5./6., 2. * np.pi * r / L0_m) / (np.pi*0.5)

                    countj += 1

                counti += 1

            # ;for i=0L,size*Zsize-1L do ZZt[i,i] = (L0_m/r0_m)**(5.D/3.D) * GAMMA(11.D/6.D) / (2.D**(5.D/6.D) * np.pi**(8.D/3.D)) * $
            val = spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * \
                    (24./5. * spe.gamma(6./5.))**(5./6.) * (2. * np.pi * 0.0001 / L0_m)**(5/6) * \
                    spe.kv( 5./6., 2. * np.pi * 0.0001 / L0_m) / (np.pi*0.5)
            for i in range(0,size*Zsize):
                ZZt[i,i] = val

            counti = 0
            countj = 0
            for i in range(0,size):
                countj = 0
                xi = size
                yi = counti
                for j in range(0,size*Zsize):
                    yj = countj // Zsize
                    xj = size-Zsize + countj - yj*Zsize

                    r = Ps * np.sqrt(np.double(xj - xi)**2 + np.double(yj - yi)**2)

                    # ;XZt[i,j] = (L0_m/r0_m)**(5./3.) * spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) *
                    XZt[i,j] = spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * \
                    	(24./5. * spe.gamma(6./5.))**(5./6.) * (2.*np.pi*r / L0_m)**(5./6.) * \
                    	spe.kv( 5./6., 2. * np.pi * r / L0_m) / (np.pi*0.5)

                    countj += 1

                counti += 1

            counti = 0
            countj = 0
            for i in range(0,size):
                countj = 0
                yi = counti // size
                xi = counti - yi*size
                for j in range(0,size):
                    if i == j:
                        countj += 1
                        continue

                    yj = countj // size
                    xj = countj - yj*size

                    r = Ps * np.sqrt(np.double(xj - xi)**2 + np.double(yj - yi)**2)

                    # ;XXt[i,j] = (L0_m/r0_m)**(5./3.) * spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) *
                    XXt[i,j] = 1. * spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * \
                    	(24./5. * spe.gamma(6./5.))**(5./6.) * (2.*np.pi*r / L0_m)**(5./6.) * \
                    	spe.kv( 5./6., 2. * np.pi * r / L0_m) / (np.pi*0.5)

                    countj += 1

                counti += 1

            # ;for i=0L,size-1L do XXt[i,i] = (L0_m/r0_m)**(5./3.) * GAMMA(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * $
            for i in range(0,size):
                XXt[i,i] = spe.gamma(11./6.) / (2.**(5./6.) * np.pi**(8./3.)) * \
                    (24./5. * spe.gamma(6./5.))**(5./6.) * (2. * np.pi * 0.0001 / L0_m)**(5/6) * \
                    spe.kv( 5./6., 2. * np.pi * 0.0001 / L0_m) / (np.pi*0.5)

            # print('XZt dims: ',np.shape(XZt))
            # print('ZZt dims: ',np.shape(ZZt))
            # matrixA = XZt # INVERT(ZZt)
            # BBt = XXt - matrixA # ZZt # TRANSPOSE(matrixA)
            self.matrixA = np.dot(XZt,np.linalg.pinv(ZZt))
            BBt = XXt - np.dot(self.matrixA, np.dot(ZZt, self.matrixA.T))
            # matrixA = idl_sharp(XZt, np.linalg.pinv(ZZt))
            # BBt = XXt - idl_sharp(matrixA, idl_sharp(ZZt, matrixA.T))

            # fits.writeto(r'C:\Home\Downloads\BBt.fits',BBt,overwrite=True)

            # SVDC, BBt, W, U, V, /DOUBLE
            U, W, V = np.linalg.svd(BBt)

            # fits.writeto(r'C:\Home\Downloads\U.fits',U,overwrite=True)
            # fits.writeto(r'C:\Home\Downloads\V.fits',V,overwrite=True)
            # fits.writeto(r'C:\Home\Downloads\W.fits',W,overwrite=True)

            SV = np.zeros([size, size])
            for k in range(0,size):
                SV[k,k] = W[k]

            # ;BBt_check = U ## SV ## TRANSPOSE(V)
            BBt_check = U @ SV @ V
            # fits.writeto(r'C:\Home\Downloads\BBt_check.fits',BBt_check,overwrite=True)

            # matrixB = TRANSPOSE(U ## SQRT(SV))
            # matrixB = (U @ np.sqrt(SV)).T
            self.matrixB = np.dot(U, np.sqrt(SV)).T

            # ;BBt_check2 = B # TRANSPOSE(B)

            # ;pm, 'BBt - BBt_check  >> ',total(BBt - BBt_check,/DOUBLE)
            # ;pm, 'BBt - BBt_check2 >> ',total(BBt - BBt_check2,/DOUBLE)

            for k in range(0,size):
                for i in range(0,size*Zsize):
                    ii = i % Zsize
                    jj = i // Zsize
                    Z[i] = self.PhaseScreen_rad[size-Zsize + ii,jj]

                np.random.seed()
                # beta_vect = np.zeros(size)
                beta_vect = np.random.normal(loc=0.0,scale=1.0,size=size)
                # for ii in range(0,size):
                #     beta_vect[ii] = rnd.gauss(0.,1.)

                beta_vect = beta_vect / np.sqrt(np.std(beta_vect))

                # ;r0 = r0_min + (r0_max - r0_min)*0.5D*(1. + cos(k*2.*np.pi / 1000.))
                # ;pm, 'r0 >> ',k,r0,r0 / (DiamTel_m)

                # tmpA = matrixA#Z
                # tmpB = (matrixB * (L0_m/r0_m)**(5./6.))#beta_vect

                tmpA = np.dot(self.matrixA, Z)
                tmpB = np.dot(self.matrixB.T * (L0_m/r0_m)**(5./6.), beta_vect)
                # tmpA = idl_sharp(matrixA, Z)
                # tmpB = idl_sharp(matrixB * (L0_m/r0_m)**(5./6.),beta_vect)

                X = tmpA + tmpB

                # NewPhase = SHIFT(PhaseScreen_rad,-1,0)
                # NewPhase[size-1L,*] = X
                NewPhase = np.roll(self.PhaseScreen_rad,-1,axis=0)
                # NewPhase[size-1,0:] = X
                # NewPhase = fourier_shift(np.fft.fftn(PhaseScreen_rad), np.array([-1, 0. ]))
                NewPhase[size-1,:] = X

                # ;wset,10
                # ;tvscl, congrid(NewPhase,winsx,winsy)

                self.PhaseScreen_rad[:,:] = NewPhase
                # fits.writeto(r'C:\Home\Downloads\PhaseScreen_rad%03d.fits' % (k),self.PhaseScreen_rad,overwrite=True)

                # ;PhaseScreen_Copy = NewPhase

            # PhaseScreen_Buffer = np.zeros([size+1,size])
            # print('A-PhaseScreen_Buffer/PhaseScreen_rad size: ',np.shape(PhaseScreen_Buffer),np.shape(PhaseScreen_rad))
            # PhaseScreen_Buffer.resize(size+1,size)
            # print('B-PhaseScreen_Buffer/PhaseScreen_rad size: ',np.shape(PhaseScreen_Buffer),np.shape(PhaseScreen_rad))
            self.PhaseScreen_Buffer[0:size,:] = self.PhaseScreen_rad

            for i in range(0,size*Zsize):
                ii = i % Zsize
                jj = i // Zsize
                Z[i] = self.PhaseScreen_rad[size-Zsize + ii,jj]

            # beta_vect = randomn(seed,size,/DOUBLE)
            # beta_vect = beta_vect / sqrt((MOMENT(beta_vect))[1])
            beta_vect = np.random.normal(loc=0.0,scale=1.0,size=size)
            beta_vect = beta_vect / np.sqrt(np.std(beta_vect))

            # ;r0 = r0_min + (r0_max - r0_min)*0.5D*(1. + cos(k*2.*np.pi / 1000.))
            # ;pm, 'r0 >> ',k,r0,r0 / (DiamTel_m)

            # tmpA = matrixA#Z
            # tmpB = (matrixB * (L0_m/r0_m)**(5./6.))#beta_vect
            tmpA = np.dot(self.matrixA,Z)
            tmpB = np.dot(self.matrixB * (L0_m/r0_m)**(5./6.), beta_vect)
            # tmpA = idl_sharp(matrixA,Z)
            # tmpB = idl_sharp(matrixB * (L0_m/r0_m)**(5./6.), beta_vect)

            X = tmpA + tmpB

            # print('HereA:')
            # print(hex(id(self.PhaseScreen_Buffer)))
            # print(np.shape(self.PhaseScreen_Buffer))
            self.PhaseScreen_Buffer[size,:] = X
            # print(np.shape(self.PhaseScreen_Buffer))

            self.CumulProgressPix[0] = 0.

            bInitialize = False

            # fits.writeto(r'C:\Home\Downloads\matrixA.fits',self.matrixA,overwrite=True)
            # fits.writeto(r'C:\Home\Downloads\matrixB.fits',self.matrixB,overwrite=True)

    # ==============================================================
    def __loop__(self, delay=0.1, screen_speed=0.1, factor=1., wf_mode='atmo_fft'):
        ''' ------------------------------------------
        Main loop: frozen screen slid over the aperture

        Options:
        ---------
        - delay: the time delay between refresh  (0.1 sec)
        -----------------------------------------  '''
        print('Entering loop! (cam: %.3e  screen: %.3e  factor: %.3e  wf_mode: %s)' % (delay, screen_speed, factor, wf_mode))

        self.offx += 2
        self.offy += 1
        self.offx = self.offx % self.csz
        self.offy = self.offy % self.csz

        subk = self.kolm2[self.offx:self.offx+self.csz,
                          self.offy:self.offy+self.csz].copy()

        prev_time_cam = time.perf_counter()
        prev_time_screen = prev_time_cam
        self.screen_speed = screen_speed
        self.cam_speed = delay

        print('Starting wafevront simulation with wf_mode=',wf_mode)

        while self.keepgoing:
            # time.sleep(delay)
            mytime = time.perf_counter()
            self.time = mytime

            diff_time_cam = mytime - prev_time_cam
            diff_time_screen = mytime - prev_time_screen

            if wf_mode == 'atmo_fft':
                self.offx += 2
                self.offy += 1
                self.offx = self.offx % self.csz
                self.offy = self.offy % self.csz

                subk = self.kolm2[self.offx:self.offx+self.csz,
                                  self.offy:self.offy+self.csz].copy()
            elif wf_mode == 'atmo_progress':
                self.progressive_phase_screen(False,self.L0_m,self.r0_m,self.LayerSpeed_mps,self.TimeStep_s)

                self.rms_i = np.std(self.PhaseScreen_rad)

                self.shm_phs.set_data(self.PhaseScreen_rad * factor + self.qstatic)

                prev_time_screen = mytime
            elif wf_mode == 'const_tilt':
                cx = factor
                cy = 0.
                # print('cx, cy >> ', cx, cy)
                subk = (self.xx * cx + \
                    self.yy * cy ) / np.max(self.xx) * np.pi * 2.
            elif wf_mode == 'oscillating':
                cx = np.cos(self.time * 2. * np.pi / (self.TimeStep_s * 50.)) * factor
                cy = np.sin(self.time * 2. * np.pi / (self.TimeStep_s * 50.)) * 0.
                subk = (self.xx * cx + \
                    self.yy * cy ) / np.max(self.xx) * np.pi * 2.
            elif wf_mode == 'circle':
                cx = np.cos(self.time * 2. * np.pi / (self.TimeStep_s * 50.)) * factor
                cy = np.sin(self.time * 2. * np.pi / (self.TimeStep_s * 50.)) * factor
                subk = (self.xx * cx + \
                    self.yy * cy ) / np.max(self.xx) * np.pi * 2.

            if self.ttc is True:
                ttx = np.sum(subk*self.xx) / self.xxnorm2
                tty = np.sum(subk*self.yy) / self.yynorm2
                subk -= ttx * self.xx + tty * self.yy

            self.rms_i = subk.std()
            self.shm_phs.set_data(subk * factor + self.qstatic)
            prev_time_screen = mytime

            time.sleep(screen_speed)
