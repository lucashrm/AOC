#!/usr/bin/env python

import numpy as np
from xaosim.shmlib import shm
from xaosim.zernike import mkzer1
from numpy.linalg import solve
from numpy.linalg import pinv
import astropy.io.fits as fits
from datetime import datetime

import sys
import time

import os

import cv2

# =============================================================================
# =============================================================================
def get_coadd(shm1, shm2=None, nav=20, reform=True):
    ''' -----------------------------------------------------------------------
    A useful tool that produces a time average of the last *nav* maps provided
    by one or two *shm* shared memory structure
    ----------------------------------------------------------------------- '''
    try:
        im_cnt1 = shm1.get_counter()
        if shm2 is not None:
            im_cnt2 = shm2.get_counter()
    except:
        print("get_coadd error: arguments is not a valid shm structure?")
        return

    myim1   = shm1.get_data(check=im_cnt1, reform=reform)
    im_cnt1 = shm1.get_counter()
    if shm2 is not None:
        myim2   = shm2.get_data(check=im_cnt2, reform=reform)
        im_cnt2 = shm2.get_counter()

    for ii in range(nav-1):
        myim1  += shm1.get_data(check=im_cnt1, reform=reform)
        im_cnt1 = shm1.get_counter()
        if shm2 is not None:
            myim2  += shm2.get_data(check=im_cnt2, reform=reform)
            im_cnt2 = shm2.get_counter()

    myim1 /= nav
    if shm2 is not None:
        myim2 /= nav
        return((myim1, myim2))
    else:
        return(myim1)

# =============================================================================
# =============================================================================

class DM(object):
    ''' -----------------------------------------------------------------------
    High level information regarding the DM
    ----------------------------------------------------------------------- '''
    def __init__(self,):
        self.nba = 97  # number of actuators for the DM 97-15
        self.dms = 11  # equivalent grid size of the DM 97-15 (11x11 - corners)
        dms = self.dms
        self.dmmask = np.ones((dms, dms), dtype=np.int_) # corner mask

        # self.dmmask[np.abs(self.xdm) + np.abs(self.ydm) > 7] = 0.0

        self.dmmask[0,0:3] = 0
        self.dmmask[0,dms-3:] = 0
        self.dmmask[dms-1,0:3] = 0
        self.dmmask[dms-1,dms-3:] = 0

        self.dmmask[1,0:2] = 0
        self.dmmask[1,dms-2:] = 0
        self.dmmask[dms-2,dms-2:] = 0
        self.dmmask[dms-2,0:2] = 0

        self.dmmask[2,0] = 0
        self.dmmask[2,dms-1] = 0
        self.dmmask[dms-3,0] = 0
        self.dmmask[dms-3,dms-1] = 0

        # CENTRAL OBSCURATION
        # xc = dms // 2
        # yc = dms // 2
        # rad0sq = (dms / 2. * 0.35)**2
        # for j in range(0,dms):
        #     for i in range(0,dms):
        #         sqdist = (i - xc)**2 + (j - yc)**2
        #         if sqdist < rad0sq:
        #             self.dmmask[j,i] = 0


        # INDEX OF nact_lin_x*nact_lin_y ARRAY TO ACTUATOR ARRAY
        self.vActMapping = np.zeros(self.nba).astype(np.int32)
        offs = 3
        self.vActMapping[0:5] = np.array(offs+np.arange(0,5)).astype(np.int32)
        offs += 3 + 2 # +8
        self.vActMapping[5:12] = np.array(offs+np.arange(5,12)).astype(np.int32)
        offs += 2 + 1 # +11
        self.vActMapping[12:21] = np.array(offs+np.arange(12,21)).astype(np.int32)
        offs += 1 # +12
        self.vActMapping[21:32] = np.array(offs+np.arange(21,32)).astype(np.int32)
        self.vActMapping[32:43] = np.array(offs+np.arange(32,43)).astype(np.int32)
        self.vActMapping[43:54] = np.array(offs+np.arange(43,54)).astype(np.int32)
        self.vActMapping[54:65] = np.array(offs+np.arange(54,65)).astype(np.int32)
        self.vActMapping[65:76] = np.array(offs+np.arange(65,76)).astype(np.int32)
        offs += 1 # +13
        self.vActMapping[76:85] = np.array(offs+np.arange(76,85)).astype(np.int32)
        offs += 1 + 2 # +16
        self.vActMapping[85:92] = np.array(offs+np.arange(85,92)).astype(np.int32)
        offs += 2 + 3 # +21
        self.vActMapping[92:97] = np.array(offs+np.arange(92,97)).astype(np.int32)

        # print('self.vActMapping[0:5] >> ', self.vActMapping[0:5])
        fits.writeto(r'R:\ActuatorsMappingIndices.fits',self.vActMapping,overwrite=True)

        # for k in range(self.nact_lin_x*self.nact_lin_y):




        xdm, ydm = np.meshgrid(np.arange(dms)-dms/2, np.arange(dms)-dms/2)
        self.xdm = xdm.T.astype('float32')
        self.ydm = ydm.T.astype('float32')
        # print('self.xdm:')
        # print(self.xdm)
        print('self.dmmask:')
        # print(self.dmmask)
        # print('product:')
        # print(self.xdm * self.dmmask)
        # self.xdm = (self.xdm * self.dmmask).astype('float32')
        # self.ydm = (self.ydm * self.dmmask).astype('float32')
        self.xdm = (self.xdm).astype('float32')
        self.ydm = (self.ydm).astype('float32')

        # NORMALIZE THE MODES...
        avg_xdm = np.mean(self.xdm[self.xdm[:] != 0.])
        avg_ydm = np.mean(self.ydm[self.ydm[:] != 0.])
        # print('avg_x/ydm: ',avg_xdm,avg_ydm)

        self.xdm[self.xdm[:] != 0.] -= avg_xdm
        self.ydm[self.ydm[:] != 0.] -= avg_ydm

        self.xdm = self.xdm / np.sqrt(np.sum(self.xdm**2))
        self.ydm = self.ydm / np.sqrt(np.sum(self.ydm**2))

        # file = r'C:\Home\Downloads\dmmask.fits'
        # os.remove(file)
        # fits.writeto(file,self.dmmask)

    # -------------------------------------------------------------------------
    def list_2_map(self, data):
        ''' -------------------------------------------------------------------
        Convert a 1D list of 97 voltages back into a 2D map for an
        easier 2D representation.
        ------------------------------------------------------------------- '''
        dms = self.dms
        res = np.zeros((dms, dms), dtype=np.float32)

        i0 = 0
        for k in xrange(dms*dms):
            if (self.dmmask.flatten()[k] > 0):
                i = k / dms
                j = k % dms
                res[i,j] = data[i0]
                i0 += 1
        return(res)

    # -------------------------------------------------------------------------
    def zer_mode_bank_2D(self, i0, i1):
        ''' -------------------------------------------------------------------
        Produce a bank of Zernike modes for the ALPAO 97-15 DM.

        Parameters:
        ----------
        - i0: index of the first zernike mode to be added
        - i1: index of the last zernike mode to be included
        ------------------------------------------------------------------- '''
        dZ = i1 - i0 + 1
        res = np.zeros((dZ, self.dms, self.dms)).astype('float32')
        print('[zer_mode_bank_2D]: self.dms, self.dmmask >> ', self.dms, np.shape(self.dmmask))
        fits.writeto(r'C:\Users\lucas\Documents\dmmask_zer.fits',self.dmmask,overwrite=True)
        for i in range(i0, i1+1):
            res[i-i0] = mkzer1(i, self.dms, self.dms/2 +2) * self.dmmask
            avg = np.average(res[i-i0], weights=self.dmmask)
            # print('Average Zernike %d: %f' % (i, avg))
            res[i-i0] -= avg
            res[i-i0] = res[i-i0] * self.dmmask
            # avg = np.average(res[i-i0], weights=self.dmmask)
            # res[i-i0] = res[i-i0] / np.sqrt(np.sum(res[i-i0]*res[i-i0]))
            res[i-i0] = res[i-i0] / np.max(np.abs(res[i-i0]))
            # print('Average Zernike after normalization %d: %f %f' % (i, avg, np.sum(res[i-i0]*res[i-i0])))
            fname = r'C:\Users\lucas\Documents\zer_mode'+str(i)+'.fits'
            fits.writeto(fname,res[i-i0],overwrite=True)

        return(res)

    # -------------------------------------------------------------------------
    def poke_mode_bank_2D(self):
        ''' -------------------------------------------------------------------
        Produce a bank of poke modes for the ALPAO 97-15 DM
        ------------------------------------------------------------------- '''
        res = np.zeros((self.nba, self.dms, self.dms)).astype('float32')

        kk = 0
        for ii in range(self.dms):
            for jj in range(self.dms):
                if self.dmmask[jj,ii] > 0.5:
                    # print('jj ,ii >> ',jj,ii)
                    res[kk, ii,jj] = 1.0
                    kk += 1
        return res

# =============================================================================
# =============================================================================

class WFC(object):
    ''' -----------------------------------------------------------------------
    Generic wavefront control class
    ----------------------------------------------------------------------- '''
    def __init__(self, shm_cor_LO, shm_cor_HO, shm_cal, nav, simu=False):
        self.DM = DM()
        self.alpao_cor_LO = shm_cor_LO # correction
        self.alpao_cor_HO = shm_cor_HO # correction
        self.alpao_cal = shm_cal # calibration
        self.nav = nav           # number of calibration coadds
        self.a0 = 0.1            # calibration amplitude
        self.simu = simu

        self.shm_comb = shm('R:\\comb.im.shm', verbose=False)

        self.modes     = np.array([self.DM.xdm, self.DM.ydm]) # (ttilt)
        self.keepgoing = False
        self.gain      = 0.001 # default gain
        self.tsleep    = 1e-6 # for tests
        self.verbose   = False
        self.abort     = False
        self.calib_on  = False
        self.loop_on   = False
        self.avgtiming = 0.
        self.timingT   = 500

        self.ttx_mode = 0.
        self.tty_mode = 0.

        self.mode_norm = np.zeros(20)
        self.zer_mode_amp = np.zeros(20)
        self.zer_mode_amp_resid = np.zeros(20)

        self.mode_select = np.ones(20)
        self.bApplyCorrection = False

        self.log_length = 100000
        self.bLogSlopes = False
        self.bLogErrors = False
        self.bLogCmds = False
    # -------------------------------------------------------------------------
    # def orthonormalize_modes(self):
    #     # GRAM-SCHMIDT DECOMPOSITION
    #     KeptModes = 0
    #     # self.OrthoMode = self.modes
    #     size = self.modes.shape[1]
    #     NewMode = np.zeros([self.nmodes,size,size])
    #     OrthoMode = np.zeros([self.nmodes,size,size])
    #     OrthoNormalMode = np.zeros([self.nmodes,size,size])
    #     RemCoeffs = np.zeros([self.nmodes,size,size])
    #     NormeOrtho = np.zeros(self.nmodes)
    #     index_kept = (np.zeros(self.nmodes) - 1).astype('int')
    #
    #     for j in range(0,self.nmodes):
    # # 		OrthoMode[j,*,*] = Mode[j,*,*]
    #         OrthoMode[j] = self.modes[j]
    #
    #         for i in range(0,j):
    # # 			if (NormeOrtho[i] ne 0.D) then begin
    # # 				;ProjMap = total(conj(Mode[j,*,*])*OrthoMode[i,*,*],/DOUBLE) / total(conj(OrthoMode[i,*,*])*OrthoMode[i,*,*],/DOUBLE)
    # # 				ProjMap = total(Mode[j,*,*]*OrthoMode[i,*,*],/DOUBLE) / total(OrthoMode[i,*,*]*OrthoMode[i,*,*],/DOUBLE)
    # # 			endif else begin
    # # 				ProjMap = 0.D
    # # 			endelse
    #             if NormeOrtho[i] != 0.:
    #                 ProjMap = np.sum(self.modes[j]*OrthoMode[i]) / np.sum(OrthoMode[i]*OrthoMode[i])
    #             else:
    #                 ProjMap = 0.
    # # ;			pm, 'ProjMap >> ',j,i,ProjMap
    # # 			RemCoeffs[j,*,*] += ProjMap * OrthoMode[i,*,*]
    #             RemCoeffs[j] += ProjMap * OrthoMode[i]
    # # 		endfor
    # #
    # # 		;OrthoMode[j,*,*] = OrthoMode[j,*,*] - RemCoeffs[j,*,*]
    # # 		OrthoMode[j,*,*] = Mode[j,*,*] - RemCoeffs[j,*,*]
    #         OrthoMode[j] = self.modes[j] - RemCoeffs[j]
    #
    # # 		;InvNormeOrthoMode = total(conj(OrthoMode[j,*,*])*OrthoMode[j,*,*],/DOUBLE)
    # # 		InvNormeOrthoMode = total(OrthoMode[j,*,*]*OrthoMode[j,*,*],/DOUBLE)
    # # 		pm, 'InvNormeOrthoMode >> ',InvNormeOrthoMode
    #         InvNormeOrthoMode = np.sum(OrthoMode[j]*OrthoMode[j])
    #         print('InvNormeOrthoMode >> %f' % (InvNormeOrthoMode))
    #
    # # ;		window,10+2*j;,xsize=size,ysize=size
    # # ;		;tvscl, (OrthoMode[j,*,*] / NormeOrtho[j] + RemCoeffs[j,*,*])
    # # ;		shade_surf, Mode[j,*,*]
    # # ;
    # # ;		window,11+2*j;,xsize=size,ysize=size
    # # ;		;tvscl, (OrthoMode[j,*,*] / NormeOrtho[j] + RemCoeffs[j,*,*])
    # # ;		shade_surf, OrthoMode[j,*,*]
    # #
    # # 		if (InvNormeOrthoMode lt 1.D-5) then begin
    # # 			pm,'Orthogonal mode ',strcompress(j,/remove_all),' is close to NULL! Setting to zero.'
    # # 			;NormeOrtho[j] = dcomplex(0.D,0.D)
    # # 			NormeOrtho[j] = 0.D
    # # 			;OrthoNormalMode[j,*,*] = dcomplex(0.D,0.D)
    # # 			OrthoNormalMode[j,*,*] = 0.D
    # # 			;nModes--
    # # 		endif else begin
    # # 			NewMode[KeptModes,*,*] =  Mode[j,*,*]
    # # 			NormeOrtho[KeptModes] = 1.D / InvNormeOrthoMode
    # # 			OrthoNormalMode[KeptModes,*,*] = OrthoMode[j,*,*] * NormeOrtho[KeptModes]
    # # 			index_kept[KeptModes] = j
    # # 			pm, 'Kept mode > ',KeptModes,j
    # # 			KeptModes++
    # # 		endelse
    #         if InvNormeOrthoMode < 1.e-5:
    #             NormeOrtho[j] = 0.
    #             OrthoNormalMode[j,0:,0:] = 0.
    #         else:
    #             NewMode[KeptModes] =  self.modes[j]
    #             NormeOrtho[KeptModes] = 1. / InvNormeOrthoMode
    #             OrthoNormalMode[KeptModes] = OrthoMode[j] * NormeOrtho[KeptModes]
    #             index_kept[KeptModes] = j
    #             KeptModes += 1
    #
    # # 	endfor
    # #
    # #
    # # 	index_kept = index_kept[where(index_kept ne -1L)]
    #
    #     index_kept.compress((index_kept != -1).flat)
    #
    # #
    # # 	for k=0,KeptModes-1L do begin
    # # 		;pm, 'Kept orthogonal mode >> ',k,total(conj(OrthoNormalMode[k,*,*])*OrthoNormalMode[k,*,*],/DOUBLE)
    # # 		pm, 'Kept orthogonal mode >> ',k,total(OrthoNormalMode[k,*,*]*OrthoNormalMode[k,*,*],/DOUBLE)
    # # 	endfor
    #     for k in range(0,KeptModes):
    #         print('Kept orthogonal mode >> %d %f' % (k, np.sum(OrthoNormalMode[k]*OrthoNormalMode[k])))
    #         fname = r'C:\Home\Downloads\ortho_mode'+str(k)+'.fits'
    #         fits.writeto(fname,OrthoNormalMode[k],overwrite=True)
    #
    #     N_OrthoModes = KeptModes
    #
    # # 	pm, 'N_OrthoModes >> ',N_OrthoModes
    # #
    # # ;stop
    # # ;goto,next
    # #
    # # 	pm,''
    # # 	for j=0,N_OrthoModes-1 do begin
    # # 		for i=0,N_OrthoModes-1 do begin
    # # 			;proj_ortho = total(conj(OrthoNormalMode[i,*,*])*OrthoNormalMode[j,*,*]) / $
    # # 			;	total(conj(OrthoNormalMode[j,*,*])*OrthoNormalMode[j,*,*])
    # # 			proj_ortho = total(OrthoNormalMode[i,*,*]*OrthoNormalMode[j,*,*]) / $
    # # 				total(OrthoNormalMode[j,*,*]*OrthoNormalMode[j,*,*])
    # # ;			pm, 'proj_ortho >> ',j,i,proj_ortho
    # #
    # # 		endfor
    # # 	endfor
    #
    #     for i in range(0,self.nmodes):
    #         OrthoNormalMode[i] = OrthoNormalMode[i] / np.sqrt(np.sum(OrthoNormalMode[i]*OrthoNormalMode[i]))
    #         print('Check ortho normal %d %f' % (i,np.sum(OrthoNormalMode[i]*OrthoNormalMode[i])))
    #
    #     # SANITY CHECK
    #     bSanityCheck = False
    #     if bSanityCheck == True:
    #         for j in range(0,N_OrthoModes):
    #             for i in range(0,N_OrthoModes):
    #                 proj_ortho = np.sum(OrthoNormalMode[i]*OrthoNormalMode[j]) / np.sum(OrthoNormalMode[j]*OrthoNormalMode[j])
    #                 print('proj_ortho >> %d %d %f' % (j,i,proj_ortho))
    #
    #
    #
    #     proj_mode = np.zeros([N_OrthoModes,N_OrthoModes])
    #
    #     for j in range(0,N_OrthoModes):
    #         for i in range(0,N_OrthoModes):
    #             proj_mode[i,j] = np.sum(NewMode[i]*OrthoNormalMode[j])
    #             # print('proj_mode >> %d %d %f' % (i,j,proj_mode[i,j]))
    #
    #     print('proj_mode:')
    #     print(proj_mode)
    #
    # 	# for j=0,N_OrthoModes-1 do begin
    # 	# 	pm, j, '     ',Mode_Piston[index_kept(j),*]
    # 	# endfor
    #
    #
    # 	# index = where(abs(proj_mode) lt 1.D-6)
    #     index = np.asarray((np.abs(proj_mode).flat < 1.e-6).nonzero())
    #     print('index:')
    #     print(index)
    # 	# if (index[0] ne -1) then begin
    # 	# 	proj_mode(index) = 0.D
    # 	# endif
    #     for k in index:
    #         proj_mode.flat[k] = 0.
    #     # proj_mode[np.nonzero(proj_mode)] = 0.
    #
    # 	# inv_proj_mode = invert(proj_mode,/Double)
    #     try:
    #         inv_proj_mode = np.linalg.pinv(proj_mode)
    #     except np.linalg.LinAlgError:
    #         pass
    #
    #     self.proj_mode = proj_mode
    #     self.inv_proj_mode = inv_proj_mode
    #
    #     print('Orthonormal Projected Modes Matrix')
    #     print(self.proj_mode)
    #     print('Projected Modes Inverse Matrix')
    #     print(self.inv_proj_mode)
    #
    #     print('Sanity Check [unity]')
    #     SanityCheck = np.matmul(self.inv_proj_mode,self.proj_mode)
    #     print(SanityCheck)
    #
    #     return OrthoNormalMode

    def orthonormalize_resp_modes(self):
        # GRAM-SCHMIDT DECOMPOSITION
        KeptModes = 0
        # self.OrthoMode = self.modes
        size = self.vRespMode.shape[1]
        NewMode = np.zeros([self.nmodes,size])
        OrthoMode = np.zeros([self.nmodes,size])
        OrthoNormalMode = np.zeros([self.nmodes,size])
        RemCoeffs = np.zeros([self.nmodes,size])
        NormeOrtho = np.zeros(self.nmodes)
        index_kept = (np.zeros(self.nmodes) - 1).astype('int')

        for j in range(0,self.nmodes):
    		# OrthoMode[j,*,*] = Mode[j,*,*]
            OrthoMode[j] = self.vRespMode[j]
            fname = r'C:\Home\Downloads\vRespMode'+str(j)+'.fits'
            fits.writeto(fname,self.vRespMode[j],overwrite=True)

            for i in range(0,j):
                if NormeOrtho[i] != 0.:
                    ProjMap = np.sum(self.vRespMode[j]*OrthoMode[i]) / np.sum(OrthoMode[i]*OrthoMode[i])
                else:
                    ProjMap = 0.
                RemCoeffs[j] += ProjMap * OrthoMode[i]

            OrthoMode[j] = self.vRespMode[j] - RemCoeffs[j]

            InvNormeOrthoMode = np.sum(OrthoMode[j]*OrthoMode[j])
            print('InvNormeOrthoMode >> %f' % (InvNormeOrthoMode))

            if InvNormeOrthoMode < 1.e-5:
                NormeOrtho[j] = 0.
                OrthoNormalMode[j,0:,0:] = 0.
            else:
                NewMode[KeptModes] =  self.vRespMode[j]
                NormeOrtho[KeptModes] = 1. / InvNormeOrthoMode
                OrthoNormalMode[KeptModes] = OrthoMode[j] * NormeOrtho[KeptModes]
                index_kept[KeptModes] = j
                KeptModes += 1

        index_kept.compress((index_kept != -1).flat)

        for k in range(0,KeptModes):
            print('Kept orthogonal mode >> %d %f' % (k, np.sum(OrthoNormalMode[k]*OrthoNormalMode[k])))
            fname = r'C:\Home\Downloads\ortho_mode'+str(k)+'.fits'
            fits.writeto(fname,OrthoNormalMode[k],overwrite=True)

        N_OrthoModes = KeptModes

        for i in range(0,self.nmodes):
            OrthoNormalMode[i] = OrthoNormalMode[i] / np.sqrt(np.sum(OrthoNormalMode[i]*OrthoNormalMode[i]))
            print('Check ortho normal %d %f' % (i,np.sum(OrthoNormalMode[i]*OrthoNormalMode[i])))

        # SANITY CHECK
        bSanityCheck = True
        if bSanityCheck == True:
            for j in range(0,N_OrthoModes):
                for i in range(0,N_OrthoModes):
                    proj_ortho = np.sum(OrthoNormalMode[i]*OrthoNormalMode[j]) / np.sum(OrthoNormalMode[j]*OrthoNormalMode[j])
                    print('proj_ortho >> %d %d %f' % (j,i,proj_ortho))



        proj_mode = np.zeros([N_OrthoModes,N_OrthoModes])

        for j in range(0,N_OrthoModes):
            for i in range(0,N_OrthoModes):
                proj_mode[i,j] = np.sum(NewMode[i]*OrthoNormalMode[j])
                # print('proj_mode >> %d %d %f' % (i,j,proj_mode[i,j]))

        print('proj_mode:')
        print(proj_mode)

        index = np.asarray((np.abs(proj_mode).flat < 1.e-6).nonzero())
        print('index:')
        print(index)

        for k in index:
            proj_mode.flat[k] = 0.
        # proj_mode[np.nonzero(proj_mode)] = 0.

        try:
            inv_proj_mode = np.linalg.pinv(proj_mode)
        except np.linalg.LinAlgError:
            pass

        self.proj_mode = proj_mode
        self.inv_proj_mode = inv_proj_mode

        print('Orthonormal Projected Modes Matrix')
        print(self.proj_mode)
        print('Projected Modes Inverse Matrix')
        print(self.inv_proj_mode)

        print('Sanity Check [unity]')
        SanityCheck = np.matmul(self.inv_proj_mode,self.proj_mode)
        print(SanityCheck)

        return OrthoNormalMode

    # -------------------------------------------------------------------------
    def get_slopes(self, nav=20, reform=True):
        ''' -------------------------------------------------------------------
        test
        ------------------------------------------------------------------- '''
        # print('Entering \"get_slopes\"\n')

        cnt = self.shm_comb.get_counter()
        #time_start = time.time()
        tmp = self.shm_comb.get_data(check=cnt, reform=True)
        #print("Timing call get_data: %.3e" % (time.time() - time_start))
        self.clockA = time.perf_counter()
        cnt = self.shm_comb.get_counter()

        x_sig = tmp[0]
        y_sig = tmp[1]
        phot  = tmp[2]

       	for ii in range(nav-1):
            tmp = self.shm_comb.get_data(check=cnt, reform=True)
            cnt = self.shm_comb.get_counter()
            x_sig += tmp[0]
            y_sig += tmp[1]
            phot  += tmp[2]
        x_sig /= nav
        y_sig /= nav
        phot /= nav

        # print('Computing done!\n')

        x_sig[phot == 0] = 0.0
        y_sig[phot == 0] = 0.0

        if self.bLogSlopes == True and self.keepgoing == True:
            # print('signal dims >> ',np.shape(signal))
            self.LogDataSlopes[self.log_counter] = tmp
            self.LogTimeStamp[self.log_counter] = self.clockA - self.clock_zero
            # log_counter += 1

        return np.concatenate((x_sig.flatten(), y_sig.flatten()))

    # -------------------------------------------------------------------------
    def reset(self,):
        self.alpao_cor_LO.set_data(0.0 * self.alpao_cor_LO.get_data())
        self.alpao_cor_HO.set_data(0.0 * self.alpao_cor_HO.get_data())
        # self.alpao_cal.set_data(0.0 * self.alpao_cor.get_data())

    # -------------------------------------------------------------------------
    def abort(self,):
        print('WFC::abort()')
        self.alpao_cor_LO.set_data(0.0 * self.alpao_cor_LO.get_data())
        self.alpao_cor_HO.set_data(0.0 * self.alpao_cor_HO.get_data())
        # self.alpao_cal.set_data(0.0 * self.alpao_cor.get_data())
        self.keepgoing = False

    # -------------------------------------------------------------------------
    def calibrate(self, a0 = 0.1, reform=False):
        ''' -------------------------------------------------------------------
        Generic calibration procedure for set of modes attached to the object
        ------------------------------------------------------------------- '''
        self.calib_on = True
        self.abort = False
        self.a0 = a0

        dm0 = self.alpao_cal.get_data() # DM starting position

        phot = self.shm_comb.get_data()[2]

        RESP = [] # empty response matrix holder

        self.nmodes = self.modes.shape[0]
        self.mode_mult = np.zeros(self.nmodes)

        print('Calibrating with %d modes...' % self.nmodes)

        self.mode_norm = np.zeros(self.nmodes)
        self.zer_mode_amp = np.zeros(self.nmodes)
        self.zer_mode_amp_resid = np.zeros(self.nmodes)
        self.mode_select = np.ones(self.nmodes,dtype='int16')

        # MODIFY THE DM MODES INTO ORTHONORMAL MODES
        # self.orthoModes = (self.orthonormalize_modes()).astype('float32')
        self.orthoModesNorm = np.ones(self.nmodes)
        # for ii in range(self.nmodes):
        #     self.orthoModesNorm[ii] = np.sqrt(np.sum((self.orthoModes[ii] * a0)**2))
        # print(type(self.orthoModes))
        # print(type(self.modes))

        # bUseOrthoModes = True
        # if bUseOrthoModes == True:
        #     self.modes = self.orthoModes
            # go over the modes to be used for this WFC loop
            # for ii in range(self.nmodes):
            #     print('ii: %d' % ii)
            #     if self.abort:
            #         self.abort = False
            #         self.calib_on = False
            #         return
            #     self.alpao_cal.set_data((dm0 + self.orthoModes[ii] * a0))
            #     fname = r'C:\Home\Downloads\orthoMode_cal'+str(ii)+'.fits'
            #     fits.writeto(fname,(self.orthoModes[ii] * a0),overwrite=True)
            #     self.mode_norm[ii] = np.sqrt(np.sum((self.orthoModes[ii] * a0)**2))
            #     time.sleep(self.tsleep) # for simulation only!
            #     sys.stdout.write("\rmode %d  %.3e\n" % (ii,self.mode_norm[ii],))
            #     sys.stdout.flush()
            #
            #     RESP.append(self.get_slopes(self.nav, reform=reform))
        # else:
        # go over the modes to be used for this WFC loop

        print('self.DM.dms >> ',self.DM.dms)
        self.vRespMode = np.zeros([self.nmodes,2*(self.DM.dms-1)**2])
        for ii in range(self.nmodes):
            print('ii: %d' % ii)
            if self.abort:
                self.abort = False
                self.calib_on = False
                return
            self.alpao_cal.set_data((dm0 + self.modes[ii] * a0))
            fname = r'C:\Home\Downloads\mode'+str(ii)+'.fits'
            fits.writeto(fname,(self.modes[ii] * a0),overwrite=True)
            self.mode_norm[ii] = np.sqrt(np.sum((self.modes[ii] * a0)**2))
            # time.sleep(1.) # for simulation only!
            sys.stdout.write("\rmode %d  %.3e\n" % (ii,self.mode_norm[ii],))
            sys.stdout.flush()
            self.vRespMode[ii,:] = self.get_slopes(self.nav, reform=reform)
            # RESP.append(self.get_slopes(self.nav, reform=reform))
            RESP.append(self.vRespMode[ii,:])

        self.alpao_cal.set_data(dm0) # back to DM starting position

        self.bUseOrthoRespModes = False
        if self.bUseOrthoRespModes:
            self.RespOrthoModes = (self.orthonormalize_resp_modes()).astype('float32')

            self.RR = np.array(self.RespOrthoModes)
            fits.writeto(r'C:\Home\Downloads\RR.fits',self.RR,overwrite=True)
            #self.RR[np.abs(self.RR) < 1e-2] = 0.0 # matrix clean up
            self.RTR = self.RR.dot(self.RR.T)
        else:
            self.RR = np.array(RESP)
            fits.writeto(r'C:\Home\Downloads\RR.fits',self.RR,overwrite=True)
            #self.RR[np.abs(self.RR) < 1e-2] = 0.0 # matrix clean up
            self.RTR = self.RR.dot(self.RR.T)
        self.calib_on = False

    # -------------------------------------------------------------------------
    def reload_cal(self, fname=None):
        ''' -------------------------------------------------------------------
        Tries to reload a previously saved calibration?
        ------------------------------------------------------------------- '''

        try:
            self.RR = fits.getdata(fname)
        except:
            print("Calibration matrix %s not available" % (fname))
            return
        self.RRinv = pinv(self.RR.T, rcond=0.1)

    # -------------------------------------------------------------------------
    def cloop(self,):
        ''' -------------------------------------------------------------------
        Generic closed-loop procedure for the modes attached to the object.
        ------------------------------------------------------------------- '''

        self.keepgoing = True
        self.loop_on = True
        self.avgtiming = 0.
        counter = 0
        time_start = time.time()
        time_end = time_start
        # time_A = 0.
        # time_B = 0.
        # time_C = 0.
        # time_D = 0.
        # time_E = 0.
        # time_F = 0.
        # time_G = 0.
        #
        # time_tmp = time_start
        self.ttx_mode = 0.
        self.tty_mode = 0.

        self.zer_modes = np.zeros(self.nmodes)

        self.zer_filter = np.zeros(self.nmodes)
        self.zer_filter[0] = 1

        # self.bLogRTData = True
        # self.log_length = 10
        self.LogDataSlopes = np.zeros([self.log_length,3,10,10],dtype='float32')
        self.LogDataErrors = np.zeros([self.log_length,self.nmodes])
        self.LogDataCmds = np.zeros([self.log_length,self.nmodes])
        self.LogDataDMProjModes = np.zeros([self.log_length,self.nmodes])
        self.LogTimeStamp = np.zeros(self.log_length,dtype='float64')

        diff_time = 0.

        now = datetime.now()

        self.clock_zero = time.perf_counter()

        self.log_counter = 1

        self.time_A = 0.
        self.time_B = 0.
        self.time_C = 0.
        self.time_D = 0.
        self.time_E = 0.
        self.time_F = 0.
        self.time_G = 0.

        while self.keepgoing:
            # time_tmp = time.time()
            signal = self.get_slopes(1, reform=False) # WFS slopes
            # if self.bLogRTData == True:
            #     # print('signal dims >> ',np.shape(signal))
            #     self.LogDataSlopes[log_counter,0:1,:,:] = signal.reshape([10,10])
            #     self.LogDataSlopes[log_counter,2,:,:] =
            # time_A += time.time() - time_tmp

            # time_tmp = time.time()
            dm0 = self.alpao_cor.get_data()        # DM shape B4 correction
            # time_B += time.time() - time_tmp

            if self.bUseOrthoRespModes == True:
                # PORJECT ERROR SIGNAL signal ONTO THE ORTHONORMAL RESPONSE BASIS
                proj_ortho_ee = self.RR.dot(signal)
                # print('proj_ortho_ee >> ',proj_ortho_ee)
                # DEPROJECT THE RESULT TO TRANSLATE INTO THE INITIAL BASIS
                ee = np.dot(proj_ortho_ee,self.inv_proj_mode)
                # ee  = self.a0 * self.RRinv.T.dot(proj_ee)    # error signal
                # print('ee >> ',ee)
            else:
                time_tmp = time.time()
                ee  = self.RRinv.dot(signal)    # error signal

            try:
                cor = np.average(self.modes, weights=ee * self.mode_select[0:self.nmodes], axis=0)
            except:
                print('exception!!')
                cor = 0.0 * self.modes[0]

            if self.bLogErrors == True:
                self.LogDataErrors[self.log_counter,:] = ee
            # time_C += time.time() - time_tmp

            # time_tmp = time.time()
            # time_D += time.time() - time_tmp

            # time_tmp = time.time()
            cor *= self.nmodes * (ee * self.mode_select[0:self.nmodes]).sum()
            # time_E += time.time() - time_tmp

            # time_tmp = time.time()
            self.global_gain = 0.998
            dm1 = self.global_gain * (dm0 - self.gain * cor.astype('float32'))
            # time_F += time.time() - time_tmp

            # time_tmp = time.time()
            if self.bApplyCorrection == True:
                self.alpao_cor.set_data(dm1)
            self.time_G += time.time() - time_tmp


            # dm1 CONTAINS THE ABSOLUTE DM SHAPE FROM ITS 'ZERO' POSITION
            # PROJECT ON ZERNIKE TO GET EACH MODE CONTRIBUTION
            for k in range(self.nmodes):
                self.zer_mode_amp[k] = np.sum(self.modes[k] * dm1) / self.mode_norm[k] * self.a0

            if self.bLogCmds == True:
                self.LogDataDMProjModes[self.log_counter,:] = self.zer_mode_amp
                # self.zer_mode_amp[k] = np.sum(self.orthoModes[k] * dm1) / self.orthoModesNorm[k]
            self.zer_mode_amp_resid[0:ee.size] = ee
            # print('modes >> ')
            # print(self.zer_mode_amp)


            '''time.sleep(self.tsleep)'''

            # self.ttx_mode += ee[0]
            # self.tty_mode += ee[1]

            counter += 1
            if counter % self.timingT == 0:
                time_end = time.time()
                self.avgtiming = (time_end - time_start)/counter
                time_start = time_end

                print("cloop average timing/frequency: %.3e %.3e %.3e %d"% (self.avgtiming, 1./self.avgtiming, self.time_G/counter, self.log_counter))
                # print(self.modes)
                # print(cor)
                # print('tip: %.3e  tilt: .3e' % (cor[0], cor[1]))
                # print("cloop average timing: %.3e %.3e %.3e %.3e %.3e %.3e %.3e"% (time_A/counter,time_B/counter,time_C/counter,time_D/counter,time_E/counter,time_F/counter,time_G/counter))
                counter = 0
                # time_A = 0.
                # time_B = 0.
                # time_C = 0.
                # time_D = 0.
                # time_E = 0.
                # time_F = 0.
                self.time_G = 0.


                if self.verbose:
                    sys.stdout.write("\rcoeffs: "+ "%+.3f " * self.nmodes % (tuple(ee)))
                    sys.stdout.flush()

            self.log_counter += 1
            # print('self.log_counter: ',self.log_counter)
            if ((self.bLogSlopes == True) or (self.bLogCmds)) and (self.log_counter >= self.log_length):
                maxmode = self.nmodes + 1
                save_folder = now.strftime("%Y.%m.%d_%Hh%Mm%Ss") + '_NUMRECS' + str(self.log_length) + '_MAXMODE' + str(maxmode) + '_CALIAMP' + str(np.float32(self.a0)) + '_LOOPGAIN' + str(np.float32(self.gain)) + '_GLGAIN' + str(np.float32(self.global_gain))
                os.makedirs('C:\\Data\\' + save_folder)

                # WRITE THE USED MODES
                for k in range(0,self.nmodes):
                    fname = 'C:\\Data\\mode'+str(k)+'.fits'
                    fits.writeto(fname,(self.modes[k] * self.a0),overwrite=True)

                if self.bLogSlopes == True:
                    hdr = fits.Header()
                    hdr['NUMRECS'] = np.float32(self.log_length)
                    hdr['MAXMODE'] = maxmode
                    hdr['CALIAMP'] = np.float32(self.a0)
                    hdr['LOOPGAIN'] = np.float32(self.gain)
                    hdr['GLGAIN'] = np.float32(self.global_gain)
                    fits.writeto('C:\\Data\\' + save_folder + '\\SlopeX.fits',self.LogDataSlopes[:,0,:,:].reshape([self.log_length,100]),header=hdr,overwrite=True)
                    fits.writeto('C:\\Data\\' + save_folder + '\\SlopeY.fits',self.LogDataSlopes[:,1,:,:].reshape([self.log_length,100]),header=hdr,overwrite=True)
                    fits.writeto('C:\\Data\\' + save_folder + '\\Illum.fits',self.LogDataSlopes[:,2,:,:].reshape([self.log_length,100]),header=hdr,overwrite=True)

                if self.bLogCmds == True:
                    fits.writeto('C:\\Data\\' + save_folder + '\\MeasuredErrors.fits',self.LogDataErrors[:,:],header=hdr,overwrite=True)
                    fits.writeto('C:\\Data\\' + save_folder + '\\DMShapeModes.fits',self.LogDataDMProjModes[:,:],header=hdr,overwrite=True)

                fits.writeto('C:\\Data\\' + save_folder + '\\timestamps.fits',self.LogTimeStamp[:],header=hdr,overwrite=True)

                self.keepgoing = False
                break

        self.loop_on = False
        print("\nWFC loop opened.")

# =============================================================================
# =============================================================================

class TT_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Tip-tilt wavefront control data structure
    ----------------------------------------------------------------------- '''

    def __init__(self, simu=False):
        cor = shm('R:\\dmdisp2.im.shm', verbose=False) # correction
        cal = shm('R:\\dmdisp6.im.shm', verbose=False) # calibration
        nav = 5
        super(TT_WFC, self).__init__(cor, cal, nav, simu)
        print('self.DM.xdm >> ',self.DM.xdm)
        self.modes = np.array([self.DM.xdm, self.DM.ydm]) # (ttilt)

    def calibrate(self, a0 = 0.1, reform=False):
        print('TT_WFC::calibrate()')
        super(TT_WFC, self).calibrate(a0=a0, reform=reform)
        self.RRinv = pinv(self.RR.T, rcond=0.1)

    def cloop(self,):
        super(TT_WFC, self).cloop()

    def reset(self,):
        super(TT_WFC, self).reset()

    def abort(self,):
        print('TT_WFC::abort()')
        super(TT_WFC, self).abort()

    def reload_cal(self, fname=None):
        super(TT_WFC, self).reload_cal(fname=fname)

# =============================================================================
# =============================================================================

class ZER_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Zernike based wavefront control data structure.
    ----------------------------------------------------------------------- '''

    def __init__(self, iz0=4, iz1=10, simu=False):
        ''' -------------------------------------------------------------------
        Class constructor.

        Specifies calibration and correction channels as well as the basis
        of modes to use to control the wavefront.

        Parameters:
        ----------
        - iz0: first index of the Zernike modes to control (default = 4)
        - iz1: last  index of the Zernike modes to control (default = 11)
        ------------------------------------------------------------------- '''
        self.iz0 = iz0
        self.iz1 = iz1
        cor_LO = shm('R:\\DM_Buffer_02.im.shm', verbose=False) # correction
        cor_HO = shm('R:\\DM_Buffer_03.im.shm', verbose=False) # correction
        cal = shm('R:\\dmdisp7.im.shm', verbose=False) # calibration
        nav = 5
        super(ZER_WFC, self).__init__(cor_LO, cor_HO, cal, nav, simu)
        print('iz0/iz1 : %d %d' %(iz0,iz1))
        self.modes = self.DM.zer_mode_bank_2D(self.iz0, self.iz1)
        self.verbose = False
        self.iz0 = iz0
        self.iz1 = iz1

    def calibrate(self, a0 = 0.1, reform=False):
        print('iz0, iz1 >> ',self.iz0,self.iz1)
        self.modes = self.DM.zer_mode_bank_2D(self.iz0, self.iz1)
        super(ZER_WFC, self).calibrate(a0=a0, reform=reform)
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        fits.writeto(r'C:\Home\Downloads\RRinv.fits',self.RRinv,overwrite=True)
        fits.writeto(r'C:\Home\Downloads\RRinvRRT.fits',np.dot(self.RRinv,self.RR.T),overwrite=True)

    def cloop(self,):
        super(ZER_WFC, self).cloop()

    def reset(self,):
        super(ZER_WFC, self).reset()

    def reload_cal(self, fname=None):
        super(ZER_WFC, self).reload_cal(fname=fname)

# =============================================================================
# =============================================================================

class ZON_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Zonal based wavefront control data structure.
    ----------------------------------------------------------------------- '''

    def __init__(self, simu=False):
        ''' -------------------------------------------------------------------
        Class constructor.

        Specifies calibration and correction channels as well as the basis
        of modes to use to control the wavefront.

        ------------------------------------------------------------------- '''

        cor = shm('R:\\dmdisp3.im.shm', verbose=False) # correction
        cal = shm('R:\\dmdisp7.im.shm', verbose=False) # calibration
        nav = 5
        super(ZON_WFC, self).__init__(cor, cal, nav, simu)
        self.modes = self.DM.poke_mode_bank_2D()
        self.verbose = False

    def calibrate(self, a0 = 0.1, reform=False):

        super(ZON_WFC, self).calibrate(a0=a0, reform=reform)
        # filtering for the pseudo-inverse matrix
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        fits.writeto(r'C:\Home\Downloads\RRinv.fits',self.RRinv,overwrite=True)
        fits.writeto(r'C:\Home\Downloads\RRinvRRT.fits',np.dot(self.RRinv,self.RR.T),overwrite=True)

    def cloop(self,):
        super(ZON_WFC, self).cloop()

    def reset(self,):
        super(ZON_WFC, self).reset()

    def reload_cal(self, fname=None):
        super(ZON_WFC, self).reload_cal(fname=fname)

# =============================================================================
# =============================================================================

if __name__ == "__main__":
    wfc = TT_WFC()
