#!/usr/bin/env python

from email import header
import numpy as np
import os, sys, time
import struct
import astropy.io.fits as fits

home = os.getenv('HOME')
print('HOME: %s' % (home))
sys.path.append(home+'\\bin\\')

ciao_home = os.getenv('CIAO_HOME')
sys.path.append(ciao_home+'\\libs\\')

import xaosim.zernike as zernike
from xaosim.shmlib import shm

if (8 * struct.calcsize("P")) == 32:
    print("Use x86 libraries.")
    sys.path.append(ciao_home+"\\libs\\Lib\\")
    from Lib.asdk import DM
else:
    print("Use x86_64 libraries.")
    sys.path.append(ciao_home+"\\libs\\Lib64\\")
    from Lib64.asdk import DM

#from asdk import DM

# ==============================
#      CONNECT TO THE DM
# ==============================
seri_num = "BOL115"              # DM serial number
alpao    = DM(seri_num)          # the reference to the DM
nbAct    = int( alpao.Get('NBOfActuator') )
alpao.Send([0.] * nbAct)         # init the DM: zeroes everywhere
dms = 11                         # DM size in elements over the pupil (actual is 11)

# ==============================
#      CREATE A DM MAP
# ==============================
mask = np.ones((dms, dms), dtype=np.float32) # to remove non-existing actuators
mask[:3,0] = np.nan
mask[0,:3] = np.nan
mask[1, 1] = np.nan
mask[-3:,0] = np.nan
mask[-1,:3] = np.nan
mask[-2,1]  = np.nan
mask[0,-3:] = np.nan
mask[1,-2:] = np.nan
mask[2,-1]  = np.nan
mask[-1,-3:] = np.nan
mask[-2,-2:] = np.nan
mask[-3,-1]  = np.nan

flat_mask = mask.flatten()

# ================================================
#  SETUP MULTI-CHANNEL SHARED MEMORY DM STRUCTURE
# ================================================

nch     = 8  # number of DM channels
dmd0    = np.zeros((dms,  dms), dtype=np.float32) # initial dm disp map
sm_cntr = np.zeros(nch+1) - 1                     # shared mem counters for channels

# create/access shared memory data structures
# -------------------------------------------

for i in range(nch):
    exec("disp%d = shm(fname='R:\\dmdisp%d.im.shm', data=dmd0, verbose=False)" % (i,i))

disp = shm('R:\\ciao_dm_disp.im.shm',  data=dmd0, verbose=False)


time.sleep(0.01)
keepgoing = True

# ===============================================================
def get_counter_channel(chn):
    ''' ----------------------------------------
    Return the current channel counter value.
    Reads from the already-opened shared memory
    data structure.
    ---------------------------------------- '''
    global cnt;
    if chn < 8:
        exec('global cnt; cnt = disp%d.get_counter()' % (chn,))
    else: # chn == 8:
        cnt = disp.get_counter()
    return cnt

# ===============================================================
def map2dm(cmap):
    ''' --------------------------------------------
    Sends a 2D command to the DM, using the mapping
    knowledge provided by the flat_mask array.
    -------------------------------------------- '''
    global flat_mask, nbAct, alpao
    temp   = cmap.flatten()
    values = [0.0] * nbAct # data to be uploaded to the DM

    jj = 0
    for i in range(cmap.size):
        if not np.isnan(flat_mask[i]):
            values[jj] = temp[i]
            jj += 1

    alpao.Send(values)

    # print('Start test loop alpao.Send()')
    # # TEST
    # for i in range(10000):
    #     alpao.Send(values)
    # print('End test loop alpao.Send()')

    return()

# ===========================================================
def list2map2D(data):
    ''' ------------------------------------------------
    Convert a list of 97 voltages back into a 2D map for
    easy 2D representation.
    ------------------------------------------------ '''
    global flat_mask, nbAct
    res = np.zeros((dms,  dms), dtype=np.float32)

    i0 = 0
    for k in range(dms*dms):
        if not np.isnan(flat_mask[k]):
            i = k // dms
            j = k % dms
            res[i,j] = data[i0]
            i0 += 1
            print('i,j,res[i,j] >> ',i,j,res[i,j])
    return(res)

# ===========================================================
def map2D2list(data2D):
    ''' ------------------------------------------------
    Convert a 2D map for easy 2D representation back
    into a list of 97 voltages.
    ------------------------------------------------ '''
    global flat_mask, nbAct, alpao
    temp   = data2D
    values = [0.0] * nbAct # data to be uploaded to the DM

    i0 = 0
    for k in range(dms*dms):
        if not np.isnan(flat_mask[k]):
            i = k // dms
            j = k % dms
            values[i0] = temp[i,j]
            i0 += 1

    return(values)

# PREVIOUS CODE
data = np.loadtxt(home + '\\bin\\zygo_flat.txt')

# zer_map = fits.getdata(ciao_home + r'data\zer_map.fits')
# hdulist = fits.open(ciao_home + r'data\zer_map.fits')
# zer_map = hdulist[0].data
# hdulist.close()

# zer_data = shm(r'C:\Users\AOC\bin\zer_dmdisp4.im.shm',verbose=False)
zer_data = shm(r'C:\Users\AOC\bin\zer_dmdisp4.im.shm',verbose=False)
# zer_data.set_data(zer_map)

act_data = shm(r'C:\Users\AOC\bin\act_dmdisp5.im.shm',verbose=False)
act_data.save_as_fits(r'C:\Users\AOC\bin\flat_map_correction_11x11.fits')

fits.writeto(r'C:\Users\AOC\bin\flat_map_correction_97.fits',np.array(map2D2list(act_data.get_data()), ndmin=2).astype(np.float32), overwrite=True)
hdul = fits.open(r'C:\Users\AOC\bin\flat_map_correction_97.fits')
hdr = hdul[0].header
hdr.set('ISDMMAP',1)
hdul.close()
fits.writeto(r'C:\Users\AOC\bin\flat_map_correction_97.fits',np.array(map2D2list(act_data.get_data()), ndmin=2).astype(np.float32), header=hdr, overwrite=True)

disp0.set_data(list2map2D(data))
disp4.set_data(zer_data.get_data())
disp5.set_data(act_data.get_data())

# zer_data.close()
# act_data.close()

# bUpdateFlat = False
# updated_flat_file = home + '\\bin\\2021.03.01_zygo_flat.txt'

# if bUpdateFlat == True:
#     # set the DM with measured Zernike and pokes offsets
#     zer_data = shm(r'C:\Users\AOC\bin\zer_dmdisp4.im.shm',verbose=False)
#     act_data = shm(r'C:\Users\AOC\bin\act_dmdisp5.im.shm',verbose=False)

#     # THESE FILES CONTAIN THE 11x11 2D MAP
#     myZer = zer_data.get_data()
#     myAct = act_data.get_data()
#     # THIS FILE ONLY CONTAINS THE 97 ACTUATORS VALUES
#     data = np.loadtxt(home + '\\bin\\zygo_flat.txt')

#     flat_map = list2map2D(data) + myZer + myAct
#     np.savetxt(updated_flat_file, map2D2list(flat_map))
# else:
#     data = np.loadtxt(updated_flat_file)
#     flat_map = list2map2D(data)

# disp0.set_data(flat_map)

updt = True
# ===========================================================
#                  MAIN EVENT LOOP
# ===========================================================
while keepgoing:

    for i in range(nch+1):
        test = get_counter_channel(i)
        if test != sm_cntr[i]:
            sm_cntr[i] = test
            updt = True

    if updt == True:
        updt = False
        global combi;
        combi = np.zeros_like(disp.get_data())
        for i in range(nch):
            exec("global combi; combi += disp%d.get_data()" % (i,))

        disp.set_data(combi)
        map2dm(combi)
        time.sleep(0.001)

    if os.path.exists(home+'\\bin\\dm_stop'):
        print('Server is about to stop!\n')
        keepgoing = False
        os.remove(home+'\\bin\\dm_stop')

    time.sleep(0.001)

# ===========================================================
#                  FINISH THINGS CLEANLY
# ===========================================================
alpao.Reset()
sys.exit()
