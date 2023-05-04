#!/usr/bin/env python3

import numpy as np
import os
import sys
import time
import pdb

aoc_home = os.getenv('HOME') + 'bin\\'
python_home = os.getenv('PYTHON_HOME')
#sys.path.append(home+'/src/python/libs/')

#from scexao_shm import shm

import astropy.io.fits as fits

from xaosim.wavefront import kolmo
from xaosim.shmlib import shm

dms = 11            # DM size in actuators (SQUARE ARRAY)
nAct = 97           # EFFECTIVE NUMBER OF ACTUATORS
wl  = 0.8           # wavelength in microns
c0 = wl/(2.0*np.pi) # from phase to disp in microns

xsz = 20*dms
ysz = 30*dms
r1 = np.random.randn(xsz, ysz)
r2 = np.random.randn(xsz, ysz)

argvs = sys.argv[1:]

help_message = '''
===============================================
          !! Turbulence simulation !!

Adds a simulated Kolmogorov phase screen to the
Channel #1 of the deformable mirror.

===============================================

Arguments of the command:
------------------------

turbul   TOTAL_RMS  TIME_DELAY  AO_CORREC

examples:
--------

>> turbul 4.0 0.1 1.0

Introduces 4.0 ADU RMS (magnitude not calibrated!)
Refreshes every 0.1 seconds
Performs no AO correction (factor = 1.0)

>> turbul 4.0 1.0 10.0

Introduces 4.0 ADU RMS
Refreshes every 1.0 seconds
Performs an AO correction by reducing low spatial
  frequencies by a factor 10.0

>> turbul 0.0
The fastest way to remove any turbulence previously
applied to the DM!
===============================================
'''

# LOAD THE MAPPING FILE TO CONVERT ACTUATOR NUMBER (97) INTO A SQUARE 11x11 IMAGE
vActMapping = list(fits.getdata(aoc_home + 'ActuatorsMappingIndices.fits'))

# # ------------------------------------------------------------------
# #  another short hand to convert numpy array into image for display
# # ------------------------------------------------------------------
# def arr2im(arr, vmin=-0.5, vmax=0.5):
#     # print(arr.shape)
#     sqarr_disp = np.zeros(dms*dms)
#     sqarr_disp[vActMapping] = arr

#     sqarr_disp = sqarr_disp.reshape((dms,dms))

#     # print('sqarr_disp >> ')
#     # print(sqarr_disp, sqarr_disp.shape)
#     # print('arr >> ')
#     # print(arr, arr.shape)

#     mmin,mmax = sqarr_disp.min(), sqarr_disp.max() # vmin, vmax #

#     img = Image.fromarray(sqarr_disp.astype('float'))
#     rimg = img.resize((zoom*dms, zoom*dms))
#     rarr = np.asarray(rimg)

#     arr2 = rarr.copy()

#     arr2 -= mmin
#     arr2 /= (mmax-mmin+1e-9)
#     test = mycmap(arr2)
#     return((255*test[:,:,:3]).astype('int'))


try:
    test = argvs[0]
except:
    print(help_message)

# ------------------------------------
try:
    myrms = float(argvs[0])
except:
    myrms = 4.0
    print("default rms value: %.3f" % (myrms,))

# ------------------------------------
try:
    mytime = float(argvs[1])
except:
    mytime = 0.1
    print("time delay between iterations: %.3f seconds" % (mytime,))

# ------------------------------------
try:
    mycorrec = float(argvs[2])
except:
    mycorrec = 1.0
    print("No partial AO correction added")

# ------------------------------------

kolm = 1.0*c0 * kolmo(r1, 5.0, 2, correc=mycorrec, rms=myrms)

#kolm.dtype = np.float32
print("Kolmo screen std = %.1f" % (kolm.std()))

#disp = shm('/dev/shm/dmdisp1.im.shm', data=kolm[:dms,:dms].astype(np.float32),
#           verbose=False)
# disp = shm('R:\\dmdisp1.im.shm', data=kolm[:dms,:dms].astype(np.float32),
#            verbose=False)

voltageMapKolmo = (kolm[:dms,:dms].flatten())[vActMapping]
print('voltageMapKolmo >> ', voltageMapKolmo.shape)
disp = shm('R:\\DM_Buffer_01.im.shm', data=voltageMapKolmo.astype(np.float32),
           verbose=False)

xdm, ydm = np.meshgrid(np.arange(dms)-dms//2, np.arange(dms)-dms//2)
dmmask   = np.ones((dms, dms), dtype=np.int) # corner mask
dmmask[np.abs(xdm) + np.abs(ydm) > 7] = 0.0

kolm2 = np.tile(kolm,(2,2))
kolm2 = kolm2[:xsz+2*dms,:ysz+2*dms]

test = shm('R:\\turb_test.im.shm', data=kolm2.astype(np.float32), verbose=False)
# test = shm('R:\\DM_Buffer_01.im.shm', data=kolm2.astype(np.float32), verbose=False)

# -------------------------------------------

keepgoing = True

offx, offy = 0, 0
while keepgoing:

    offx += 1
    offy += 3
    offx = offx % xsz
    offy = offy % ysz

    #kolm2 = np.tile(kolm,(2,2))
    #kolm2 = kolm2[:xsz+2*dms,:ysz+2*dms]
    k1 = kolm2[offx:offx+dms, offy:offy+dms].copy()
    k1 -= k1.mean()
    mydata = (k1*dmmask).flatten()[vActMapping].astype(np.float32)
    print('(k1*dmmask).flatten()[vActMapping] >> ', mydata.shape)
    disp.set_data(mydata)

    # -------------------------------------------------
    # uncomment the following three lines to visualize
    # how the kolmogorov phase screen is scanned, using
    # >> shmview /tmp/turb_test.im.shm
    # ...
    # removed by default because CPU intensive for
    # nothing
    # -------------------------------------------------

    # kolm2[offx:offx+dms,offy:offy+dms] *= 2.0
    # test.set_data(kolm2.astype(np.float32))
    # kolm2[offx:offx+dms,offy:offy+dms] /= 2.0
    # -----------------------------------------
    time.sleep(mytime)

    if os.path.exists('.\\turb_stop'):
        keepgoing = False
        os.remove('.\\turb_stop')
