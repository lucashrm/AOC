#!/usr/bin/env python

import pygame, sys
from pygame.locals import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import pyfits as pf
import astropy.io.fits as pf
import pdb
import threading
import mmap
import struct
import os

import time
import datetime
from scipy.interpolate import griddata
from scipy.signal      import medfilt2d
from numpy.linalg      import solve
from scipy.ndimage.interpolation import rotate

# =================================================
# =================================================


home = os.getenv('HOME')
ciao_home = os.getenv('CIAO_HOME')

#sys.path.append(home+"/src/python/libs/")
#import zernike
#from scexao_shm import shm

import xaosim.zernike as zernike
from xaosim.shmlib import shm

# =================================================
# =================================================


hmsg = """
-------------------------------------------------------
           _____               _ _
          |__  /___ _ __ _ __ (_) | _____
            / // _ \ '__| '_ \| | |/ / _ \
           / /|  __/ |  | | | | |   <  __/
          /____\___|_|  |_| |_|_|_|\_\___|


     Control the ALPAO DM using Zernike modes

Frantz.
-------------------------------------------------------
"""

plt.ion()
plt.show()

# ------------------------------------------------------------------
#                       global variables
# ------------------------------------------------------------------
zoom   = 1                  # zoom for the display of the images
mycmap = cm.jet             # colormap
dms    = 11                 # DM size
z0, z1 = 4, 17              # Zernike indices
z0, z1 = 2, 15              # Zernike indices
dz = z1 - z0 + 1

wmodes = ["MANUAL", "AUTO", "RECAL"]
umode = 0
wmode = wmodes[umode]
zreso = 0.001 # resolution for the amplitude of the modes
amin, amax = -1.0, 1.0#0.5, 0.5

z_vals = np.zeros(dz)      # zernike coefficients

sync = True
# -----------------------


# -----------------------
#   set up the window
# -----------------------
adelay = 200 # action delay (keys and joystick actions)
pygame.init()
pygame.key.set_repeat(adelay, 50)
last_joy_updt = pygame.time.get_ticks()
joy_pressed = 0

FPS = 10                        # frames per second setting
fpsClock = pygame.time.Clock()  # start the pygame clock!
XW, YW = 500, 900               # window dimensions
BW = 5                          # border width

screen = pygame.display.set_mode((XW, YW), 0, 32)
pygame.display.set_caption('Zernike')

chn = 4                         # DM channel for this program
disp = shm('R:\\dmdisp%d.im.shm' % (chn,), verbose=False)
# save_disp = shm(r'C:\Home\Programmation\Projets\CIAO_win\data\dmdisp%d.im.shm' % (chn,), verbose=False)

# ===============================================================
def combine_modes():
    # print('[combine_modes] z_vals = ',z_vals)
    a = np.zeros((dms, dms), dtype='float32')
    if sync:
        for i in range(z0, z1+1):
            a += z_vals[i-z0] * zernike.mkzer1(i, dms, 6, True)#False)
        #print a.max(), a.min()
        disp.set_data(a)

# ===============================================================
def draw_cursor(surf, t1, t2, active):
    xc, yc, hh, ww = t1
    vmin, vmax, vcur = t2
    x0, y0 = xc - ww/2, yc - hh/2
    x1 = x0 + ww * (vcur-vmin) / (vmax-vmin) - ww/20
    y1 = yc-hh

    if (active):
        BTCOL = GOLD
    else:
        BTCOL = BLUE
    pygame.draw.rect(surf, BTCOL, (x0,y0,ww, hh), 0)
    pygame.draw.rect(surf, FGCOL, (x0,y0,ww, hh), 1)

    pygame.draw.rect(surf, BTCOL, (x1,y1,hh, 2*hh), 0)
    pygame.draw.rect(surf, FGCOL, (x1,y1,hh, 2*hh), 1)

    hdl = font1.render("%3d nm" % (1000 * vcur,),  True, FGCOL)

    zn = hdl.get_rect()
    zn.center = (xc, yc)
    surf.blit(hdl, zn)

# ------------------------------------------------------------------
#              !!! now we are in business !!!!
# ------------------------------------------------------------------

WHITE = (255, 255, 255)
GREEN = (  0, 153,  76)
BLUE  = (  0,   0, 153)
RED   = (255,   0,   0)
BLK   = (  0,   0,   0)
DDGRY = ( 50,  50,  50)
DGRAY = ( 80,  80,  80)
GOLD  = (255, 215,   0)
ORANG = (255, 128,   0)

FGCOL = WHITE  # foreground color (text)
BGCOL = BLK    # background color
BTCOL = BLUE   # *button* color

background = pygame.Surface(screen.get_size())
background = background.convert()
background.fill(BLK)


# ----------------------------
#          labels
# ----------------------------
font1 = pygame.font.SysFont("default",   28)
font2 = pygame.font.SysFont("monospace", 16)

hdline1 = font1.render("MANUAL", True, FGCOL)#, DDGRY)
hdline2 = font1.render("AUTO",   True, FGCOL)
hdline3 = font1.render("RECAL",  True, FGCOL)

zone1 = hdline1.get_rect()
zone2 = hdline2.get_rect()
zone3 = hdline3.get_rect()

ZH = zone1.h                    # text height
FW, FH = XW-2*BW, YW-2*BW-2*ZH  # frame dimensions

zone1.center = (BW+FW/6, BW+ZH)
zone2.center = (BW+FW/2, BW+ZH)
zone3.center = (BW+5*FW/6, BW+ZH)

pygame.draw.rect(screen, DDGRY, (BW,BW,FW/3, 2*ZH), 0)
pygame.draw.rect(screen, DDGRY, (BW+FW/3,BW,FW/3, 2*ZH), 0)
pygame.draw.rect(screen, DDGRY, (BW+2*FW/3,BW,FW/3, 2*ZH), 0)

screen.blit(hdline1, zone1)
screen.blit(hdline2, zone2)
screen.blit(hdline3, zone3)

# -------------------------------
#        rectangles
# -------------------------------

blabels = ['SYNC', 'RESET', 'MEMO', 'SAVE', 'LOAD', 'RANDOM']
bxs = [2*BW+FW/6, 2*BW+FW/2, 2*BW+5*FW/6, 2*BW+FW/6, 2*BW+FW/2, 2*BW+5*FW/6]
by0 = YW - 4*ZH
bys = [by0, by0, by0, by0 + 2*ZH, by0 + 2*ZH, by0 + 2*ZH]

for i in range(6):
    exec('b_lbl%d = font1.render(blabels[i], True, FGCOL)' % (i,))
    exec('b_rec%d = b_lbl%d.get_rect()' % (i,i,))
    exec('b_rec%d.center = (bxs[i], bys[i]+0.75*ZH)' % (i,))
    exec('screen.blit(b_lbl%d, b_rec%d), ' % (i,i,))

xws, yws = 100, 100
imin, imax = 0, 0
surf_live_img   = pygame.surface.Surface((xws, yws))
rect1           = surf_live_img.get_rect()
rect1.center    = (5+xws/2, 5+yws/2)


zlabel = ['focus', 'astig1', 'astig2', 'coma1', 'coma2',
          'trfl1', 'trfl2', 'spheric', 'quad1', 'quad2',
          'secast1', 'secast2', 'unknown1', 'unknown2']

zlabel = ['tip', 'tilt', 'focus', 'astig1', 'astig2', 'coma1', 'coma2',
          'trfl1', 'trfl2', 'spheric', 'quad1', 'quad2',
          'secast1', 'secast2', 'unknown1', 'unknown2']

mode_actv = 0 # index to keep track of what mode is actuated

# =======================================================
# =======================================================
while True: # the main game loop
    clicked = False

    # =====================================
    for event in pygame.event.get():

        if event.type == QUIT:
            pygame.quit()

            # close shared memory access
            # --------------------------
            print("The program has ended normally.")
            disp.close()
            sys.exit()

        elif event.type == KEYDOWN:
            mmods = pygame.key.get_mods()

            if event.key == K_s:
                print("want to save?")
                if (mmods & KMOD_LCTRL):
                    #saveim = True
                    print("saving!")

            if event.key == K_ESCAPE:
                pygame.quit()
                # close shared memory access
                # --------------------------
                print("The program has ended normally.")
                disp.close()
                sys.exit()

            if event.key == K_c:
                mmods = pygame.key.get_mods()

                if (mmods & KMOD_LSHIFT):
                    t = threading.Thread(target=calibrate_wfs, args=())
                    t.start()
                else:
                    plot_cross = True - plot_cross

            # ============= ARROWS ==================
            if event.key == K_UP:
                if (umode == 0):
                    mode_actv = (mode_actv - 1) % dz

            if event.key == K_DOWN:
                if (umode == 0):
                    mode_actv = (mode_actv + 1) % dz

            if event.key == K_RIGHT:
                if (umode == 0):
                    mult = 1.0
                    if (mmods & KMOD_LSHIFT):
                        mult = 10.0
                    temp = z_vals[mode_actv] + mult * zreso
                    if (temp < amax):
                        z_vals[mode_actv] = temp
                    combine_modes()

            if event.key == K_LEFT:
                if (umode == 0):
                    mult = 1.0
                    if (mmods & KMOD_LSHIFT):
                        mult = 10.0
                    temp = z_vals[mode_actv] - mult * zreso
                    if (temp > amin):
                        z_vals[mode_actv] = temp
                    combine_modes()
            # ========================================

            if event.key == K_SPACE:
                toto = True

            if event.key == K_h:
                print(hmsg)

            # ------ TAB used to switch modes ------
            if event.key == K_TAB:
                umode = (umode + 1) % 3
                wmode = wmodes[umode]

        # ---------------------------------------------
        #             joystick buttons
        # ---------------------------------------------
        elif event.type == pygame.JOYBUTTONDOWN:
            if event.button == 5:
                umode = (umode + 1) % 3
                wmode = wmodes[umode]

            if event.button == 4:
                umode = (umode - 1) % 3
                wmode = wmodes[umode]

            if event.button == 8: # SYNCHRONIZE (SELECT BUTTON)
                if (umode == 0):
                    sync = True - sync
                    combine_modes()

            if event.button == 9: # RESET (START BUTTON)
                z_vals = np.zeros(dz)
                combine_modes()

        # ---------------------------------------------
        #             joystick axis motion
        # ---------------------------------------------
        elif event.type == pygame.JOYAXISMOTION:
            joy.compute()
            print(joy.x1, joy.y1)
            print("coucou 1")

        # ---------------------------------------------
        #             joystick hats
        # ---------------------------------------------
        elif event.type == pygame.JOYHATMOTION:
            #if event.value == (0,0):
                #last_joy_updt = pygame.time.get_ticks()
                #joy_pressed = 0
            if event.value == (1,0):
                if (umode == 0):
                    temp = z_vals[mode_actv] + 10 * zreso
                    if (temp < amax):
                        z_vals[mode_actv] = temp
                    combine_modes()

            elif event.value == (-1,0):
                if (umode == 0):
                    temp = z_vals[mode_actv] - 10 * zreso
                    if (temp > amin):
                        z_vals[mode_actv] = temp
                    combine_modes()

            elif event.value == (0,-1):
                if (umode == 0):
                    mode_actv = (mode_actv + 1) % dz

            elif event.value == (0, 1):
                if (umode == 0):
                    mode_actv = (mode_actv - 1) % dz

        # ---------------------------------------------
        #               MOUSE actions
        # ---------------------------------------------
        elif event.type == MOUSEMOTION:
            mx, my = event.pos

        elif event.type == MOUSEBUTTONDOWN:
            mx, my = event.pos
            clicked = True

            # --------- CLICK used to switch modes ------
            for i in range(3):
                exec('mytest = zone%d.collidepoint(mx, my)' % (i+1,))
                if (mytest):
                    umode = i
                    wmode = wmodes[umode]

            # ---------- manual mode use case ----------
            if (umode == 0):
                for i in range(6):
                    exec('mytest = b_rec%d.collidepoint(mx, my)' % (i,))
                    if (mytest):
                        if (i == 0):   # SYNC
                            sync = True - sync
                            combine_modes()

                        elif (i == 1): # RESET
                            z_vals = np.zeros(dz)
                            combine_modes()
                        elif (i == 2):
                            print("NOT IMPLEMENTED YET") # MEMO
                        elif (i == 3): # SAVE
                            # print("NOT IMPLEMENTED YET")
                            disp.save_as_fits(ciao_home + r'data\zer_map.fits')
                        elif (i == 4): # LOAD
                            print("NOT IMPLEMENTED YET")
                        else:          # EMPTY
                            z_vals = 0.3 * amax * np.random.randn(dz)
                            combine_modes()

    # ---------------------------------------------
    pygame.draw.rect(screen, DDGRY, (BW,BW+2*ZH,FW, FH), 0)
    pygame.draw.rect(screen, FGCOL, (BW,BW+2*ZH,FW, FH), 2)

    pygame.draw.rect(screen, DGRAY, (BW,BW,FW, 2*ZH), 0)
    pygame.draw.rect(screen, FGCOL, (BW,BW,FW, 2*ZH), 1)


    # ============== ACTIVATE THE MODE FUNCTIONALITIES =================
    if wmode == "AUTO":
        # ---------------------------------------------------------
        pygame.draw.rect(screen, DDGRY, (BW+FW/3,BW,FW/3, 2*ZH), 0)
        pygame.draw.rect(screen, FGCOL, (BW+FW/3,BW,FW/3, 2*ZH), 2)
        # ---------------------------------------------------------

    elif wmode == "RECAL":
        # ---------------------------------------------------------
        pygame.draw.rect(screen, DDGRY, (BW+2*FW/3,BW, FW/3, 2*ZH), 0)
        pygame.draw.rect(screen, FGCOL, (BW+2*FW/3,BW, FW/3, 2*ZH), 2)
        # ---------------------------------------------------------

    else:
        # ---------------------------------------------------------
        pygame.draw.rect(screen, DDGRY, (BW,BW,FW/3, 2*ZH), 0)
        pygame.draw.rect(screen, FGCOL, (BW,BW,FW/3, 2*ZH), 2)

        for i in range(z0, z1+1):
            yy = 100 + 50 * (i - z0)
            hdl = font1.render(zlabel[i-z0],  True, FGCOL)

            zone = hdl.get_rect()
            zone.center = (BW+FW/6, yy)
            screen.blit(hdl, zone)

            active = False
            if mode_actv == i-z0:
                active = True
            draw_cursor(screen, (BW+FW/2, yy, 22, 200),
                        (amin, amax, z_vals[i-z0]), active)

        if sync:
            pygame.draw.rect(screen, GREEN,
                             (bxs[0]-FW/6, bys[0], FW/3-2*BW, 1.5*ZH), 0)
        else:
            pygame.draw.rect(screen, RED,
                             (bxs[0]-FW/6, bys[0], FW/3-2*BW, 1.5*ZH), 0)

        pygame.draw.rect(screen, FGCOL,
                         (bxs[0]-FW/6, bys[0], FW/3-2*BW, 1.5*ZH), 2)
        screen.blit(b_lbl0, b_rec0)

        for i in range(1,6):
            pygame.draw.rect(screen, ORANG,
                             (bxs[i]-FW/6, bys[i], FW/3-2*BW, 1.5*ZH), 0)
            pygame.draw.rect(screen, FGCOL,
                             (bxs[i]-FW/6, bys[i], FW/3-2*BW, 1.5*ZH), 2)
            exec('screen.blit(b_lbl%d, b_rec%d), ' % (i,i,))

        # ---------------------------------------------------------

    # ============= FINISH DRAWING THE WINDOW + TABS ==============
    screen.blit(hdline1, zone1)
    screen.blit(hdline2, zone2)
    screen.blit(hdline3, zone3)


    pygame.display.flip()

    fpsClock.tick(FPS)

pygame.quit()
sys.exit()
