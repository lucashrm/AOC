import astropy.io.fits as fits
import numpy as np
import ctypes
import _ctypes
import time

class testClass():
    def thread_call():
        lib = ctypes.CDLL(r'C:\RSI\SourcesC\CIAO_win\DEV_CIAO\x64\Release\DLLfromPython.dll')
        lib.calc_SH_data.restype = ctypes.c_int
        lib.calc_SH_data.argtypes = [ctypes.c_int,
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32),
            ctypes.c_int,
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.float32),
            np.ctypeslib.ndpointer(dtype=np.float32)]
            
        lib.trigger_calc_SH_data.restype = ctypes.c_int
        lib.trigger_calc_SH_data.argtypes = [ctypes.c_int,
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32),
            ctypes.c_int,
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.float32),
            np.ctypeslib.ndpointer(dtype=np.float32)]
            
        lib.StartThread_ThreadManager.restype = ctypes.c_int
        lib.StartThread_ThreadManager.argtypes = [ctypes.c_int]
            
        lib.TerminateThread_ThreadManager.restype = ctypes.c_int
        lib.TerminateThread_ThreadManager.argtypes = []
            
        lib.trigger_calc_SH_data_mt.restype = ctypes.c_int
        lib.trigger_calc_SH_data_mt.argtypes = []
            
        sizex = 128
        sizey = 128
        vImg = np.arange(0, sizex*sizey, 1, np.int32)
        nx = 11
        ny = 11
        v_xx = np.arange(0, nx, 1, np.int32)
        v_yy = np.arange(0, ny, 1, np.int32)
        SH_xc = np.arange(0, nx*ny, 1, np.float32)
        SH_yc = np.arange(0, nx*ny, 1, np.float32)

        nThreads = 12
        retval = lib.StartThread_ThreadManager(nThreads)
        print('retval : %d' %(retval))

        #for i in np.arange(10):
        retval = lib.trigger_calc_SH_data_mt()


        #print(SH_xc)

        #vImg = np.array(fits.getdata('c:\\home\\live_img.fits',0),dtype=np.int32).reshape(sizex*sizey,1)
        #v_xx = np.array(fits.getdata('c:\\home\\xx.fits',0),dtype=np.int32).reshape(nx,1)
        #v_yy = np.array(fits.getdata('c:\\home\\yy.fits',0),dtype=np.int32).reshape(ny,1)

        #print(v_xx)
        #print(v_yy)

        #retval = lib.calc_SH_data(sizex,sizey,vImg,nx,ny,v_xx,v_yy,SH_xc,SH_yc)

        clockA = time.perf_counter()
        #retval = lib.trigger_calc_SH_data(sizex,sizey,vImg,nx,ny,v_xx,v_yy,SH_xc,SH_yc)
        clockB = time.perf_counter()
        print(clockB - clockA)

        _ctypes.FreeLibrary(lib._handle)



        #clockA = time.perf_counter()
        #fib.array_test(100,myarray)
        #clockB = time.perf_counter()
        #print(clockB - clockA)