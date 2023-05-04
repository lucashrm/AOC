from sys import version_info as _alpao_wrapper_version_info

# Since SWIG cannot generate wrapper for python3.lib

if _alpao_wrapper_version_info < (3, 8, 0): 
    from Lib.asdk37 import *
elif _alpao_wrapper_version_info < (3, 9, 0): 
    from Lib.asdk38 import *
else:
    print('Not supported Python version.');


