import ctypes
import numpy as np
import gwatpy.config as cf

rlib = ctypes.cdll.LoadLibrary(cf.LIB)

###############################################
rlib.DTOA_DETECTOR_py.argtypes = \
    [ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_char_p,
    ctypes.c_char_p,
]
rlib.DTOA_DETECTOR_py.restype = ctypes.c_double
###############################################

def DTOA_DETECTOR_py( RA, DEC, GMST_rad, detector1, detector2):
    returnval = rlib.DTOA_DETECTOR_py(
        RA,
        DEC,
        GMST_rad,
        detector1.encode('utf-8'),
        detector2.encode('utf-8'))
    return returnval
def populate_noise_py(frequencies,detector, length, integration_time):
    f=rlib.populate_noise_py
    f.argtypes=[\
        ctypes.POINTER(ctypes.c_double),\
        ctypes.c_char_p,ctypes.POINTER(ctypes.c_double),\
        ctypes.c_int,\
        ctypes.c_double \
    ]
    
    freq = np.asarray(frequencies).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    detect = detector.encode("utf-8")
    noise_root = (ctypes.c_double*length)()
    f(freq, detect,noise_root,ctypes.c_int(length),ctypes.c_double(integration_time))
    return(noise_root)
def get_detector_parameters_py(detector):
    f = rlib.get_detector_parameters
    f.argtypes=[ctypes.c_char_p,ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(9*ctypes.c_double)]
    LAT = ctypes.c_double()
    LON = ctypes.c_double()
    LOC = (ctypes.c_double*3)()
    D = (ctypes.c_double*9)()
    f(detector.encode("utf-8"),ctypes.byref(LAT),ctypes.byref(LON),LOC,D )
    return LAT.value, LON.value, np.asarray(LOC), np.asarray(D).reshape(3,3)
