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


rlib.detector_response_equatorial_py.argtypes = \
    [ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double)
]
rlib.detector_response_equatorial_py.restype = ctypes.c_void_p

def detector_response_equatorial_py(detector, ra,dec, psi,gmst, active_polarizations):
    response_functions = (ctypes.c_double*6)()
    d_type = ctypes.c_double * (6)
    b_type = ctypes.c_bool *(6)
    
    rlib.detector_response_equatorial_py(detector.encode('utf-8'),ra, dec,psi,  gmst, b_type(*(np.ascontiguousarray(active_polarizations,dtype=ctypes.c_bool))), response_functions)
    return np.array(response_functions)


def DTOA_DETECTOR_py( RA, DEC, GMST_rad, detector1, detector2):
    returnval = rlib.DTOA_DETECTOR_py(
        RA,
        DEC,
        GMST_rad,
        detector1.encode('utf-8'),
        detector2.encode('utf-8'))
    return returnval
def populate_noise_py(frequencies,detector, integration_time):
    f=rlib.populate_noise_py
    f.argtypes=[\
        ctypes.POINTER(ctypes.c_double),\
        ctypes.c_char_p,ctypes.POINTER(ctypes.c_double),\
        ctypes.c_int,\
        ctypes.c_double \
    ]
    length = len(frequencies) 
    d_dtype = ctypes.c_double*(len(frequencies))
    detect = detector.encode("utf-8")
    noise_root = (ctypes.c_double*length)()
    f(d_dtype (*np.ascontiguousarray(frequencies,dtype=ctypes.c_double)), detect,noise_root,ctypes.c_int(len(frequencies)),ctypes.c_double(integration_time))
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
