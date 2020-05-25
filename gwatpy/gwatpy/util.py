import ctypes
import gwatpy.config as cf

rlib = ctypes.cdll.LoadLibrary(cf.LIB)

c = 299792458.
T_year = 31557600.
MPC_SEC = 3.085677581491367278913937957796471611e22/c
def calculate_chirpmass_py(mass1,mass2):
    f=rlib.calculate_chirpmass_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(mass1),ctypes.c_double(mass2),ctypes.byref(cm))
    return(cm.value)

def calculate_eta_py(mass1,mass2):
    f=rlib.calculate_eta_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(mass1),ctypes.c_double(mass2),ctypes.byref(cm))
    return(cm.value)

def calculate_mass1_py(chirpmass,eta):
    f=rlib.calculate_mass1_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(chirpmass),ctypes.c_double(eta),ctypes.byref(cm))
    return(cm.value)

def calculate_mass2_py(chirpmass,eta):
    f=rlib.calculate_mass2_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(chirpmass),ctypes.c_double(eta),ctypes.byref(cm))
    return(cm.value)

def DL_from_Z_py(z, COSMOLOGY):
    f=rlib.DL_from_Z_py
    f.argtypes=[\
        ctypes.c_double, \
        ctypes.c_char_p,\
        ctypes.POINTER(ctypes.c_double) \
    ]
    COSMO = COSMOLOGY.encode("utf-8")
    cm = ctypes.c_double()
    f(ctypes.c_double(z),COSMO,ctypes.byref(cm))
    return(cm.value)
