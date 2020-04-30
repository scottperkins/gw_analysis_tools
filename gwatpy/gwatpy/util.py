import ctypes

rlib = ctypes.cdll.LoadLibrary("../build/install/lib/libgwat.so")

def calculate_chirpmass_py(mass1,mass2):
    f=rlib.calculate_chirpmass_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(29),ctypes.c_double(20),ctypes.byref(cm))
    return(cm.value)
