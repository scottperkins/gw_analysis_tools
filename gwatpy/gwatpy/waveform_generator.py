import ctypes
import numpy as np
import gwatpy.config as cf
import gwatpy.util as gpu

#rlib = ctypes.cdll.LoadLibrary(cf.LIB)
rlib = ctypes.CDLL(cf.LIB)

####################################################################

rlib.time_waveform_full_py.argtypes=\
    [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), 
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p]
rlib.time_waveform_full_py.restype=ctypes.c_int

rlib.fourier_waveform_py.argtypes=\
    [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), 
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p]
rlib.fourier_waveform_py.restype=ctypes.c_int

rlib.time_detector_response_py.argtypes=\
    [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), 
    ctypes.POINTER(ctypes.c_double), \
    ctypes.c_char_p, ctypes.c_char_p,ctypes.c_void_p]
rlib.time_detector_response_py.restype=ctypes.c_int

rlib.fourier_detector_response_py.argtypes=\
    [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), 
    ctypes.POINTER(ctypes.c_double), \
    ctypes.c_char_p, ctypes.c_char_p,ctypes.c_void_p]
rlib.fourier_detector_response_py.restype=ctypes.c_int
####################################################################

def time_waveform_generator(times, generation_method, parameters):
    length = len(times)
    array_type = ctypes.c_double*length
    t = np.ascontiguousarray(times, dtype=ctypes.c_double)
    wf_plus_real = (ctypes.c_double*length)()
    wf_plus_imag = (ctypes.c_double*length)()
    wf_cross_real = (ctypes.c_double*length)()
    wf_cross_imag = (ctypes.c_double*length)()
    wf_x_real = (ctypes.c_double*length)()
    wf_x_imag = (ctypes.c_double*length)()
    wf_y_real = (ctypes.c_double*length)()
    wf_y_imag = (ctypes.c_double*length)()
    wf_b_real = (ctypes.c_double*length)()
    wf_b_imag = (ctypes.c_double*length)()
    wf_l_real = (ctypes.c_double*length)()
    wf_l_imag = (ctypes.c_double*length)()
    
    rlib.time_waveform_full_py(
        array_type(*t), 
        length, 
        wf_plus_real, 
        wf_plus_imag, 
        wf_cross_real, 
        wf_cross_imag, 
        wf_x_real, 
        wf_x_imag, 
        wf_y_real, 
        wf_y_imag, 
        wf_b_real, 
        wf_b_imag, 
        wf_l_real, 
        wf_l_imag, 
        generation_method.encode('utf-8'),
        parameters.obj)
    wf_plus_real = np.asarray(wf_plus_real)
    wf_plus_imag = np.asarray(wf_plus_imag)
    wf_cross_real = np.asarray(wf_cross_real)
    wf_cross_imag = np.asarray(wf_cross_imag)
    wf_x_real = np.asarray(wf_x_real)
    wf_x_imag = np.asarray(wf_x_imag)
    wf_y_real = np.asarray(wf_y_real)
    wf_y_imag = np.asarray(wf_y_imag)
    wf_b_real = np.asarray(wf_b_real)
    wf_b_imag = np.asarray(wf_b_imag)
    wf_l_real = np.asarray(wf_l_real)
    wf_l_imag = np.asarray(wf_l_imag)
    
    return wf_plus_real + 1j*wf_plus_imag, wf_cross_real + 1j*wf_cross_imag, \
            wf_x_real + 1j * wf_x_imag, wf_y_real + 1j * wf_y_imag, \
            wf_b_real + 1j * wf_b_imag, wf_l_real + 1j * wf_l_imag, 


def waveform_generator(frequencies, generation_method, parameters):
    length = len(frequencies)
    array_type = ctypes.c_double*length
    f = np.ascontiguousarray(frequencies, dtype=ctypes.c_double)
    wf_plus_real = (ctypes.c_double*length)()
    wf_plus_imag = (ctypes.c_double*length)()
    wf_cross_real = (ctypes.c_double*length)()
    wf_cross_imag = (ctypes.c_double*length)()
    
    rlib.fourier_waveform_py(
        array_type(*f), 
        length, 
        wf_plus_real, 
        wf_plus_imag, 
        wf_cross_real, 
        wf_cross_imag, 
        generation_method.encode('utf-8'),
        parameters.obj)
    wf_plus_real = np.asarray(wf_plus_real)
    wf_plus_imag = np.asarray(wf_plus_imag)
    wf_cross_real = np.asarray(wf_cross_real)
    wf_cross_imag = np.asarray(wf_cross_imag)
    
    return wf_plus_real + 1j*wf_plus_imag, wf_cross_real + 1j*wf_cross_imag



def time_response_generator(times, detector,generation_method, parameters):
    length = len(times)
    array_type = ctypes.c_double*length
    t = np.ascontiguousarray(times, dtype=ctypes.c_double)
    response_real = (ctypes.c_double*length)()
    response_imag = (ctypes.c_double*length)()
    
    rlib.time_detector_response_py(
        array_type(*t), 
        length, 
        response_real, 
        response_imag, 
        detector.encode('utf-8'),
        generation_method.encode('utf-8'),
        parameters.obj)
    response_real = np.asarray(response_real)
    response_imag = np.asarray(response_imag)
    
    return response_real + 1j*response_imag

def response_generator(frequencies, detector,generation_method, parameters):
    length = len(frequencies)
    array_type = ctypes.c_double*length
    f = np.ascontiguousarray(frequencies, dtype=ctypes.c_double)
    response_real = (ctypes.c_double*length)()
    response_imag = (ctypes.c_double*length)()
    
    rlib.fourier_detector_response_py(
        array_type(*f), 
        length, 
        response_real, 
        response_imag, 
        detector.encode('utf-8'),
        generation_method.encode('utf-8'),
        parameters.obj)
    response_real = np.asarray(response_real)
    response_imag = np.asarray(response_imag)
    
    return response_real + 1j*response_imag

def waveform_generator(frequencies, generation_method, parameters):
    length = len(frequencies)
    array_type = ctypes.c_double*length
    f = np.ascontiguousarray(frequencies, dtype=ctypes.c_double)
    wf_plus_real = (ctypes.c_double*length)()
    wf_plus_imag = (ctypes.c_double*length)()
    wf_cross_real = (ctypes.c_double*length)()
    wf_cross_imag = (ctypes.c_double*length)()
    
    rlib.fourier_waveform_py(
        array_type(*f), 
        length, 
        wf_plus_real, 
        wf_plus_imag, 
        wf_cross_real, 
        wf_cross_imag, 
        generation_method.encode('utf-8'),
        parameters.obj)
    wf_plus_real = np.asarray(wf_plus_real)
    wf_plus_imag = np.asarray(wf_plus_imag)
    wf_cross_real = np.asarray(wf_cross_real)
    wf_cross_imag = np.asarray(wf_cross_imag)
    
    return wf_plus_real + 1j*wf_plus_imag, wf_cross_real + 1j*wf_cross_imag
