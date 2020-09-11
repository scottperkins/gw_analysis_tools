import ctypes
import numpy as np
import gwatpy.config as cf

#rlib = ctypes.cdll.LoadLibrary(cf.LIB)
rlib = ctypes.CDLL(cf.LIB)

####################################################################
rlib.gen_params_base_py.argtypes = [ctypes.c_double,ctypes.c_double]
rlib.gen_params_base_py.restype = ctypes.c_void_p

rlib.gen_params_base_py_destructor.argtypes = [ctypes.c_void_p]
rlib.gen_params_base_py_destructor.restype = ctypes.c_void_p

rlib.fourier_waveform_py.argtypes=\
    [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), 
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
    ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p]
rlib.fourier_waveform_py.restype=ctypes.c_int
####################################################################

def waveform_generator(frequencies, generation_method, parameters):
    length = len(frequencies)
    array_type = ctypes.c_double*length
    f = np.ascontiguousarray(frequencies, dtype=ctypes.c_double)
    wf_plus_real = np.ascontiguousarray(np.zeros(length), dtype=ctypes.c_double)
    wf_plus_imag = np.ascontiguousarray(np.zeros(length), dtype=ctypes.c_double)
    wf_cross_real = np.ascontiguousarray(np.zeros(length), dtype=ctypes.c_double)
    wf_cross_imag = np.ascontiguousarray(np.zeros(length), dtype=ctypes.c_double)
    
    rlib.fourier_waveform_py(
        array_type(*f), 
        length, 
        array_type(*wf_plus_real), 
        array_type(*wf_plus_imag), 
        array_type(*wf_cross_real), 
        array_type(*wf_cross_imag), 
        generation_method.encode('utf-8'),
        parameters.obj)
    
    return wf_plus_real + 1j*wf_plus_imag, wf_cross_real + 1j*wf_cross_imag
class gen_params:
    mass1=10
    mass2=10 
    spin1_x = 0
    spin1_y = 0
    spin1_z = 0
    spin2_x = 0
    spin2_y = 0
    spin2_z = 0
    Luminosity_Distance = 100
    phiRef = 0
    tc = 3
    psi=0
    incl_angle=np.pi
    RA=0
    DEC=0
    gmst= 0
    f_ref= 20
    theta_l=0
    phi_l=0
    theta=0
    phi=0
    cosmology="PLANCK15"
    equatorial_orientation=False
    horizon_coord = False
    NSflag1=False
    NSflag2=False
    dep_postmerger=False
    shift_time=True
    shift_phase=True
    sky_average=False
    LISA_alpha0=0
    LISA_phi0=0
    Nmod_phi = 0
    Nmod_sigma = 0
    Nmod_beta = 0
    Nmod_alpha = 0
    phii = []
    sigmai = []
    betai = []
    alphai = []
    delta_phi= []
    delta_sigma= []
    delta_beta= []
    delta_alpha= []
    Nmod=0
    bppe=[]
    betappe=[]

    
    def __init__(self,**kwargs):
        
        if "mass1" in kwargs:
            self.mass1 = kwargs["mass1"] 
        if "mass2" in kwargs:
            self.mass2 = kwargs["mass2"] 
        if "spin1_x" in kwargs:
            self.spin1_x = kwargs["spin1_x"] 
        if "spin1_y" in kwargs:
            self.spin1_y = kwargs["spin1_y"] 
        if "spin1_z" in kwargs:
            self.spin1_z = kwargs["spin1_z"] 
        if "spin2_x" in kwargs:
            self.spin2_x = kwargs["spin2_x"] 
        if "spin2_y" in kwargs:
            self.spin2_y = kwargs["spin2_y"] 
        if "spin2_z" in kwargs:
            self.spin2_z = kwargs["spin2_z"] 
        self.obj = rlib.gen_params_base_py(self.mass1,self.mass2)
    def __del__(self):
        rlib.gen_params_base_py_destructor(self.obj)

