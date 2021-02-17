import numpy as np
import ctypes
import gwatpy.config as cf

rlib = ctypes.cdll.LoadLibrary(cf.LIB)
############################################################

rlib.gen_params_base_py.argtypes = \
    [ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_char_p,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_bool,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
]
rlib.gen_params_base_py.restype = ctypes.c_void_p


rlib.t_0PN_py.argtypes = \
    [ctypes.c_double,
    ctypes.c_double]
rlib.t_0PN_py.restype = ctypes.c_double

rlib.f_0PN_py.argtypes = \
    [ctypes.c_double,
    ctypes.c_double]
rlib.f_0PN_py.restype = ctypes.c_double

rlib.gen_params_base_py_destructor.argtypes = [ctypes.c_void_p]
rlib.gen_params_base_py_destructor.restype = ctypes.c_void_p

############################################################

c = 299792458.
T_year = 31557600.
MPC_SEC = 3.085677581491367278913937957796471611e22/c
MSOL_SEC = 4.925491025543575903411922162094833998e-6
def calculate_chirpmass_py(mass1,mass2):
    f=rlib.calculate_chirpmass_py
    f.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)]
    cm = ctypes.c_double()
    f(ctypes.c_double(mass1),ctypes.c_double(mass2),ctypes.byref(cm))
    return(cm.value)

def t_0PN_py(f ,chirpmass):
    return rlib.t_0PN_py(f,chirpmass)
def f_0PN_py(t ,chirpmass):
    return rlib.f_0PN_py(t,chirpmass)

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

class gen_params:
    mass1=10
    mass2=10 
    spin1 = [0,0,0]
    spin2 = [0,0,0]
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
    x0 = 0
    obj=None

    
    def __init__(self,**kwargs):
        
        if "mass1" in kwargs:
            self.mass1 = kwargs["mass1"] 
        if "mass2" in kwargs:
            self.mass2 = kwargs["mass2"] 
        if "spin1" in kwargs:
            self.spin1 = kwargs["spin1"] 
        if "spin2" in kwargs:
            self.spin2 = kwargs["spin2"] 
        if "Luminosity_Distance" in kwargs:
            self.Luminosity_Distance = kwargs["Luminosity_Distance"] 
        if "phiRef" in kwargs:
            self.phiRef = kwargs["phiRef"] 
        if "tc" in kwargs:
            self.tc = kwargs["tc"] 
        if "psi" in kwargs:
            self.psi = kwargs["psi"] 
        if "incl_angle" in kwargs:
            self.incl_angle = kwargs["incl_angle"] 
        if "RA" in kwargs:
            self.RA = kwargs["RA"] 
        if "DEC" in kwargs:
            self.DEC = kwargs["DEC"] 
        if "gmst" in kwargs:
            self.gmst = kwargs["gmst"] 
        if "f_ref" in kwargs:
            self.f_ref = kwargs["f_ref"] 
        if "theta_l" in kwargs:
            self.theta_l = kwargs["theta_l"] 
        if "phi_l" in kwargs:
            self.phi_l = kwargs["phi_l"] 
        if "theta" in kwargs:
            self.theta = kwargs["theta"] 
        if "phi" in kwargs:
            self.phi = kwargs["phi"] 
        if "cosmology" in kwargs:
            self.cosmology = kwargs["cosmology"] 
        if "equatorial_orientation" in kwargs:
            self.cosmology = kwargs["equatorial_orientation"] 
        if "horizon_coord" in kwargs:
            self.horizon_coord = kwargs["horizon_coord"] 
        if "NSflag1" in kwargs:
            self.NSflag1 = kwargs["NSflag1"] 
        if "NSflag2" in kwargs:
            self.NSflag2 = kwargs["NSflag2"] 
        if "dep_postmerger" in kwargs:
            self.dep_postmerger = kwargs["dep_postmerger"] 
        if "shift_time" in kwargs:
            self.shift_time = kwargs["shift_time"] 
        if "shift_phase" in kwargs:
            self.shift_phase = kwargs["shift_phase"] 
        if "sky_average" in kwargs:
            self.sky_average = kwargs["sky_average"] 
        if "LISA_alpha0" in kwargs:
            self.LISA_alpha0 = kwargs["LISA_alpha0"] 
        if "LISA_phi0" in kwargs:
            self.LISA_phi0 = kwargs["LISA_phi0"] 
        if "Nmod_phi" in kwargs:
            self.Nmod_phi = kwargs["Nmod_phi"] 
        if "Nmod_sigma" in kwargs:
            self.Nmod_sigma = kwargs["Nmod_sigma"] 
        if "Nmod_beta" in kwargs:
            self.Nmod_beta = kwargs["Nmod_beta"] 
        if "Nmod_alpha" in kwargs:
            self.Nmod_alpha = kwargs["Nmod_alpha"] 
        if "phii" in kwargs:
            self.phii = kwargs["phii"] 
        if "sigmai" in kwargs:
            self.sigmai = kwargs["sigmai"] 
        if "betai" in kwargs:
            self.betai = kwargs["betai"] 
        if "alphai" in kwargs:
            self.alphai = kwargs["alphai"] 
        if "delta_phi" in kwargs:
            self.delta_phi = kwargs["delta_phi"] 
        if "delta_sigma" in kwargs:
            self.delta_sigma = kwargs["delta_sigma"] 
        if "delta_beta" in kwargs:
            self.delta_beta = kwargs["delta_beta"] 
        if "delta_alpha" in kwargs:
            self.delta_alpha = kwargs["delta_alpha"] 
        if "Nmod" in kwargs:
            self.Nmod = kwargs["Nmod"] 
        if "bppe" in kwargs:
            self.bppe = kwargs["bppe"] 
        if "x0" in kwargs:
            self.x0 = kwargs["x0"] 
        if "betappe" in kwargs:
            self.betappe = kwargs["betappe"] 
        spin_array_type = 3*ctypes.c_double
        phii_array_type = self.Nmod_phi*ctypes.c_int
        sigmai_array_type = self.Nmod_sigma*ctypes.c_int
        betai_array_type = self.Nmod_beta*ctypes.c_int
        alphai_array_type = self.Nmod_alpha*ctypes.c_int

        delta_phi_array_type = self.Nmod_phi*ctypes.c_double
        delta_sigma_array_type = self.Nmod_sigma*ctypes.c_double
        delta_beta_array_type = self.Nmod_beta*ctypes.c_double
        delta_alpha_array_type = self.Nmod_alpha*ctypes.c_double

        bppe_type = self.Nmod * ctypes.c_double;
        betappe_type = self.Nmod * ctypes.c_double;

        self.obj = rlib.gen_params_base_py(
            self.mass1,
            self.mass2,
            spin_array_type(*(np.ascontiguousarray(self.spin1,dtype=ctypes.c_double))),
            spin_array_type(*(np.ascontiguousarray(self.spin2,dtype=ctypes.c_double))),
            self.Luminosity_Distance,
            self.incl_angle,
            self.RA,
            self.DEC,
            self.psi,
            self.gmst,
            self.f_ref,
            self.theta_l,
            self.phi_l,
            self.theta,
            self.phi,
            self.cosmology.encode('utf-8'),
            self.equatorial_orientation,
            self.horizon_coord,
            self.NSflag1,
            self.NSflag2,
            self.dep_postmerger,
            self.shift_time,
            self.shift_phase,
            self.sky_average,
            self.LISA_alpha0,
            self.LISA_phi0,
            self.Nmod_phi,
            self.Nmod_sigma,
            self.Nmod_beta,
            self.Nmod_alpha,
            phii_array_type(*(np.ascontiguousarray(self.phii,dtype=ctypes.c_int))),
            sigmai_array_type(*(np.ascontiguousarray(self.sigmai,dtype=ctypes.c_int))),
            betai_array_type(*(np.ascontiguousarray(self.betai,dtype=ctypes.c_int))),
            alphai_array_type(*(np.ascontiguousarray(self.alphai,dtype=ctypes.c_int))),
            delta_phi_array_type(*(np.ascontiguousarray(self.delta_phi,dtype=ctypes.c_double))),
            delta_sigma_array_type(*(np.ascontiguousarray(self.delta_sigma,dtype=ctypes.c_double))),
            delta_beta_array_type(*(np.ascontiguousarray(self.delta_beta,dtype=ctypes.c_double))),
            delta_alpha_array_type(*(np.ascontiguousarray(self.delta_alpha,dtype=ctypes.c_double))),
            self.Nmod,
            bppe_type(*(np.ascontiguousarray(self.bppe,dtype=ctypes.c_double))),
            betappe_type(*(np.ascontiguousarray(self.betappe,dtype=ctypes.c_double)))
        )


    def __del__(self):
        rlib.gen_params_base_py_destructor(self.obj)

