from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex
import os
#scriptpath = os.path.dirname(os.path.realpath(__file__)) 
import sys
from libcpp.string cimport string
#sys.path.append(scriptpath)

cdef extern from "util.h" :
    double c,G, MSOL_SEC,MPC_SEC,T_year, AXIAL_TILT, AU_SEC
    double DL_from_Z(double Z, string cosmology)
    double calculate_mass1(double chirpmass, double eta)
    double calculate_mass2(double chirpmass, double eta)
cdef extern from "GWATConfig.h" :
    cdef string GWAT_ROOT_DIRECTORY


