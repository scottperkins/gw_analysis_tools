##@file 
#
#File that contains cython code to wrap the c++ library

## @package waveform_generator_ext
#
#Python wrapper for the waveform generation in waveform_generator.cpp


from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex
import os
import sys
cimport util

cpy = c
Gpy = G
MPC_SECpy = MPC_SEC
MSOL_SECpy = MSOL_SEC
T_yearpy = T_year
GWAT_ROOT_DIRECTORY_PY = GWAT_ROOT_DIRECTORY 
def DL_from_Z_py(double Z, string cosmology):
    return DL_from_Z(Z,cosmology)
def calculate_mass1_py(double chirpmass, double eta):
    return calculate_mass1(chirpmass,eta)
def calculate_mass2_py(double chirpmass, double eta):
    return calculate_mass2(chirpmass,eta)
