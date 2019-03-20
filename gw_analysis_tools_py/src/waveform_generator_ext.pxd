from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex
import os
scriptpath = os.path.dirname(os.path.realpath(__file__)) 
import sys
sys.path.append(scriptpath)

cdef extern from "util.h" :
    struct gen_params:
        double mass1
        double mass2
        double Luminosity_Distance
        double spin1[3]
        double spin2[3]
        double phic
        double tc
        int bppe
        double betappe
        double incl_angle
        double theta
        double phi
        bool NSflag
        double f_ref
        double phiRef

cdef extern from "waveform_generator.h" :
    int fourier_waveform(double *frequencies,
				int length,
                                double *waveform_real,
                                double *waveform_imag,
                                string generation_method,
                                gen_params *parameters)
    int fourier_waveform(double *frequencies,
				int length,
                                double *waveform_plus_real,
                                double *waveform_plus_imag,
                                double *waveform_cross_real,
                                double *waveform_cross_imag,
                                string generation_method,
                                gen_params *parameters)

    int fourier_amplitude(double *frequencies, 
    			int length,
    			double *amplitude, 
    			string generation_method,
    			gen_params *parameters
    			)
    
    int fourier_phase(double *frequencies, 
    			int length,
    			double *phase, 
    			string generation_method,
    			gen_params *parameters
    			)
cdef extern from "waveform_util.h" :
    double data_snr_maximized_extrinsic(double *frequencies,
                                int length,
                                double *data_real,
                                double *data_imag,
                                string detector,
                                string generation_method,
                                gen_params param
                                )

cdef class gen_params_py:
    cdef gen_params params

