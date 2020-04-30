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
#scriptpath = os.path.dirname(os.path.realpath(__file__)) 
import sys
#sys.path.append(scriptpath)
#cimport gwanalysistools_def 
cimport waveform_generator_ext

## Python wrapper for the generation parameters structure, as defined in util.cpp
#
cdef class gen_params_py:
    ## Constructor
    #  All quantities are detector frame
    #  @param mass1 Mass of the larger body (solar masses)
    #  @param mass2 Mass of the smaller body (solar masses)
    #  @param DL Luminosity Distance to the source (Mpc)
    #  @param spin1 Dimensionless spin of the larger body
    #  @param spin2 Dimensionless spin of the smaller body
    #  @param phic phase of the waveform at coalescence
    #  @param tc time of coalescence
    #  @param bppe ppE parameter b of all modifications
    #  @param betappe ppE parameter beta for all modifications
    #  @param Nmod Number of modifications to GR (length of bppe and betappe)
    #  @param theta azimuthal angle in spherical coordinates in the detector frame (rad)
    #  @param phi polar angle in spherical coordinates in the detector frame (rad)
    #  @param incl_angle Inclination angle of the binary at f_ref (rad)
    #  @param f_ref Reference frequency
    #  @param phiRef phase of the waveform at the reference frequency
    #  @param NSflag Boolean flag to indicate one or both components are neutron stars 
    def __init__(self, double mass1, double mass2, double DL, spin1, spin2,double phic,double tc,int[::1] bppe,double[::1] betappe,int Nmod,double theta,double phi,double incl_angle,double f_ref,double phiRef,bool NSflag):
        self.params.mass1 = mass1
        self.params.mass2 = mass2
        self.params.Luminosity_Distance = DL
        self.params.spin1 = spin1
        self.params.spin2 = spin2
        self.params.tc = tc
        self.params.bppe = &bppe[0]
        self.params.betappe = &betappe[0]
        self.params.Nmod = Nmod
        self.params.incl_angle = incl_angle
        self.params.theta = theta
        self.params.phi = phi
        self.params.f_ref = f_ref
        self.params.phiRef = phiRef
        self.params.NSflag1 = NSflag
        self.params.NSflag2 = NSflag
         
##Computes the waveform in Fourier space
# @param frequencies The array of frequencies to use
# @param generation_method Method to use for the waveform generation
# @param gen_params_py Parameters of the binary
def fourier_waveform_py(double[::1] frequencies ,
                                string generation_method,
                                gen_params_py parameters
                                ):
    cdef double[::1] waveform_real = np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    cdef double[::1] waveform_imag= np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    waveform_generator_ext.fourier_waveform(&frequencies[0],
				frequencies.size,
                                &waveform_real[0],
                                &waveform_imag[0],
                                generation_method,
                                &parameters.params)
    cdef np.ndarray[np.complex128_t,ndim=1] waveform = np.zeros((frequencies.size),dtype=np.complex128)
    cdef int i = 0
    while i<frequencies.size:
        waveform[i] = waveform_real[i] + 1j*waveform_imag[i]
        i = i +1
    return waveform

def fourier_amplitude_py(double[::1] frequencies ,
                                string generation_method,
                                gen_params_py parameters
                                ):
    cdef double[::1] amplitude = np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    waveform_generator_ext.fourier_amplitude(&frequencies[0],
				frequencies.size,
                                &amplitude[0],
                                generation_method,
                                &parameters.params)
    return amplitude
def fourier_phase_py(double[::1] frequencies ,
                                string generation_method,
                                gen_params_py parameters
                                ):
    cdef double[::1] phase = np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    waveform_generator_ext.fourier_phase(&frequencies[0],
				frequencies.size,
                                &phase[0],
                                generation_method,
                                &parameters.params)
    return phase
def fourier_waveform_polarizations_py(double[::1] frequencies ,
                                string generation_method,
                                gen_params_py parameters
                                ):
    cdef double[::1] waveform_plus_real = np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    cdef double[::1] waveform_plus_imag= np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    cdef double[::1] waveform_cross_real = np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    cdef double[::1] waveform_cross_imag= np.ascontiguousarray(np.zeros((frequencies.size),dtype=np.float64))
    waveform_generator_ext.fourier_waveform(&frequencies[0],
				frequencies.size,
                                &waveform_plus_real[0],
                                &waveform_plus_imag[0],
                                &waveform_cross_real[0],
                                &waveform_cross_imag[0],
                                generation_method,
                                &parameters.params)
    cdef np.ndarray[np.complex128_t,ndim=1] waveform_plus = np.zeros((frequencies.size),dtype=np.complex128)
    cdef np.ndarray[np.complex128_t,ndim=1] waveform_cross = np.zeros((frequencies.size),dtype=np.complex128)
    cdef int i = 0
    while i<frequencies.size:
        waveform_plus[i] = waveform_plus_real[i] + 1j*waveform_plus_imag[i]
        waveform_cross[i] = waveform_cross_real[i] + 1j*waveform_cross_imag[i]
        i = i +1
    return waveform_plus, waveform_cross
#def data_snr_maximized_extrinsic_py(double[::1] frequencies,
#                                double[::1] data_real,
#                                double[::1] data_imag,
#                                string detector,
#                                string generation_method,
#                                gen_params_py parameters):
#    cdef double snr = waveform_generator_ext.data_snr_maximized_extrinsic(&frequencies[0],
#                                        frequencies.size,
#                                        &data_real[0],
#                                        &data_imag[0],
#                                        detector,
#                                        generation_method,
#                                        parameters.params)
#    return snr

