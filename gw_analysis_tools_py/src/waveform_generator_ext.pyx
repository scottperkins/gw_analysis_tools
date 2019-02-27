from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex


cdef extern from "waveform_generator.h" :
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

    int fourier_waveform(double *frequencies,
				int length,
                                double *waveform_real,
                                double *waveform_imag,
                                string generation_method,
                                gen_params *parameters)
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

    def __init__(self, mass1,mass2, DL, spin1, spin2, phic, tc, bppe, betappe, theta,phi,incl_angle, NSflag):
        self.params.mass1 = mass1
        self.params.mass2 = mass2
        self.params.Luminosity_Distance = DL
        self.params.spin1 = spin1
        self.params.spin2 = spin2
        self.params.phic = phic
        self.params.tc = tc
        self.params.bppe = bppe
        self.params.betappe = betappe
        self.params.incl_angle = incl_angle
        self.params.theta = theta
        self.params.phi = phi
        self.params.NSflag = NSflag
         
#cdef np.ndarray[np.complex128_t,ndim=1] fourier_waveform_py(double[::1] frequencies ,
def fourier_waveform_py(double[::1] frequencies ,
                                string generation_method,
                                gen_params_py parameters):
    cdef double[::1] waveform_real = np.zeros((frequencies.size),dtype=np.float64)
    cdef double[::1] waveform_imag= np.zeros((frequencies.size),dtype=np.float64)
    fourier_waveform(&frequencies[0],
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

def data_snr_maximized_extrinsic_py(double[::1] frequencies,
                                double[::1] data_real,
                                double[::1] data_imag,
                                string detector,
                                string generation_method,
                                gen_params_py parameters):
    cdef double snr = data_snr_maximized_extrinsic(&frequencies[0],
                                        frequencies.size,
                                        &data_real[0],
                                        &data_imag[0],
                                        detector,
                                        generation_method,
                                        parameters.params)
    return snr

