from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex

cdef extern from "mcmc_routines.h" :
    struct fftw_outline:
        pass
    double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
                                fftw_outline *plan)
    double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
                                double Luminosity_Distance,
                                double theta,
                                double phi,
                                double iota,
				bool NSflag,
                                fftw_outline *plan)
    double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				int length,
				gen_params *params,
				string detector,
				string generation_method,
				fftw_outline *plan
				)
    void initiate_likelihood_function(fftw_outline *plan,int length)
    void deactivate_likelihood_function(fftw_outline *plan) 
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
    int fourier_waveform(double *frequencies,
				int length,
                                double *waveform_plus_real,
                                double *waveform_plus_imag,
                                double *waveform_cross_real,
                                double *waveform_cross_imag,
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
