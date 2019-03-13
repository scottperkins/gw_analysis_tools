from numpy cimport ndarray
import numpy as np
from libcpp cimport bool
from libcpp cimport complex
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex
cimport numpy as np
cimport cython
import os
scriptpath = os.path.dirname(os.path.realpath(__file__)) 
import sys
sys.path.append(scriptpath)
cimport gwanalysistools_def 
#from waveform_generator_ext import gen_params_py

cdef class gen_params_py:
    cdef gwanalysistools_def.gen_params params

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
				gwanalysistools_def.gen_params *params,
				string detector,
				string generation_method,
				fftw_outline *plan
				)
    void initiate_likelihood_function(fftw_outline *plan,int length)
    void deactivate_likelihood_function(fftw_outline *plan) 

@cython.auto_pickle(True)
cdef class fftw_outline_py:
    cdef fftw_outline plan
    cdef int N
    def __init__(self,N):
        self.N = N
        initiate_likelihood_function(&self.plan,N)
    def __reduce__(self):
        return (self.__class__,(self.N,))
## Documentation
# hello
def maximized_coal_log_likelihood_IMRPhenomD_py(double[::1] frequencies ,
				double[::1] real_data ,
				double[::1] imag_data ,
				double[::1] noise ,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
                                fftw_outline_py plan):
    return maximized_coal_log_likelihood_IMRPhenomD(&frequencies[0],
				frequencies.size,
				&real_data[0],
				&imag_data[0],
				&noise[0],
				SNR,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				NSflag,
                                &plan.plan )
def maximized_coal_log_likelihood_IMRPhenomD_Full_Param_py(double[::1] frequencies ,
				double[::1] real_data ,
				double[::1] imag_data ,
				double[::1] noise ,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
                                double Luminosity_Distance,
                                double theta,
                                double phi,
                                double iota,
				bool NSflag,
                                fftw_outline_py plan):
    return maximized_coal_log_likelihood_IMRPhenomD_Full_Param(&frequencies[0],
				frequencies.size,
				&real_data[0],
				&imag_data[0],
				&noise[0],
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
                                Luminosity_Distance,
                                theta,
                                phi,
                                iota,
				NSflag,
                                &plan.plan)
def maximized_coal_Log_Likelihood_py(double[::1] data_real, 
				double[::1] data_imag,
				double[::1] psd,
				double[::1] frequencies,
				gen_params_py params,
				string detector,
				string generation_method,
				fftw_outline_py plan):
    return maximized_coal_Log_Likelihood(&data_real[0],
                                &data_imag[0],
                                &psd[0],
                                &frequencies[0],
                                frequencies.size,
                                &params.params,
                                detector,
                                generation_method,
                                &plan.plan)
def initiate_likelihood_function_py(fftw_outline_py plan,int length):
    initiate_likelihood_function(&plan.plan,length)
def deactivate_likelihood_function_py(fftw_outline_py plan):
    deactivate_likelihood_function(&plan.plan)
