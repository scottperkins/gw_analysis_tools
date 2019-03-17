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
sys.path.append(scriptpath+"/src")
cimport mcmc_routines_ext
cimport waveform_generator_ext
#from waveform_generator_ext import gen_params_py


@cython.auto_pickle(True)
cdef class fftw_outline_py:
    cdef mcmc_routines_ext.fftw_outline plan
    cdef int N
    def __init__(self,N):
        self.N = N
        mcmc_routines_ext.initiate_likelihood_function(&self.plan,N)
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
    return mcmc_routines_ext.maximized_coal_log_likelihood_IMRPhenomD(&frequencies[0],
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
    return mcmc_routines_ext.maximized_coal_log_likelihood_IMRPhenomD_Full_Param(&frequencies[0],
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
				waveform_generator_ext.gen_params_py params,
				string detector,
				string generation_method,
				fftw_outline_py plan):
    #print("Inside pyx file")
    #print(params.params.mass1)
    #print(params.params.mass2)
    #print(params.params.Luminosity_Distance)
    #print(params.params.phic)
    #print(params.params.tc)
    #print(params.params.bppe)
    #print(params.params.betappe)
    #print(params.params.incl_angle)
    #print(params.params.theta)
    #print(params.params.phi)
    #cdef double m1 = params.params.mass1
    #cdef double m2 = params.params.mass2
    #cdef double symmetric_mass_ratio =(m1*m2)/ (m1+m2)**2.
    #cdef double chirpmass = (m1 + m2) * symmetric_mass_ratio**(3./5)
    #print(m1,m2)
    #return mcmc_routines_ext.maximized_coal_log_likelihood_IMRPhenomD_Full_Param(&frequencies[0],
    #    			frequencies.size,
    #    			&data_real[0],
    #    			&data_imag[0],
    #    			&psd[0],
    #    			chirpmass,	
    #    			symmetric_mass_ratio, 
    #    			params.params.spin1[2],
    #    			params.params.spin2[2],
    #                            params.params.Luminosity_Distance,
    #                            params.params.theta,
    #                            params.params.phi,
    #                            params.params.incl_angle,
    #    			params.params.NSflag,
    #                            &plan.plan)
    return mcmc_routines_ext.maximized_coal_Log_Likelihood(&data_real[0],
                                &data_imag[0],
                                &psd[0],
                                &frequencies[0],
                                frequencies.size,
                                &params.params,
                                detector,
                                generation_method,
                                &plan.plan)
def initiate_likelihood_function_py(fftw_outline_py plan,int length):
    mcmc_routines_ext.initiate_likelihood_function(&plan.plan,length)
def deactivate_likelihood_function_py(fftw_outline_py plan):
    mcmc_routines_ext.deactivate_likelihood_function(&plan.plan)
