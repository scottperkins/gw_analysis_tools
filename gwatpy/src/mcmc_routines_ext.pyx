## @file
#
#  File that wraps the code in mcmc_gw.cpp, mcmc_sampler.cpp, mcmc_sampler_internals.cpp, autocorrelation.cpp
from numpy cimport ndarray
import numpy as np
from libcpp cimport bool
from libcpp cimport complex
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.complex cimport complex
from libc.stdlib cimport malloc, free
cimport numpy as np
cimport cython
import os
#scriptpath = os.path.dirname(os.path.realpath(__file__)) 
import sys
#sys.path.append(scriptpath)
#sys.path.append(scriptpath+"/src")
cimport mcmc_routines_ext
#cimport waveform_generator_ext
#import waveform_generator_ext
#from waveform_generator_ext import gen_params_py


@cython.auto_pickle(True)
cdef class fftw_outline_py:
    cdef mcmc_routines_ext.fftw_outline plan
    cdef int N
    def __init__(self,N):
        self.N = N
        mcmc_routines_ext.allocate_FFTW_mem_forward(&self.plan,N)
    def __reduce__(self):
        return (self.__class__,(self.N,))

def write_auto_corr_file_from_data_file_py(
        string autocorr_filename,
        string datafile,
        int length,
        int dimension,
        int num_segments,
        double target_corr,
        int num_threads,
        bool cumulative):
    return write_auto_corr_file_from_data_file(autocorr_filename, datafile,length, dimension, num_segments, target_corr,num_threads,cumulative)

@cython.boundscheck(False)
@cython.wraparound(False)
def write_auto_corr_file_from_data_py(
        string autocorr_filename,
        double[:,::1] data,
        int length,
        int dimension,
        int num_segments,
        double target_corr,
        int num_threads,
        bool cumulative):
    #Not ideal -- have to wrap the memview in a real c++ array
    cdef double **temparr = <double **>malloc(sizeof(double*) * length)
    for i in np.arange(length):
        temparr[i] = &data[i,0]
    write_auto_corr_file_from_data(autocorr_filename, temparr,length, dimension, num_segments, target_corr,num_threads,cumulative)
    #for i in np.arange(length):
    #    free(temparr[i])
    free(temparr)

    
#def maximized_coal_log_likelihood_IMRPhenomD_py(double[::1] frequencies ,
#				double[::1] real_data ,
#				double[::1] imag_data ,
#				double[::1] noise ,
#				double SNR,
#				double chirpmass,	
#				double symmetric_mass_ratio, 
#				double spin1,
#				double spin2,
#				bool NSflag,
#                                fftw_outline_py plan):
#    return mcmc_routines_ext.maximized_coal_log_likelihood_IMRPhenomD(&frequencies[0],
#				frequencies.size,
#				&real_data[0],
#				&imag_data[0],
#				&noise[0],
#				SNR,
#				chirpmass,	
#				symmetric_mass_ratio, 
#				spin1,
#				spin2,
#				NSflag,
#                                &plan.plan )
#def maximized_coal_log_likelihood_IMRPhenomD_Full_Param_py(double[::1] frequencies ,
#				double[::1] real_data ,
#				double[::1] imag_data ,
#				double[::1] noise ,
#				double chirpmass,	
#				double symmetric_mass_ratio, 
#				double spin1,
#				double spin2,
#                                double Luminosity_Distance,
#                                double theta,
#                                double phi,
#                                double iota,
#				bool NSflag,
#                                fftw_outline_py plan):
#    return mcmc_routines_ext.maximized_coal_log_likelihood_IMRPhenomD_Full_Param(&frequencies[0],
#				frequencies.size,
#				&real_data[0],
#				&imag_data[0],
#				&noise[0],
#				chirpmass,	
#				symmetric_mass_ratio, 
#				spin1,
#				spin2,
#                                Luminosity_Distance,
#                                theta,
#                                phi,
#                                iota,
#				NSflag,
#                                &plan.plan)
#def maximized_coal_Log_Likelihood_py(double[::1] data_real, 
#				double[::1] data_imag,
#				double[::1] psd,
#				double[::1] frequencies,
#				waveform_generator_ext.gen_params_py params,
#				string detector,
#				string generation_method,
#				fftw_outline_py plan):
#    return mcmc_routines_ext.maximized_coal_Log_Likelihood(&data_real[0],
#                                &data_imag[0],
#                                &psd[0],
#                                &frequencies[0],
#                                frequencies.size,
#                                &params.params,
#                                detector,
#                                generation_method,
#                                &plan.plan)
def allocate_FFTW_mem_forward_py(fftw_outline_py plan,int length):
    mcmc_routines_ext.allocate_FFTW_mem_forward(&plan.plan,length)
def deallocate_FFTW_mem_py(fftw_outline_py plan):
    mcmc_routines_ext.deallocate_FFTW_mem(&plan.plan)
