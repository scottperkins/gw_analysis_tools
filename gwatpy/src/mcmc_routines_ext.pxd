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
cimport waveform_generator_ext
#from waveform_generator_ext import gen_params

cdef extern from "autocorrelation.h" :
    void write_auto_corr_file_from_data_file( 
                string autocorr_filename,
                string datafile,
                int length,
                int dimension,
                int num_segments,
                double target_corr,
                int num_threads)

    void write_auto_corr_file_from_data( 
                string autocorr_filename,
                double **data,
                int length,
                int dimension,
                int num_segments,
                double target_corr,
                int num_threads)

cdef extern from "mcmc_gw.h" :
    struct fftw_outline:
        pass
    #double Log_Likelihood(double *data_real,
    #    double *data_imag,
    #    double *psd,
    #    size_t length,
    #    gen_params *params,
    #    string detector,
    #    string generation_method,
    #    fftw_outline *plan) 
        
#    double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
#				size_t length,
#				double *real_data,
#				double *imag_data,
#				double *noise,
#				double SNR,
#				double chirpmass,	
#				double symmetric_mass_ratio, 
#				double spin1,
#				double spin2,
#				bool NSflag,
#                                fftw_outline *plan)
#    double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
#				size_t length,
#				double *real_data,
#				double *imag_data,
#				double *noise,
#				double chirpmass,	
#				double symmetric_mass_ratio, 
#				double spin1,
#				double spin2,
#                                double Luminosity_Distance,
#                                double theta,
#                                double phi,
#                                double iota,
#				bool NSflag,
#                                fftw_outline *plan)
#    double maximized_coal_Log_Likelihood(double *data_real, 
#				double *data_imag,
#				double *psd,
#				double *frequencies,
#				size_t length,
#				waveform_generator_ext.gen_params *params,
#				string detector,
#				string generation_method,
#				fftw_outline *plan
#				)
    void allocate_FFTW_mem_forward(fftw_outline *plan,int length)
    void deallocate_FFTW_mem(fftw_outline *plan) 
    
