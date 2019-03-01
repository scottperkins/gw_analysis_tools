## @package mcmc_routines_ext
#Documentation
from numpy cimport ndarray
from libcpp cimport bool
from libcpp cimport complex
cimport numpy as np
cimport cython


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
#def maximized_coal_log_likelihood_IMRPhenomD_preinitialized_py(double SNR,
#				double chirpmass,	
#				double symmetric_mass_ratio, 
#				double spin1,
#				double spin2,
#				bool NSflag):
#    return maximized_coal_log_likelihood_IMRPhenomD(&freqs[0],
#				length,
#				&re_data[0],
#				&im_data[0],
#				&noise_c[0],
#				SNR,
#				chirpmass,	
#				symmetric_mass_ratio, 
#				spin1,
#				spin2,
#				NSflag)
def initiate_likelihood_function_py(fftw_outline_py plan,int length):
    initiate_likelihood_function(&plan.plan,length)
def deactivate_likelihood_function_py(fftw_outline_py plan):
    deactivate_likelihood_function(&plan.plan)
