#ifndef MCMC_ROUTINES_H
#define MCMC_ROUTINES_H
#include <complex>
#include <fftw3.h>
#include "util.h"
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
/*! \file 
 */

//########################################
//MCMC global variables - facilitate wrapping 
//of likelihood function 
//extern double **mcmc_noise=null;
//extern std::complex<double> **mcmc_data=null;
//extern double **mcmc_frequencies=null;
//extern std::string *mcmc_detectors=null;
//extern std::string *generation_method=null;
//extern int *mcmc_data_length = null;
//extern double **mcmc_noise;
//extern std::complex<double> **mcmc_data;
//extern double **mcmc_frequencies;
//extern std::string *mcmc_detectors;
//extern std::string *generation_method;
//extern int *mcmc_data_length ;
static double **mcmc_noise=NULL;
static std::complex<double> **mcmc_data=NULL;
static double **mcmc_frequencies=NULL;
static std::string *mcmc_detectors=NULL;
static std::string mcmc_generation_method="";
static int *mcmc_data_length=NULL ;
static fftw_outline *mcmc_fftw_plans=NULL ;
static int mcmc_num_detectors=2;
static double mcmc_gps_time=0;
static gsl_interp_accel **mcmc_accels = NULL;
static gsl_spline **mcmc_splines = NULL;
//extern const double **mcmc_noise=NULL;
//extern const std::complex<double> **mcmc_data=NULL;
//extern const double **mcmc_frequencies=NULL;
//extern const std::string *mcmc_detectors=NULL;
//extern const std::string *generation_method=NULL;
//extern const int *mcmc_data_length=NULL ;
//########################################


double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double SNR,
				double chirpmass,	
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan);
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
				bool NSflag);

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
				fftw_outline *plan);
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan);

double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag);

double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan);


double maximized_Log_Likelihood_aligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan
				);
double Log_Likelihood(std::complex<double> *data,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);

double maximized_Log_Likelihood_unaligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *hplus,
				std::complex<double> *hcross,
				size_t length,
				fftw_outline *plan
				);
double maximized_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan,
				double *tc,
				double *phic
				);
double maximized_coal_Log_Likelihood_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan,
				double *tc,
				double *phic
				);
double Log_Likelihood_internal(std::complex<double> *data,
			double *psd,
			double *frequencies,
			std::complex<double> *detector_response,
			int length,
			fftw_outline *plan
			);
void MCMC_MH_GW(double ***output,
			int dimension,
			int N_steps,
			int chain_N,
			double *initial_pos,
			double *chain_temps,
			int swp_freq,
			double(*log_prior)(double *param, int dimension, int chain_id),
			int numThreads,
			bool pool,
			bool show_prog,
			int num_detectors,
			std::complex<double> **data,
			double **noise_psd,
			double **frequencies,
			int *data_length,
			double gps_time,
			std::string *detector,
			std::string generation_method,
			std::string statistics_filename,
			std::string chain_filename,
			std::string auto_corr_filename
			);
void MCMC_fisher_wrapper(double *param, int dimension, double **output, int chain_id) ;
double MCMC_likelihood_wrapper(double *param, int dimension, int chain_id) ;
#endif
