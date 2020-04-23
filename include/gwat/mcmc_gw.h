#ifndef MCMC_GW_H
#define MCMC_GW_H
#include <complex>
#include <fftw3.h>
#include "util.h"
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
/*! \file 
 *
 * Header file for the Graviational Wave specific MCMC routines
 */

//########################################
//MCMC global variables - facilitate wrapping 
//of likelihood function 
static double **mcmc_noise=NULL;
static std::complex<double> **mcmc_data=NULL;
static double **mcmc_frequencies=NULL;
static std::string *mcmc_detectors=NULL;
static std::string mcmc_generation_method="";
static int *mcmc_data_length=NULL ;
static fftw_outline *mcmc_fftw_plans=NULL ;
static int mcmc_num_detectors=2;
static double mcmc_gps_time=0;
static double mcmc_gmst=0;
static int mcmc_max_dim;
static int mcmc_min_dim;
static int mcmc_Nmod;
static int mcmc_Nmod_max;
static int *mcmc_bppe;
static gsl_interp_accel **mcmc_accels = NULL;
static gsl_spline **mcmc_splines = NULL;
static bool mcmc_log_beta;
static bool mcmc_intrinsic;
static bool mcmc_save_waveform;
static int mcmc_deriv_order=2;
static gsl_rng **mcmc_rvec;
//########################################

struct MCMC_modification_struct
{
	int ppE_Nmod = 0; //ppE
	int *bppe = NULL; //ppE
	int gIMR_Nmod_phi = 0; //gIMR
	int *gIMR_phii = NULL; //gIMR
	int gIMR_Nmod_sigma = 0; //gIMR
	int *gIMR_sigmai = NULL; //gIMR
	int gIMR_Nmod_beta = 0; //gIMR
	int *gIMR_betai = NULL; //gIMR
	int gIMR_Nmod_alpha = 0; //gIMR
	int *gIMR_alphai = NULL; //gIMR
};
static MCMC_modification_struct *mcmc_mod_struct;

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
				gen_params_base<double> *params,
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
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				//gen_params *params,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				);
double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				double *template_real,
				double *template_imag,
				fftw_outline *plan
				);
double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
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

double MCMC_likelihood_wrapper_SKYSEARCH(double *param, int dimension, int chain_id,void *parameters);
void SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	int *bppe,
	std::complex<double> *hplus,
	std::complex<double> *hcross,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void continue_RJPTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int max_dim,
	int min_dim,
	int N_steps,
	int swp_freq,
	double(*log_prior)(double *param, int *status, int dimension, int chain_id,void *parameters),
	void (*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dim, int chain_id, double step_width,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	int *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename
	);
void RJPTMCMC_MH_GW(double ***output,
	int ***status,
	int max_dim,
	int min_dim,
	int N_steps,
	int chain_N,
	double *initial_pos,
	int *initial_status,
	double *seeding_var,	
	double *chain_temps,
	int swp_freq,
	double(*log_prior)(double *param, int *status, int dimension, int chain_id,void *parameters),
	void (*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, int max_dim, int chain_id, double step_width,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod_max,
	int *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file);

void PTMCMC_MH_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
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
	int Nmod,
	int *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(std::string checkpoint_file_start,
			double **output,
			int N_steps,
			int max_chain_N_thermo_ensemble,
			double *chain_temps,
			int swp_freq,
			int t0,
			int nu,
			int corr_threshold,
			int corr_segments,
			double corr_converge_thresh,
			double corr_target_ac,
			std::string chain_distribution_scheme,
			double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
			int numThreads,
			bool pool,
			bool show_prog,
			int num_detectors,
			std::complex<double> **data,
			double **noise_psd,
			double **frequencies,
			int *data_length,
			double gps_time,
			std::string *detectors,
			int Nmod,
			int *bppe,
			std::string generation_method,
			std::string statistics_filename,
			std::string chain_filename,
			std::string likelihood_log_filename,
			std::string checkpoint_filename
			);
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void PTMCMC_MH_dynamic_PT_alloc_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	int *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	);
void continue_PTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int dimension,
	int N_steps,
	int swp_freq,
	double(*log_prior)(double *param, int dimension, int chain_id,void *parameters),
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
	int Nmod,
	int *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename,
	bool tune
	);
void PTMCMC_method_specific_prep(std::string generation_method, int dimension,double *seeding_var, bool local_seeding);

void RJPTMCMC_method_specific_prep(std::string generation_method, int max_dim, int min_dim,double *seeding_var, bool local_seeding);

double MCMC_likelihood_extrinsic(bool save_waveform, 
	gen_params_base<double> *parameters,
	std::string generation_method, 
	int *data_length, 
	double **frequencies, 
	std::complex<double> **data, 
	double **psd, 
	std::string *detectors, 
	fftw_outline *fftw_plans, 
	int num_detectors, 
	double RA, 
	double DEC,
	double gps_time);

void MCMC_fisher_wrapper(double *param, int dimension, double **output, int chain_id,void *parameters) ;

std::string MCMC_prep_params(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method);

double MCMC_likelihood_wrapper(double *param, int dimension, int chain_id,void *parameters) ;

double RJPTMCMC_likelihood_wrapper(double *param, int *status,int max_dim, int chain_id,void *parameters) ;

void RJPTMCMC_fisher_wrapper(double *param, int *status, int min_dim, double **output, int chain_id,void *parameters);
void RJPTMCMC_RJ_proposal(double *current_param, double *proposed_params, int *current_status, int *proposed_status,int max_dim, int chain_id, double step_width,void *parameters) ;
#endif
