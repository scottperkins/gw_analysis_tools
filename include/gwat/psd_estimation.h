#ifndef PSD_ESTIMATION_H
#define PSD_ESTIMATION_H

#include "mcmc_sampler.h"

struct bayesline_sampling_struct
{
	int N_L_MIN = 0;	
	int N_L_MAX = 0;	
	int N_S_MIN = 0;	
	int N_S_MAX = 0;	
	std::complex<double> *data;
	double *data_mod_sq;
	double *frequencies;
	int data_length;
	double signal_length;
	double deltaf_factor = 50;
	double MIN_FREQ = 10;
	double MAX_FREQ = 4096;
	double MIN_SN_AMP = 1e-30;
	double MAX_SN_AMP = 1;
	double MIN_L_AMP = 1e-30;
	double MAX_L_AMP = 1;
	double MIN_Q = 1e-30;
	double MAX_Q = 1;
};

struct PSD_output
{
	double *SN;
	double T;
	int length;
};

void bayesline_psd_estimation(
	std::complex<double> *data_stream , 
	double *frequencies,
	int length,
	double T,
	int N_L_MIN, 
	int N_L_MAX, 
	int N_S_MIN, 
	int N_S_MAX,
	double *initial_pos,
	int *initial_status,
	double *seeding_var,
	std::string chain_allocation_scheme,
	int samples,	
	int chain_N,
	int max_chain_ensemble,
	double *chain_temps,
	int swap_freq,
	int t0,
	int nu,	
	int max_chunk_size,
	int threads,
	bool pool,
	bool show_prog,
	std::string chain_file,
	std::string stat_file,
	std::string checkpoint_file,
	PSD_output *output
	);

void bayesline_RJ_proposal(double *current_pos, double *prop_pos,int *current_status, int *prop_status,int *current_model, int *prop_model, double *MH_correction,mcmc_data_interface *interface, void *parameters);
double bayesline_prior(double *pos, int *status, int model, mcmc_data_interface *interface, void *parameters);
double bayesline_likelihood(double *pos, int *status, int model, mcmc_data_interface *interface, void *parameters);
void smooth_component(double *pos,int NS, bayesline_sampling_struct *p,double *frequencies, int length, double *SN );
void order_list(double *pos, int NS, double *pts, double *SNs);
void lorentzian_component(double *pos, int NL, bayesline_sampling_struct *p, double *frequencies, int length, double *SN);
double lorentzian(double pt, double amp, double q, double f,double deltaf);
#endif
