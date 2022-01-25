#ifndef MCMC_GW_EXTENDED_H
#define MCMC_GW_EXTENDED_H

#include "util.h"
#include "mcmc_gw.h"

struct mcmcVariables 
{
	double **mcmc_noise=nullptr;
	double *mcmc_init_pos=nullptr;
	std::complex<double> **mcmc_data=nullptr;
	double **mcmc_frequencies=nullptr;
	std::string *mcmc_detectors=nullptr;
	std::string mcmc_generation_method="";
	int *mcmc_data_length = nullptr;
	fftw_outline *mcmc_fftw_plans=nullptr;
	int mcmc_num_detectors  ;
	double mcmc_gps_time = 0;
	double mcmc_gmst = 0;
	int mcmc_max_dim;
	int mcmc_min_dim;
	MCMC_modification_struct *mcmc_mod_struct=nullptr;
	bool mcmc_save_waveform;	
	bool mcmc_intrinsic=false;	
	int mcmc_deriv_order = 4;
	bool mcmc_log_beta = false;
};

#endif
