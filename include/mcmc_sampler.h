#ifndef MCMC_SAMPLER_H
#define MCMC_SAMPLER_H
#include <iostream>
/* \file 
 * Header file for mcmc_sampler
 */

void mcmc_step_threaded(int j);
void mcmc_swap_threaded(int i, int j);

void MCMC_MH(	double ***output, 
		int dimension, 	
		int N_steps,	
		int chain_N,	
		double *initial_pos, 	
		double *chain_temps,	
		int swp_freq,	
		double (*log_prior)(double *param, int dimension),	
		double (*log_likelihood)(double *param, int dimension),	
		void (*fisher)(double *param, int dimension, double **fisher),
		int numThreads,
		bool pool,
		bool show_prog,
		std::string statistics_filename,
		std::string chain_filename,
		std::string auto_corr_filename
		);
#endif
