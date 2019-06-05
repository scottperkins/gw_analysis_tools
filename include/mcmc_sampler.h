#ifndef MCMC_SAMPLER_H
#define MCMC_SAMPLER_H
#include <iostream>
#include <functional>
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
		double *seeding_var,
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
void MCMC_MH(	double ***output, 
		int dimension, 	
		int N_steps,	
		int chain_N,	
		double *initial_pos, 	
		double *seeding_var,
		double *chain_temps,	
		int swp_freq,	
		double (*log_prior)(double *param, int dimension, int chain_id),	
		double (*log_likelihood)(double *param, int dimension, int chain_id),	
		void (*fisher)(double *param, int dimension, double **fisher, int chain_id),
		int numThreads,
		bool pool,
		bool show_prog,
		std::string statistics_filename,
		std::string chain_filename,
		std::string auto_corr_filename
		);
void MCMC_MH_internal(	double ***output, 
		int dimension, 	
		int N_steps,	
		int chain_N,	
		double *initial_pos, 	
		double *seeding_var,
		double *chain_temps,	
		int swp_freq,	
		//double (*log_prior)(double *param, int dimension, int chain_id),	
		//double (*log_likelihood)(double *param, int dimension, int chain_id),	
		//void (*fisher)(double *param, int dimension, double **fisher, int chain_id),
		std::function<double(double*,int,int)> log_prior,
		std::function<double(double*,int,int)> log_likelihood,
		std::function<void(double*,int,double**,int)>fisher,
		int numThreads,
		bool pool,
		bool show_prog,
		std::string statistics_filename,
		std::string chain_filename,
		std::string auto_corr_filename
		);
#endif
