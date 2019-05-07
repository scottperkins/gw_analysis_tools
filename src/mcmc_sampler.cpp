#include "mcmc_sampler.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "util.h"
#include "mcmc_sampler_internals.h"
#include <omp.h>

#ifndef _OPENMP
#define omp ignore
#endif

#include <unistd.h>

/*!\file 
 *
 * Source file for generic MCMC sampler. Sub routines that are 
 * application agnostic are housed in mcmc_sampler_internals
 */

//Random number variables
const gsl_rng_type *T;
gsl_rng * r;

void testmcmc()
{
	int x = 0, n = 1000;
	for (x =0;x<n;x++){
		usleep(500);
		printProgress((double)x/n);	
	}
	std::cout<<std::endl;
	
}


/*\brief Generic sampler, where the likelihood, prior are parameters supplied by the user
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
 */
void MCMC_MH(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
		int dimension, 	/**< dimension of the parameter space being explored*/
		int N_steps,	/**< Number of total steps to be taken, per chain*/
		int chain_N,	/**< Number of chains*/
		double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
		double *chain_temps,	/**<Double array of temperatures for the chains*/
		int swp_freq,	/**< the frequency with which chains are swapped*/
		double (*log_prior)(double *param, int dimension),	/**<Funcion pointer for the log_prior*/
		double (*log_likelihood)(double *param, int dimension),	/**<Function pointer for the log_likelihood*/
		void (*fisher)(double *param, int dimension, double **fisher),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
		std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
		std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
		std::string auto_corr_filename/**< Filename to output auto correlation in some interval, if empty string, not output*/
		)
{
	omp_set_num_threads(15);
	//random number generator initialization
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r = gsl_rng_alloc(T);

	//Array holding the probability of each 
	//type of step - Gaussian, differential evolution, MMALA, Fisher
	sampler sampler;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		sampler.fisher_exist = false;
	}
	else 
		sampler.fisher_exist = true;
	
	
	//Construct sampler structure
	sampler.lp = log_prior;
	sampler.ll = log_likelihood;
	sampler.fish = fisher;
	sampler.swp_freq = swp_freq;
	sampler.chain_temps = chain_temps;
	sampler.chain_N = chain_N;
	sampler.dimension = dimension;
	sampler.r = r;
	sampler.history_length = 100;
	
	allocate_sampler_mem(&sampler);
	for (int chain_index; chain_index<sampler.chain_N; chain_index++)
		assign_probabilities(&sampler, chain_index);
	

	//int i,j, k=0;
	//int j, k=0;
	int  k=0;
	int swp_accepted=0, swp_rejected=0;
	int step_accepted=0, step_rejected=0;
	
	//Assign initial position to start chains
	//Currently, just set all chains to same initial position
	for (int j=0;j<sampler.chain_N;j++){
		for (int i = 0; i<sampler.dimension; i++)
		{
			sampler.de_primed[j]=false;
			output[j][0][i] = initial_pos[i];
		}
	}
	
	//Sampler Loop
	//#pragma omp parallel //num_threads(4)
	{
	while (k<N_steps-1){
		int cutoff ;
		if( N_steps-k <= sampler.swp_freq) cutoff = N_steps-k-1;	
		else cutoff = sampler.swp_freq;	
		#pragma omp parallel for reduction(+:step_accepted) reduction(+:step_rejected)
		for (int j=0; j<chain_N; j++)
		{
			for (int i = 0 ; i< cutoff;i++)
			{
				int success;
				success = mcmc_step(&sampler, output[j][k+i], output[j][k+i+1],j);	
				if(success==1){step_accepted+=1;}
				else{step_rejected+=1;}
				update_history(&sampler,output[j][k+i+1], j);
			}
			if(!sampler.de_primed[j]) 
			{
				if ((k+cutoff)>sampler.history_length)
				{
					sampler.de_primed[j]=true;
					assign_probabilities(&sampler,j);	
				}
			}
		}
		k+= cutoff;
		chain_swap(&sampler, output, k, &swp_accepted, &swp_rejected);
		printProgress((double)k/N_steps);	
	}
	}
	std::cout<<std::endl;
	double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	accepted_percent = (double)(step_accepted)/(step_accepted+step_rejected);
	rejected_percent = (double)(step_rejected)/(step_accepted+step_rejected);
	std::cout<<"Accepted percentage of steps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of steps (all chains): "<<rejected_percent<<std::endl;
	std::cout<<"NANS (all chains): "<<sampler.nan_counter<<std::endl;
	
	//###########################################################
	//Auto-correlation
	double **temp = (double **) malloc(sizeof(double*)*N_steps);
	for (int i = 0 ; i< sampler.dimension; i++){
		temp[i] = (double *)malloc(sizeof(double)*N_steps);
		for(int j =0; j< N_steps; j++){
			temp[i][j] = output[0][j][i];
		}
	}
	int segments = 10;
	double stepsize = (double)N_steps/segments;
	double lengths[segments];
	double ac[sampler.dimension][segments];
	for (int i =0; i<segments;i++)
		lengths[i]=stepsize*(1. + i);
	
	#pragma omp parallel for 
	for (int i = 0 ; i< sampler.dimension; i++)
	{
		for(int l =0; l<segments; l++){
			for(int j =0; j< lengths[l]; j++){
				temp[i][j] = output[0][j][i];
			}
			ac[i][l]=auto_correlation(temp[i],lengths[l], 0.01);
			//std::cout<<"autocorrelation: "<<ac[i][l]<<std::endl;
		}
	}	
	
	for (int i = 0 ; i< sampler.dimension; i++)
	{
		for(int l =0; l<segments; l++){
			std::cout<<"AUTO-CORR ("<<i<<","<<l<<"): "<<ac[i][l]<<std::endl;
		}
	}
	for (int i = 0 ; i< sampler.dimension; i++)
		free(temp[i]);
	free(temp);
	deallocate_sampler_mem(&sampler);
}


