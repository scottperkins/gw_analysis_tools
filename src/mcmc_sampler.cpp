#include "mcmc_sampler.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "util.h"
#include "mcmc_sampler_internals.h"

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
		void (*fisher)(double *param, int dimension, double **fisher)	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
		)
{
	//random number generator initialization
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r = gsl_rng_alloc(T);

	//Array holding the probability of each 
	//type of step - Gaussian, differential evolution, MMALA, Fisher
	sampler sampler;
	bool fisher_bool;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL)
		fisher_bool = false;
	else 
		fisher_bool = true;
	
	if(!fisher_bool)//Obviously must add up to 1
	{
		//sampler.step_prob[0]=.7;
		//sampler.step_prob[1]=.3;
		//Testing
		sampler.step_prob[0]=1.;
		sampler.step_prob[1]=0.;
		sampler.step_prob[2]=0;
		sampler.step_prob[3]=0;
	}
	else
	{
		sampler.step_prob[0]=.3;
		sampler.step_prob[1]=.2;
		sampler.step_prob[2]=.3;
		sampler.step_prob[3]=.2;

	}
	//Split probabilities into boundaries for if-else loop
	sampler.prob_boundaries[0] = sampler.step_prob[0];
	sampler.prob_boundaries[1] = sampler.step_prob[1]+sampler.prob_boundaries[0];
	sampler.prob_boundaries[2] = sampler.step_prob[2]+sampler.prob_boundaries[1];
	sampler.prob_boundaries[3] = sampler.step_prob[3]+sampler.prob_boundaries[2];
	
	//Construct sampler structure
	sampler.lp = log_prior;
	sampler.ll = log_likelihood;
	sampler.fish = fisher;
	sampler.swp_freq = swp_freq;
	sampler.chain_temps = chain_temps;
	sampler.chain_N = chain_N;
	sampler.dimension = dimension;
	sampler.r = r;

	int i,j, k=0;
	int swp_accepted, swp_rejected;

	while (k<N_steps){
		for (j=0; j<chain_N; j++)
		{
			for (i = 0 ; i< sampler.swp_freq;i++)
			{
				mcmc_step(&sampler, output[j][i], output[j][i+k]);	
			}
		}
		k+= sampler.swp_freq;
		chain_swap(&sampler, output, k, &swp_accepted, &swp_rejected);
		printProgress((double)k/N_steps);	
	}
}
		
