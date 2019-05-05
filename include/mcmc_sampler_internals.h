#ifndef MCMC_SAMPLER_INTERNALS_H
#define MCMC_SAMPLER_INTERNALS_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/*! \file
 * Internal functions of the generic MCMC sampler (nothing specific to GW)
 */

typedef double (*log_prior)(double *param, int dimension);
typedef double (*log_likelihood)(double *param, int dimension);
typedef void (*fisher)(double *param, int dimension,double **fisher);

/*! Structure storing everything that defines an instance of the sampler
 */
struct sampler
{
	double step_prob[4];
	double prob_boundaries[4];
	double *chain_temps;
	double swp_freq;
	int chain_N;
	int dimension;
	log_prior lp;
	log_likelihood ll;
	fisher fish;
	gsl_rng * r;
		
};

int mcmc_step(sampler *sampler, double *current_param,double *next_param);

void gaussian_step(sampler *sampler, double *current_param,double *proposed_param);

void fisher_step(sampler *sampler,double *current_param, double *proposed_param);

void mmala_step(sampler *sampler,double *current_param, double *proposed_param);

void diff_ev_step(sampler *sampler,double *current_param, double *proposed_param);

void chain_swap(sampler *sampler, double ***output, int step_num,int *swp_accepted, int *swp_rejected);

int single_chain_swap(sampler *sampler, double *swappe, double *swapper,int Te_num, int Tr_num);

#endif
