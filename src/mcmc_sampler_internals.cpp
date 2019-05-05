#include "mcmc_sampler_internals.h"
#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/*! \file
 * File containing definitions for all the internal, generic mcmc subroutines
 */

/*! \brief interface function between the sampler and the internal step functions
 */
int mcmc_step(sampler *sampler, double *current_param, double *next_param, int chain_number)
{
	//Random number to determine type of step
	double alpha = gsl_rng_uniform(sampler->r);

	double proposed_param[sampler->dimension];

	
	if (alpha<sampler->prob_boundaries[0])
	{
		gaussian_step(sampler, current_param, proposed_param);
	}
	else if (alpha<sampler->prob_boundaries[1])
	{
		diff_ev_step(sampler, current_param, proposed_param);
	}
	else if (alpha<sampler->prob_boundaries[2])
	{
		mmala_step(sampler, current_param, proposed_param);
	}
	else 
	{
		fisher_step(sampler, current_param, proposed_param);
	}
	
	//Calculate log_likelihood and log prior
	double current_ll = sampler->ll(current_param, sampler->dimension);
	current_ll = (current_ll )/sampler->chain_temps[chain_number];
	double current_lp = sampler->lp(current_param, sampler->dimension);
	double proposed_ll = sampler->ll(proposed_param, sampler->dimension);
	proposed_ll = (proposed_ll )/sampler->chain_temps[chain_number];
	double proposed_lp = sampler->lp(proposed_param, sampler->dimension);

	//Calculate MH ratio
	double MH_ratio = std::exp(current_ll-proposed_ll+current_lp - proposed_lp);

	int i;
	if (MH_ratio>1.)
	{
		for ( i=0;i<sampler->dimension; i ++)
		{
			next_param[i] = proposed_param[i];
			return 1;
		}
	}	
	else 
	{
		//Random number to determine step acceptance
		double beta = gsl_rng_uniform(sampler->r);
		if(MH_ratio< beta){
			for ( i=0;i<sampler->dimension; i ++)
			{
				next_param[i] = current_param[i];
			}
			return -1;
		}	
		else
		{
			for ( i=0;i<sampler->dimension; i ++)
			{
				next_param[i] = proposed_param[i];
			}
			return 1;
		}		
	
	}
}
	

/*! \brief Straight gaussian step
 */
void gaussian_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param /**< [out] Proposed position in parameter space*/
		)
{
	int i ;
	double alpha = gsl_rng_uniform(sampler->r);
	for (i=0;i<sampler->dimension;i++){
		proposed_param[i] = gsl_ran_gaussian(sampler->r, alpha)+current_param[i];
	}
}

/*!\brief Fisher informed gaussian step
 */
void fisher_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param /**< [out] Proposed position in parameter space*/
		)
{
	
}

/*!\brief MMALA informed step
 */
void mmala_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param /**< [out] Proposed position in parameter space*/
		)
{
	
}

/*!\brief differential evolution informed step
 */
void diff_ev_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param /**< [out] Proposed position in parameter space*/
		)
{
	
}

/*! \brief subroutine to perform chain comparison for parallel tempering
 */
void chain_swap(sampler *sampler, /**<sampler struct*/
		double ***output, /**<output vector containing chains*/
		int step_num,	  /**<current step number*/ 
		int *swp_accepted,
		int *swp_rejected
		)
{
	for (int i =1; i < sampler->chain_N; i ++)
	{
		int success = single_chain_swap(sampler, output[i-1][step_num],output[i][step_num], i-1, i);
		if(success==1)
			*swp_accepted++;
		else
			*swp_rejected++;
	
	}
}

/*! \brief subroutine to actually swap two chains
 */
int single_chain_swap(sampler *sampler, /**< sampler structure*/
			double *swappe, /**< parameter position of chain that could be changed*/
			double *swapper, /**< chain that is not swapped, but provides parameters to be swapped by the other chain*/
			int Te_num,	/**<number of chain swappe in chain_temps*/
			int Tr_num	/**<number of chain swapper in chain_temps*/
			)
{
	double swappe_ll =  sampler->ll(swappe, sampler->dimension);
	double swapper_ll =  sampler->ll(swapper, sampler->dimension);
	//double swappe_lp =  sampler->lp(swappe, sampler->dimension);
	//double swapper_lp =  sampler->lp(swapper, sampler->dimension);
	double Te = sampler->chain_temps[Te_num];
	double Tr = sampler->chain_temps[Tr_num];
	double MH_ratio = std::exp(swappe_ll/Tr + swapper_ll/Te - swappe_ll/Te - swapper_ll/Tr);
	if(MH_ratio>1.)
	{
		for(int i =0; i < sampler->dimension;i++)
		{
			swappe[i] = swapper[i];
		}
		return 1;
	}
	else
	{
		double alpha = gsl_rng_uniform(sampler->r);
		if(MH_ratio<alpha){
			return -1;
		}
		else
		{
			for(int i =0; i < sampler->dimension;i++)
			{
				swappe[i] = swapper[i];
			}
			return 1;
		}
			
	}
}

