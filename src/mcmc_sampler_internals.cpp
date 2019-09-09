#include "mcmc_sampler_internals.h"
#include "autocorrelation.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <eigen3/Eigen/Eigen>
#include <limits>
#include <iomanip>
#include <fftw3.h>
#include <stdio.h>

/*! \file
 * File containing definitions for all the internal, generic mcmc subroutines
 */

/*! \brief interface function between the sampler and the internal step functions
 */
int mcmc_step(sampler *sampler, double *current_param, double *next_param, int *current_status, int *next_status, int chain_number)
{
	//Random number to determine type of step
	double alpha = gsl_rng_uniform(sampler->rvec[chain_number]);
	

	double proposed_param[sampler->max_dim];
	int proposed_status[sampler->max_dim];

	int step;	
	//Determine which step to take and calculate proposed coord.
	if (alpha<sampler->prob_boundaries[chain_number][0])
	{
		gaussian_step(sampler, current_param, proposed_param, current_status, proposed_status, chain_number);
		sampler->num_gauss[chain_number]+=1;
		step = 0;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][1])
	{
		diff_ev_step(sampler, current_param, proposed_param, current_status, proposed_status, chain_number);
		sampler->num_de[chain_number]+=1;
		step = 1;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][2])
	{
		mmala_step(sampler, current_param, proposed_param, current_status, proposed_status);
		sampler->num_mmala[chain_number]+=1;
		step= 2;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][3])
	{
		fisher_step(sampler, current_param, proposed_param, current_status, proposed_status, chain_number);
		sampler->num_fish[chain_number]+=1;
		step = 3;
	}
	else 
	{
		RJ_step(sampler, current_param, proposed_param, current_status, proposed_status, chain_number);
		sampler->num_RJstep[chain_number]+=1;
		step = 4;
	}
	
	double current_lp = sampler->lp(current_param, current_status,sampler->max_dim, chain_number);
	double proposed_lp = sampler->lp(proposed_param,proposed_status, sampler->max_dim, chain_number);
	double current_ll=0, proposed_ll=0;
	double MH_ratio;
	double power;

	//Check if proposed step is out of range by calculating prior
	//if out of range, reject step
	if(current_lp == limit_inf || proposed_lp == limit_inf){
		MH_ratio = limit_inf;
	}
	else{
		//Calculate log_likelihood and log prior
		current_ll = sampler->current_likelihoods[chain_number];
		proposed_ll = sampler->ll(proposed_param, proposed_status,sampler->max_dim,chain_number);
		proposed_ll = (proposed_ll )/sampler->chain_temps[chain_number];
		//Calculate MH ratio
		if(std::isnan(proposed_ll)){
			MH_ratio = limit_inf;
		}
		else{
			MH_ratio = -current_ll+proposed_ll-current_lp + proposed_lp;
		}
	}
	if(sampler->chain_temps[chain_number]==1 && sampler->chain_pos[chain_number]%1000==0){
		std::cout<<proposed_ll<<" "<<current_ll<<std::endl;
		std::cout<<proposed_lp<<" "<<current_lp<<std::endl;
		std::cout<<chain_number<<" "<<sampler->chain_pos[chain_number]<<std::endl;
	}
	int i;
	//Random number to determine step acceptance
	double beta = log(gsl_rng_uniform(sampler->rvec[chain_number]));
	if(MH_ratio< beta){
		for ( i=0;i<sampler->max_dim; i ++)
		{
			next_param[i] = current_param[i];
			next_status[i] = current_status[i];
		}
		assign_ct_m(sampler, step,chain_number);
		return -1;
	}	
	else
	{
		for ( i=0;i<sampler->max_dim; i ++)
		{
			next_param[i] = proposed_param[i];
			next_status[i] = proposed_status[i];
		}
		assign_ct_p(sampler, step, chain_number);
		sampler->current_likelihoods[chain_number] = proposed_ll;
		return 1;
	}		
	
}
	

/*! \brief Straight gaussian step
 */
void gaussian_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int chain_id
		)
{
	double alpha = sampler->randgauss_width[chain_id][0];
	for (int i=0;i<sampler->max_dim;i++){
		if(current_status[i] == 1){
			proposed_param[i] = gsl_ran_gaussian(sampler->rvec[chain_id], alpha)+current_param[i];
		}
		else{
			proposed_param[i] = 0;
		}
		proposed_status[i] = current_status[i];
	}
}

/*!\brief Fisher informed gaussian step
 */
void fisher_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int chain_index
		)
{
	if(!sampler->RJMCMC || sampler->min_dim ==0){
		//Check whether or not we need to update the fisher
		if(sampler->fisher_update_ct[chain_index]==sampler->fisher_update_number)
			update_fisher(sampler, current_param, current_status,chain_index);	

		//update the count of steps since last fisher update
		sampler->fisher_update_ct[chain_index] += 1;
		
		//beta determines direction to step in eigen directions
		int beta = (int)((sampler->max_dim)*(gsl_rng_uniform(sampler->rvec[chain_index])));
		
		double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
					 sampler->randgauss_width[chain_index][3]);

		double scaling = 0.0;
		//ensure the steps aren't ridiculous
		if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
		else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}

		else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
					sampler->chain_temps[chain_index];}
		//Take step
		for(int i =0; i< sampler->max_dim;i++)
		{
			proposed_param[i] = current_param[i] +
				alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
			proposed_status[i] = current_status[i];
		}
	}
	//If RJPTMCMC and there's a base model, use the fisher for the base model, and gaussian steps for the modifications
	else {
		//Check whether or not we need to update the fisher
		if(sampler->fisher_update_ct[chain_index]==sampler->fisher_update_number)
			update_fisher(sampler, current_param, current_status,chain_index);	

		//update the count of steps since last fisher update
		sampler->fisher_update_ct[chain_index] += 1;
		
		//beta determines direction to step in eigen directions
		int beta = (int)((sampler->min_dim)*(gsl_rng_uniform(sampler->rvec[chain_index])));
		
		double alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
					 sampler->randgauss_width[chain_index][3]);

		double scaling = 0.0;
		//ensure the steps aren't ridiculous
		if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
		else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}

		else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
					sampler->chain_temps[chain_index];}
		//Take step
		for(int i =0; i< sampler->min_dim;i++)
		{
			proposed_param[i] = current_param[i] +
				alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
			proposed_status[i] = current_status[i];
		}
		//Generate new step for gaussian steps, using gaussian width	
		alpha = gsl_ran_gaussian(sampler->rvec[chain_index],
					 sampler->randgauss_width[chain_index][0]);
		for(int i =sampler->min_dim; i< sampler->max_dim;i++)
		{
			if(current_status[i] == 1){
				proposed_param[i] = alpha+current_param[i];
				proposed_status[i] = current_status[i];
			}
			else{
				proposed_param[i] = 0;
				proposed_status[i] = current_status[i];
			}
		}
		
	}

}


void update_fisher(sampler *sampler, double *current_param, int *param_status, int chain_index)
{
	int local_dim = sampler->max_dim;
	if(sampler->RJMCMC && sampler->min_dim !=0) local_dim=sampler->min_dim;
	//Fisher calculation
	double **fisher=(double **)malloc(sizeof(double*)*local_dim);	
	for (int i =0; i<local_dim;i++){
		fisher[i] = (double*)malloc(sizeof(double)*local_dim);
	}
	sampler->fish(current_param, param_status,local_dim, fisher,chain_index);

	//Convert to 1D array for Eigen
	double *oneDfisher=(double *)malloc(sizeof(double)*local_dim*local_dim);
	for (int i =0; i<local_dim;i++){
		for (int j = 0; j<local_dim; j++){
			oneDfisher[local_dim*i+j] = fisher[i][j];///
		}
		
	}
	
	//Find eigen vectors and eigen values
	Eigen::Map<Eigen::MatrixXd> m(oneDfisher,local_dim,local_dim);
 	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m);
 	Eigen::MatrixXd eigen_vecs = eigensolver.eigenvectors();
 	Eigen::VectorXd eigen_vals = eigensolver.eigenvalues();
	
	//##############################################################
	//Numerical matrix inversion can be tricky - catch nans here and replace
	//with gaussian step just to not have the program crash because of 
	//one position in parameter space
	//
	//To see impact of nans, variable is stored in sampler->nan_counter
	//##############################################################
	int nansum = 0;
	for(int j = 0 ; j<local_dim; j++){
		for(int i =0; i<local_dim; i++)
			nansum+= std::isnan(eigen_vecs.col(j)(i));	
		nansum+= std::isnan(eigen_vals(j));
	}
	if(!nansum){
		for (int i =0; i < local_dim; i++)
		{
			for(int j = 0; j<local_dim; j++)
			{
				sampler->fisher_vecs[chain_index][i][j] = eigen_vecs.col(i)(j);
			}
			sampler->fisher_vals[chain_index][i]=eigen_vals[i];
		}
		sampler->fisher_update_ct[chain_index]=0;
	}
	else{ 
		sampler->fisher_update_ct[chain_index]=sampler->fisher_update_number-1;
		sampler->nan_counter[chain_index]+=1;
	}

	for (int i =0; i<local_dim;i++){
		free(fisher[i]);
	}
	free(fisher);
	free(oneDfisher);
}

/*!\brief MMALA informed step -- Currently not supported
 */
void mmala_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status
		)
{
	
}

/*!\brief differential evolution informed step
 *
 * Differential evolution uses the past history of the chain to inform the proposed step:
 *
 * Take the difference of two random, accepted previous steps and step along that with some step size,
 * determined by a gaussian
 */
void diff_ev_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int *current_status,
		int *proposed_status,
		int chain_id
		)
{
	double *eff_history_coord = new double[sampler->max_dim];
	int *eff_history_status = new int[sampler->max_dim];
	double *eff_history_coord2 = new double[sampler->max_dim];
	int *eff_history_status2 = new int[sampler->max_dim];
	if(sampler->RJMCMC){
		//Pick a history member
		int i = (int)((sampler->history_length-1)
			*(gsl_rng_uniform(sampler->rvec[chain_id])));
		//Second position ID
		int j;
		do{
			j=(int)((sampler->history_length-1)*
				(gsl_rng_uniform(sampler->rvec[chain_id])));	
		}while(j==i);

		//We'll almost definitely need to do some transformation to it

		//Since there's no way to feasibly store a substantial 
		//population of history items for every combination of dimensions:
		//I'll assume the additional dimensions have small effects on the model
		//described by dimensions 0-min_dim, so I pick an element from 
		//the history file that has the smallest deviation from the base model, 
		//but still has the additional dimension. This means that the extra dimensions
		//beyond min_dim are indepent, but if the affects are indeed small, the benefit 
		//of DE should still be recovered (this is only a proposal, after all)
		if(sampler->min_dim != 0){
			RJ_smooth_history(sampler, current_param,current_status, i, eff_history_coord, 
				eff_history_status,chain_id);
			RJ_smooth_history(sampler, current_param, current_status, j, eff_history_coord2, 
				eff_history_status2,chain_id);
		}
		//For models that are composed of discrete models (ie model A or B 
		//and min_dim = 0), I just need to continue picking history members 
		//until I get two that match the correct model.If two elements that 
		//contain a given model, we just use a gaussian step and abandon DE 
		//for this step.
		else{

		}
		
	}
	//Regular PTMCMC
	else{
		//First position ID
		int i = (int)((sampler->history_length-1)*(gsl_rng_uniform(sampler->rvec[chain_id])));
		//Second position ID
		int j;
		do{
			j=(int)((sampler->history_length-1)*(gsl_rng_uniform(sampler->rvec[chain_id])));	
		}while(j==i);

		for(int k = 0 ; k <sampler->max_dim; k++){
			eff_history_coord[k] = sampler->history[chain_id][i][k];
			eff_history_coord2[k] = sampler->history[chain_id][j][k];
			eff_history_status[k] = sampler->param_status[chain_id][0][k];
			eff_history_status2[k] = sampler->param_status[chain_id][0][k];
			//eff_history_status[k] = sampler->history_status[chain_id][i][k];
			//eff_history_status2[k] = sampler->history_status[chain_id][j][k];
		}	
			
	}

	
	double alpha = .1;
	double beta = gsl_rng_uniform(sampler->rvec[chain_id]);
	if(beta<.9)
		alpha=gsl_ran_gaussian(sampler->rvec[chain_id],sampler->randgauss_width[chain_id][1]);
	for (int k = 0; k<sampler->max_dim; k++)
	{
//		proposed_param[k] = current_param[k] + alpha*
			//(sampler->history[chain_id][i][k]-sampler->history[chain_id][j][k]);
		if(current_status[k] == 1){
			proposed_param[k] = current_param[k] + alpha *
				(eff_history_coord[k] - eff_history_coord2[k]);
		}
		else{
			proposed_param[k] = 0;
		}
		proposed_status[k]=eff_history_status[k];
	}
	delete [] eff_history_coord;
	delete [] eff_history_status;
	delete [] eff_history_coord2;
	delete [] eff_history_status2;
}

void RJ_smooth_history(sampler *sampler, /**<Current sampler */
	double *current_param,/**<Current parameters to match*/
	int *current_param_status,/**<Current parameters to match*/
	int base_history_id, /**<Original history element*/
	double *eff_history_coord, /**<[out] Modified history coord*/
	int *eff_history_status, /**<[out] Modified History status*/
	int chain_id/**<Chain ID of the current chain*/
	)
{
	//Check which dimensions are included in the history element
	int *bad_ids = new int[sampler->max_dim-sampler->min_dim];
	int id_count = 0;
	for (int j = sampler->min_dim ; j<sampler->max_dim; j++){
		//if(current_param_status[j] != sampler->history_status[chain_id][base_history_id][j]){
		if(current_param_status[j] == 1 
			&& sampler->history_status[chain_id][base_history_id][j] != 1){

			bad_ids[id_count] = j;	
			id_count++;
		}
	}
	//Copy over all history params, even if they are zero to ensure no memory errors
	for(int j = 0 ; j<sampler->max_dim; j ++){
		if(current_param_status[j] == 1 ){
			eff_history_coord[j] = sampler->history[chain_id][base_history_id][j];
			eff_history_status[j] = 1;
		}
		else{
			eff_history_coord[j] =0;
			eff_history_status[j] =0;
		}

	}
	//Go back, and overwrite the params not part of the history params
	int min_dev_id = 0;
	double min_dev = 0;
	double current_dev = 0;
	bool at_least_one = false; 
	//For each bad_id, loop through history and find the smallest deviation
	for(int k = 0 ; k<id_count; k++){
		for(int l = 0 ; l<sampler->history_length; l++){
			if(sampler->history_status[chain_id][l][bad_ids[k]]==1){
				for(int j=0 ;j<sampler->min_dim; j++){
						current_dev+= (sampler->history[chain_id][l][j] - sampler->history[chain_id][base_history_id][j])/sampler->history[chain_id][base_history_id][j];
				}
				current_dev/= sampler->min_dim;
				if(!at_least_one){
					min_dev = current_dev;		
					min_dev_id = l;		
					at_least_one = true;
				}
				else{
					if(current_dev<min_dev){
						min_dev = current_dev;		
						min_dev_id = l;		
					}
				}
			}
			current_dev = 0 ;
		}
		if(at_least_one){
			eff_history_coord[k] = sampler->history[chain_id][min_dev_id][k];
		}
		else{
			//Gaussian step for this dimension
			double alpha = gsl_ran_gaussian(sampler->rvec[chain_id],1);
			eff_history_coord[k] = alpha+current_param[k];
		}
		at_least_one = false;
	}
	
	
	delete [] bad_ids;

}

/*! \brief RJ-proposal step for trans-dimensional MCMCs
 *
 * This extra step may seem unnecessary, I'm just adding it in in case the extra flexibility is useful in the future for preprocessing of the chain before sending it to the user's RJ_proposal
 */
void RJ_step(sampler *sampler, /**< sampler*/
	double *current_param, /**< current coordinates in parameter space*/
	double *proposed_param, /**<[out] Proposed coordinates in parameter space*/
	int *current_status, /**< Current status of parameters*/
	int *proposed_status, /**<[out] Proposed status of parameters*/
	int chain_number/**< chain mumber*/
	)
{
	sampler->rj(current_param, proposed_param, current_status, proposed_status, sampler->max_dim,chain_number, sampler->randgauss_width[chain_number][4]);
}

/*! \brief subroutine to perform chain comparison for parallel tempering
 *
 * The total output file is passed, and the chains are swapped sequentially
 *
 * This is the routine for ``Deterministic'' sampling (parallel or sequential, 
 * but not pooled)
 */
void chain_swap(sampler *sampler, /**<sampler struct*/
		double ***output, /**<output vector containing chains*/
		int ***param_status,/**< Parameter status*/
		int step_num,	  /**<current step number*/ 
		int *swp_accepted,
		int *swp_rejected
		)
{
	for (int i =0; i < sampler->chain_N-1; i ++)
	{
		int success = single_chain_swap(sampler, output[i][step_num],output[i+1][step_num],
			param_status[i][step_num],param_status[i+1][step_num], i, i+1);
		if(success==1){
			sampler->swap_accept_ct[i]+=1;
			sampler->swap_accept_ct[i+1]+=1;
			//PT dynamics
			if(sampler->PT_alloc)
				sampler->A[i+1] = 1;
		}
		else{
			sampler->swap_reject_ct[i]+=1;
			sampler->swap_reject_ct[i+1]+=1;
			//PT dynamics
			if(sampler->PT_alloc)
				sampler->A[i+1] = 0;
		}
	
	}
}

/*! \brief subroutine to actually swap two chains
 *
 * This is the more general subroutine, which just swaps the two chains passed to the function
 */
int single_chain_swap(sampler *sampler, /**< sampler structure*/
			double *chain1, /**< parameter position of chain that could be changed*/
			double *chain2, /**< chain that is not swapped, but provides parameters to be swapped by the other chain*/
			int *chain1_status,/**<Parameter status array for chain1*/
			int *chain2_status,/**<Parameter status array for chain2*/
			int T1_index,	/**<number of chain swappe in chain_temps*/
			int T2_index	/**<number of chain swapper in chain_temps*/
			)
{
	//Unpack parameters
	double T1 = sampler->chain_temps[T1_index];
	double T2 = sampler->chain_temps[T2_index];
	double ll1 =  T1*sampler->current_likelihoods[T1_index];
	double ll2 =  T2*sampler->current_likelihoods[T2_index];
	double pow = (ll1-ll2)/T2 - (ll1-ll2)/T1;
	double MH_ratio;
	MH_ratio = pow;
	//Averaging the two random numbers from each chains seed
	double alpha = log( (gsl_rng_uniform(sampler->rvec[T1_index])+gsl_rng_uniform(sampler->rvec[T2_index]))/2.);
	if (MH_ratio<alpha)
	{
		return -1;
	}	
	else
	{
		double temp[sampler->max_dim];
		int tempstat[sampler->max_dim];
		for(int i =0; i < sampler->max_dim;i++)
		{
			temp[i] = chain1[i];
			tempstat[i] = chain1_status[i];
			chain1[i] = chain2[i];
			chain1_status[i] = chain2_status[i];
			chain2[i]=temp[i];
			chain2_status[i]=tempstat[i];
		}
		double templl = sampler->current_likelihoods[T1_index];
		sampler->current_likelihoods[T1_index] = 
				T2/T1 * sampler->current_likelihoods[T2_index];
		sampler->current_likelihoods[T2_index] = T1/T2*templl;
		return 1;
	}

}

/*! \brief update and initiate probabilities for each variety of step
 *
 * Type 0: Gaussian step
 *
 * Type 1: Differential Evolution step
 *
 * Type 2: MMALA step (currently not supported)
 *
 * Type 3: Fisher step
 */
void assign_probabilities(sampler *sampler, int chain_index)
{

	//non-RJPTMCMC
	if(!sampler->RJMCMC){
		sampler->step_prob[chain_index][4]=.0;
		//no fisher and de not ready
		if(!sampler->fisher_exist && !sampler->de_primed[chain_index])//Obviously must add up to 1
		{
			sampler->step_prob[chain_index][0]=1.;
			sampler->step_prob[chain_index][1]=0.;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=0;
		}
		//fisher available, but de not yet ready
		else if (sampler->fisher_exist && !sampler->de_primed[chain_index])
		{
			//sampler->step_prob[chain_index][0]=.3;
			//sampler->step_prob[chain_index][1]=0;
			//sampler->step_prob[chain_index][2]=.3;
			//sampler->step_prob[chain_index][3]=.4;
			//Testing
			sampler->step_prob[chain_index][0]=.1;
			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.9;

		}
		//No fisher, but de ready
		else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
		{
			
			sampler->step_prob[chain_index][0]=.2;
			sampler->step_prob[chain_index][1]=.8;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.0;

		}
		//all methods available
		else
		{
			sampler->step_prob[chain_index][0]=.05;
			sampler->step_prob[chain_index][1]=.2;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.75;

		}
	}
	else{
		//no fisher and de not ready
		if(!sampler->fisher_exist && !sampler->de_primed[chain_index])//Obviously must add up to 1
		{
			sampler->step_prob[chain_index][0]=.8;
			sampler->step_prob[chain_index][1]=0.;
			sampler->step_prob[chain_index][2]=0;
			sampler->step_prob[chain_index][3]=0;
			sampler->step_prob[chain_index][4]=.2;
		}
		//fisher available, but de not yet ready
		else if (sampler->fisher_exist && !sampler->de_primed[chain_index])
		{
			sampler->step_prob[chain_index][0]=.05;
			sampler->step_prob[chain_index][1]=0;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.75;
			sampler->step_prob[chain_index][4]=.2;

		}
		//No fisher, but de ready
		else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
		{
			
			sampler->step_prob[chain_index][0]=.1;
			sampler->step_prob[chain_index][1]=.7;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.0;
			sampler->step_prob[chain_index][4]=.2;

		}
		//all methods available
		else
		{
			sampler->step_prob[chain_index][0]=.05;
			sampler->step_prob[chain_index][1]=.15;
			sampler->step_prob[chain_index][2]=.0;
			sampler->step_prob[chain_index][3]=.6;
			sampler->step_prob[chain_index][4]=.2;

		}
	}
	//Split probabilities into boundaries for if-else loop
	sampler->prob_boundaries[chain_index][0] = sampler->step_prob[chain_index][0];
	sampler->prob_boundaries[chain_index][1] = sampler->step_prob[chain_index][1]+sampler->prob_boundaries[chain_index][0];
	sampler->prob_boundaries[chain_index][2] = sampler->step_prob[chain_index][2]+sampler->prob_boundaries[chain_index][1];
	sampler->prob_boundaries[chain_index][3] = sampler->step_prob[chain_index][3]+sampler->prob_boundaries[chain_index][2];
	sampler->prob_boundaries[chain_index][4] = sampler->step_prob[chain_index][4]+sampler->prob_boundaries[chain_index][3];
}	

/*! \brief Copies contents of one chain to another
 *
 * Transfers id_source in samplerptr_source to id_dest samplerptr_dest 
 *
 * NOTE: This copies the VALUE, not the reference. This could be expensive, so use with caution
 *
 * id_dest is ERASED
 *
 * samplerptr_dest and samplerptr_source MUST have the same dimension, the same sampling details (like having or not having a fisher) etc
 * 	
 * samplerptr_dest must be previously allocated properly
 *
 * As output is the largest transfer by far, the transfer_output flag can be used to allow the user to handle that manually.
 */
void transfer_chain(sampler *samplerptr_dest,sampler *samplerptr_source, int id_dest, int id_source, bool transfer_output)
{	
	//Position and output and likelihood and temp
	samplerptr_dest->chain_temps[id_dest] = 
		samplerptr_source->chain_temps[id_source];
	samplerptr_dest->chain_pos[id_dest] = samplerptr_source->chain_pos[id_source];
	if(transfer_output){
		for(int i =0 ;i <=samplerptr_dest->chain_pos[id_dest]; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->output[id_dest][i][j] = 
					samplerptr_source->output[id_source][i][j];
				//Might only need to do this for RJMCMC
				if(samplerptr_dest->RJMCMC){
					samplerptr_dest->param_status[id_dest][i][j] = 
						samplerptr_source->param_status[id_source][i][j];
				}
			}
		}
	}
	samplerptr_dest->current_likelihoods[id_dest] = samplerptr_source->current_likelihoods[id_source];

	//Histories
	samplerptr_dest->de_primed[id_dest] 
		= samplerptr_source->de_primed[id_source];
	samplerptr_dest->current_hist_pos[id_dest] 
		= samplerptr_source->current_hist_pos[id_source];
	//If History not filled, only copy over till current pos
	if(!samplerptr_dest->de_primed[id_dest]){
		for(int i =0 ;i <=samplerptr_dest->current_hist_pos[id_dest]; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->history[id_dest][i][j] = 
					samplerptr_source->history[id_source][i][j];
			}
		}
	}
	//Else, copy over whole history
	else{
		for(int i =0 ;i <samplerptr_dest->history_length; i++){
			for(int j =0; j<samplerptr_dest->max_dim; j++){
				samplerptr_dest->history[id_dest][i][j] = 
					samplerptr_source->history[id_source][i][j];
			}
		}
	}

	//Step parameters
	samplerptr_dest->gauss_last_accept_ct[id_dest] 
		= samplerptr_source->gauss_last_accept_ct[id_source]; 
	samplerptr_dest->gauss_last_reject_ct[id_dest] 
		= samplerptr_source->gauss_last_reject_ct[id_source]; 
	samplerptr_dest->de_last_accept_ct[id_dest] 
		= samplerptr_source->de_last_accept_ct[id_source]; 
	samplerptr_dest->de_last_reject_ct[id_dest] 
		= samplerptr_source->de_last_reject_ct[id_source]; 
	samplerptr_dest->fish_last_accept_ct[id_dest] 
		= samplerptr_source->fish_last_accept_ct[id_source]; 
	samplerptr_dest->fish_last_reject_ct[id_dest] 
		= samplerptr_source->fish_last_reject_ct[id_source]; 

	samplerptr_dest->gauss_accept_ct[id_dest] 
		= samplerptr_source->gauss_accept_ct[id_source]; 
	samplerptr_dest->gauss_reject_ct[id_dest] 
		= samplerptr_source->gauss_reject_ct[id_source]; 
	samplerptr_dest->de_accept_ct[id_dest] 
		= samplerptr_source->de_accept_ct[id_source]; 
	samplerptr_dest->de_reject_ct[id_dest] 
		= samplerptr_source->de_reject_ct[id_source]; 
	samplerptr_dest->fish_accept_ct[id_dest] 
		= samplerptr_source->fish_accept_ct[id_source]; 
	samplerptr_dest->fish_reject_ct[id_dest] 
		= samplerptr_source->fish_reject_ct[id_source]; 
	samplerptr_dest->mmala_accept_ct[id_dest] 
		= samplerptr_source->mmala_accept_ct[id_source]; 
	samplerptr_dest->mmala_reject_ct[id_dest] 
		= samplerptr_source->mmala_reject_ct[id_source]; 
	for(int i =0 ;i<samplerptr_dest->types_of_steps; i++){
		samplerptr_dest->step_prob[id_dest][i] 
			= samplerptr_source->step_prob[id_source][i];
		samplerptr_dest->prob_boundaries[id_dest][i] 
			= samplerptr_source->prob_boundaries[id_source][i];
		samplerptr_dest->randgauss_width[id_dest][i] 
			= samplerptr_source->randgauss_width[id_source][i];
	}
	//Fisher
	if(samplerptr_dest->fisher_exist){
		if(samplerptr_source->fisher_update_ct[id_source] != samplerptr_source->fisher_update_number){
			for (int i =0 ; i<samplerptr_source->max_dim; i++){
				for (int j =0 ; j<samplerptr_source->max_dim; j++){
					samplerptr_dest->fisher_vecs[id_dest][i][j] = 
						samplerptr_source->fisher_vecs[id_source][i][j];
				}
				samplerptr_dest->fisher_vals[id_dest][i] = 
					samplerptr_source->fisher_vals[id_source][i];
			}
			samplerptr_dest->fisher_update_ct[id_dest] = samplerptr_source->fisher_update_ct[id_source];
		}
		else{
			samplerptr_dest->fisher_update_ct[id_dest] = samplerptr_dest->fisher_update_number;
		}
	}

	//Sampling parameters
	samplerptr_dest->waiting[id_dest] = samplerptr_source->waiting[id_source];
	samplerptr_dest->priority[id_dest] = samplerptr_source->priority[id_source];
	samplerptr_dest->ref_chain_status[id_dest] = samplerptr_source->ref_chain_status[id_source];
	samplerptr_dest->nan_counter[id_dest] = samplerptr_source->nan_counter[id_source];
	samplerptr_dest->num_gauss[id_dest] = samplerptr_source->num_gauss[id_source];
	samplerptr_dest->num_fish[id_dest] = samplerptr_source->num_fish[id_source];
	samplerptr_dest->num_de[id_dest] = samplerptr_source->num_de[id_source];
	samplerptr_dest->num_mmala[id_dest] = samplerptr_source->num_mmala[id_source];
	samplerptr_dest->swap_accept_ct[id_dest] = samplerptr_source->swap_accept_ct[id_source];
	samplerptr_dest->swap_reject_ct[id_dest] = samplerptr_source->swap_reject_ct[id_source];
	samplerptr_dest->step_accept_ct[id_dest] = samplerptr_source->step_accept_ct[id_source];
	samplerptr_dest->step_reject_ct[id_dest] = samplerptr_source->step_reject_ct[id_source];
	if(samplerptr_source->log_ll){
		for(int i =0 ; i<samplerptr_source->chain_pos[id_source]; i++)
			samplerptr_dest->ll_lp_output[id_dest][i][0] 
				= samplerptr_source->ll_lp_output[id_source][i][0];
	}
	if(samplerptr_source->log_lp){
		for(int i =0 ; i<samplerptr_source->chain_pos[id_source]; i++)
			samplerptr_dest->ll_lp_output[id_dest][i][1] 
				= samplerptr_source->ll_lp_output[id_source][i][1];
	}
	if(samplerptr_source->PT_alloc)
		samplerptr_dest->A[id_dest] = samplerptr_source->A[id_source];
}

/*! \brief Checks the status of a sampler for the stochastic sampling
 *
 * Just loops through the ref_chain_status variables
 */
bool check_sampler_status(sampler *samplerptr)
{
	for(int i =0 ; i<samplerptr->chain_N; i++){
		if(!samplerptr->ref_chain_status[i])
			return false;
	}
	return true;
}
/*! \brief Updates the step widths, shooting for 20% acceptance ratios for each type of step
 */
void update_step_widths(sampler *samplerptr, int chain_id)
{
	int j = chain_id;
	//update stepsize to maximize step efficiency
	//increases in stepsizes of 10%
	double frac, acc, rej;
	if(samplerptr->chain_pos[j]%samplerptr->check_stepsize_freq[j] == 0){
		//Gaussian
		if(samplerptr->step_prob[j][0]!= 0){
			acc = samplerptr->gauss_accept_ct[j] - samplerptr->gauss_last_accept_ct[j];	
			rej = samplerptr->gauss_reject_ct[j] - samplerptr->gauss_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][0] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][0] *=1.1;	
			}
			samplerptr->gauss_last_accept_ct[j]=samplerptr->gauss_accept_ct[j];
			samplerptr->gauss_last_reject_ct[j]=samplerptr->gauss_reject_ct[j];
		}	
		//de
		if(samplerptr->step_prob[j][1]!= 0){
			acc = samplerptr->de_accept_ct[j] - samplerptr->de_last_accept_ct[j];	
			rej = samplerptr->de_reject_ct[j] - samplerptr->de_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1] *=1.1;	
			}
			samplerptr->de_last_accept_ct[j]=samplerptr->de_accept_ct[j];
			samplerptr->de_last_reject_ct[j]=samplerptr->de_reject_ct[j];
		}	
		//fisher
		if(samplerptr->step_prob[j][3]!= 0){
			acc = samplerptr->fish_accept_ct[j] - samplerptr->fish_last_accept_ct[j];	
			rej = samplerptr->fish_reject_ct[j] - samplerptr->fish_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3] *=1.1;	
			}
			samplerptr->fish_last_accept_ct[j]=samplerptr->fish_accept_ct[j];
			samplerptr->fish_last_reject_ct[j]=samplerptr->fish_reject_ct[j];
		}	
		//RJ
		if(samplerptr->update_RJ_width){
			if(samplerptr->step_prob[j][4]!= 0){
				acc = samplerptr->RJstep_accept_ct[j] - samplerptr->RJstep_last_accept_ct[j];	
				rej = samplerptr->RJstep_reject_ct[j] - samplerptr->RJstep_last_reject_ct[j];	
				frac = acc / (acc + rej);
				if(frac<samplerptr->min_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][4] *=.9;	
				}
				else if(frac>samplerptr->max_target_accept_ratio[j]){
					samplerptr->randgauss_width[j][4] *=1.1;	
				}
				samplerptr->RJstep_last_accept_ct[j]=samplerptr->RJstep_accept_ct[j];
				samplerptr->RJstep_last_reject_ct[j]=samplerptr->RJstep_reject_ct[j];
			}	
		}
		
	}
}

void allocate_sampler_mem(sampler *sampler)
{
	gsl_rng_env_setup();
	const gsl_rng_type *T= gsl_rng_default;

	int i;
	sampler->step_prob = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->prob_boundaries = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->de_primed = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->waiting = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->current_hist_pos = (int *)malloc(sizeof(int ) * sampler->chain_N);
	sampler->chain_pos = (int *)malloc(sizeof(int ) * sampler->chain_N);

	sampler->fish_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->mmala_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->mmala_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->RJstep_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->RJstep_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fisher_update_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * sampler->chain_N);

	sampler->nan_counter = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_gauss = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_fish = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_de = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_mmala = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->num_RJstep = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->priority = (int *)malloc(sizeof(int) * sampler->chain_N);

	sampler->current_likelihoods = (double *)malloc(sizeof(double) * sampler->chain_N);

	sampler->check_stepsize_freq = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->max_target_accept_ratio = (double *)malloc(sizeof(double) * sampler->chain_N);
	sampler->min_target_accept_ratio = (double *)malloc(sizeof(double) * sampler->chain_N);
	sampler->gauss_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->gauss_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->fish_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->de_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	if(sampler->update_RJ_width){
		sampler->RJstep_last_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
		sampler->RJstep_last_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	}
	sampler->randgauss_width = allocate_2D_array(sampler->chain_N, sampler->types_of_steps); //Second dimension is types of steps

	//RJ parameters -- initialize status array with 1's for now, then repopulate with initial position
	if(sampler->RJMCMC){
		//sampler->param_status = allocate_3D_array_int(sampler->chain_N, sampler->N_steps,sampler->max_dim);	
		sampler->history_status = allocate_3D_array_int(sampler->chain_N, sampler->history_length,sampler->max_dim);	
	}
	else{
		//Need one copy per chain at least -- just in case of race conditions
		sampler->param_status = (int ***)malloc(sizeof(int **)*sampler->chain_N);
		sampler->history_status = (int ***)malloc(sizeof(int **)*sampler->chain_N);
	}

	sampler->ref_chain_status = (bool *) malloc(sizeof(bool)*sampler->chain_N);

	for (i =0; i<sampler->chain_N; i++)
	{
		//designate T=1 chains as reference chains
		if(sampler->chain_temps[i] == 1){
			sampler->ref_chain_status[i] = false;
		}
		else{
			sampler->ref_chain_status[i] = true;
		}
		//initial value set to one's
		if(sampler->RJMCMC){
			//for(int j = 0 ; j< sampler->max_dim ;j ++){
			//	sampler->param_status[i][0][j] = 1;
			//}
		}
		else{
			sampler->param_status[i] = (int **)malloc(sizeof(int *)*sampler->N_steps);
			sampler->param_status[i][0] = (int *)malloc(sizeof(int )*sampler->max_dim);
			for(int k = 0 ; k<sampler->max_dim; k++){
				sampler->param_status[i][0][k] = 1;
			}
			for(int k = 1 ; k<sampler->N_steps; k ++){
				sampler->param_status[i][k] = sampler->param_status[i][0];
			}
			sampler->history_status[i] = (int **)malloc(sizeof(int *)*sampler->history_length);
			sampler->history_status[i][0] = (int *)malloc(sizeof(int )*sampler->max_dim);
			for(int k = 0 ; k<sampler->max_dim; k++){
				sampler->history_status[i][0][k] = 1;
			}
			for(int k = 1 ; k<sampler->history_length; k ++){
				sampler->history_status[i][k] = sampler->history_status[i][0];
			}
		}

		sampler->step_prob[i] = (double *)malloc(sizeof(double)*sampler->types_of_steps);
		sampler->prob_boundaries[i] = (double *)malloc(sizeof(double)*sampler->types_of_steps);
		sampler->de_primed[i] = false;
		sampler->current_hist_pos[i] = 0;
		sampler->fish_accept_ct[i]=0;
		sampler->fish_reject_ct[i]=0;
		sampler->de_accept_ct[i]=0;
		sampler->de_reject_ct[i]=0;
		sampler->gauss_accept_ct[i]=0;
		sampler->gauss_reject_ct[i]=0;
		sampler->mmala_accept_ct[i]=0;
		sampler->mmala_reject_ct[i]=0;
		sampler->RJstep_accept_ct[i]=0;
		sampler->RJstep_reject_ct[i]=0;
		sampler->fisher_update_ct[i] = sampler->fisher_update_number;
		sampler->chain_pos[i]=0;
		sampler->waiting[i]=true;
		sampler->swap_accept_ct[i]=0;
		sampler->swap_reject_ct[i]=0;
		sampler->step_accept_ct[i]=0;
		sampler->step_reject_ct[i]=0;

		sampler->nan_counter[i]=0;
		sampler->num_gauss[i]=0;
		sampler->num_fish[i]=0;
		sampler->num_de[i]=0;
		sampler->num_mmala[i]=0;
		sampler->num_RJstep[i]=0;
	
		sampler->priority[i] = 1; //Default priority

		sampler->rvec[i] = gsl_rng_alloc(T);

		//Seed differently
		gsl_rng_set(sampler->rvec[i] , i+1);
	
		sampler->check_stepsize_freq[i] = 500;
		//max probability is a function of the temperature -- higher temp 
		//are allowed to step more
		//sampler->max_target_accept_ratio[i] = .90-.15/sampler->chain_temps[i];
		//sampler->min_target_accept_ratio[i] = .60;
		sampler->max_target_accept_ratio[i] = .60-.2/sampler->chain_temps[i];
		sampler->min_target_accept_ratio[i] = .2;
		sampler->gauss_last_accept_ct[i] = 0.;
		sampler->gauss_last_reject_ct[i] = 0.;
		sampler->fish_last_accept_ct[i] = 0.;
		sampler->fish_last_reject_ct[i] = 0.;
		sampler->de_last_accept_ct[i] = 0.;
		sampler->de_last_reject_ct[i] = 0.;
		if(sampler->update_RJ_width){
			sampler->RJstep_last_accept_ct[i] = 0.;
			sampler->RJstep_last_reject_ct[i] = 0.;
		}

		//Initial width size for all chains, all steps is 1.
		sampler->randgauss_width[i][0]=.01;
		sampler->randgauss_width[i][1]=.05;
		sampler->randgauss_width[i][2]=.05;
		sampler->randgauss_width[i][3]=.5;
		//For RJPTMCMC, this may not be used, but it'll be available
		sampler->randgauss_width[i][4]=.5;

	}		
	sampler->history = allocate_3D_array(sampler->chain_N, 
				sampler->history_length, sampler->max_dim);
	sampler->fisher_vecs = allocate_3D_array(sampler->chain_N, 
				sampler->max_dim, sampler->max_dim);
	sampler->fisher_vals = allocate_2D_array(sampler->chain_N, sampler->max_dim);
	

	//Trouble Shooting:
	if(sampler->log_ll || sampler->log_lp){
		sampler->ll_lp_output = allocate_3D_array(sampler->chain_N, 
			sampler->N_steps, 2);
	}
}

void deallocate_sampler_mem(sampler *sampler)
{
	int i;
	for (i =0; i<sampler->chain_N; i++)
	{
		free(sampler->step_prob[i]);
		free(sampler->prob_boundaries[i]); 
		gsl_rng_free(sampler->rvec[i]);

	}		
	free(sampler->step_prob); 
	free(sampler->prob_boundaries); 
	free(sampler->de_primed);
	free(sampler->current_hist_pos);
	free(sampler->fish_accept_ct);
	free(sampler->fish_reject_ct);
	free(sampler->de_accept_ct);
	free(sampler->de_reject_ct);
	free(sampler->gauss_accept_ct);
	free(sampler->gauss_reject_ct);
	free(sampler->mmala_accept_ct);
	free(sampler->mmala_reject_ct);
	free(sampler->RJstep_accept_ct);
	free(sampler->RJstep_reject_ct);
	free(sampler->chain_pos);
	free(sampler->waiting);
	free(sampler->swap_accept_ct);
	free(sampler->swap_reject_ct);
	free(sampler->step_accept_ct);
	free(sampler->step_reject_ct);

	free(sampler->nan_counter);
	free(sampler->num_gauss);
	free(sampler->num_fish);
	free(sampler->num_de);
	free(sampler->num_mmala);
	free(sampler->num_RJstep);

	free(sampler->current_likelihoods);
	free(sampler->priority);

	deallocate_3D_array(sampler->history,sampler->chain_N, 
				sampler->history_length, sampler->max_dim);
	deallocate_3D_array(sampler->fisher_vecs, sampler->chain_N, sampler->max_dim, sampler->max_dim);
	deallocate_2D_array(sampler->fisher_vals, sampler->chain_N, sampler->max_dim);
 
	free(sampler->fisher_update_ct);
	free(sampler->rvec);

	free(sampler->check_stepsize_freq);
	free(sampler->gauss_last_accept_ct);
	free(sampler->gauss_last_reject_ct);
	free(sampler->de_last_accept_ct);
	free(sampler->de_last_reject_ct);
	free(sampler->fish_last_accept_ct);
	free(sampler->fish_last_reject_ct);
	if(sampler->update_RJ_width){
		free(sampler->RJstep_last_accept_ct);
		free(sampler->RJstep_last_reject_ct);
	}
	free(sampler->max_target_accept_ratio);
	free(sampler->min_target_accept_ratio);
	deallocate_2D_array(sampler->randgauss_width,sampler->chain_N, sampler->types_of_steps);
	if(sampler->RJMCMC){
		//deallocate_3D_array(sampler->param_status,sampler->chain_N, sampler->N_steps, sampler->max_dim);
		deallocate_3D_array(sampler->history_status,sampler->chain_N, sampler->history_length, sampler->max_dim);
	}
	else{
		for(int j = 0; j<sampler->chain_N; j++){
			free(sampler->param_status[j][0]);
			free(sampler->param_status[j]);
		}
		free(sampler->param_status);
		for(int j = 0; j<sampler->chain_N; j++){
			free(sampler->history_status[j][0]);
			free(sampler->history_status[j]);
		}
		free(sampler->history_status);
	}
	free(sampler->ref_chain_status);
	

	//Trouble shooting
	if(sampler->log_ll || sampler->log_lp){
		deallocate_3D_array(sampler->ll_lp_output,
			sampler->chain_N, sampler->N_steps, 2);
	}
	
}

void update_history(sampler *sampler, double *new_params, int *new_param_status, int chain_index)
{
	if(sampler->current_hist_pos[chain_index] < sampler->history_length-1)
	{
		sampler->current_hist_pos[chain_index]+=1;
	}
	else
	{
		sampler->current_hist_pos[chain_index] = 0;
	}
	for (int i =0; i < sampler->max_dim; i++)
	{
		sampler->history[chain_index][sampler->current_hist_pos[chain_index]][i] =
			new_params[i];
		sampler->history_status[chain_index][sampler->current_hist_pos[chain_index]][i] =
			new_param_status[i];
	}
		

}


void write_stat_file(sampler *sampler, 
		std::string filename
		)
{
	int rejected_swps=0, accepted_swps=0;
	int *accepted_steps=(int *)malloc(sizeof(int)*sampler->chain_N);
	int *rejected_steps=(int *)malloc(sizeof(int)*sampler->chain_N);
		
	for (int i =0;i<sampler->chain_N; i++)
	{
		accepted_swps+=sampler->swap_accept_ct[i];
		rejected_swps+=sampler->swap_reject_ct[i];
		accepted_steps[i]=sampler->step_accept_ct[i];
		rejected_steps[i]=sampler->step_reject_ct[i];
	}
	double total_swps= accepted_swps + rejected_swps;
	double accepted_swp_fraction = (double)accepted_swps/(total_swps);		
	double rejected_swp_fraction = (double)rejected_swps/(total_swps);		
	
	std::ofstream out_file;
	out_file.open(filename);	
	out_file.precision(4);
	//File variables
	int width= 80;
	int third = (int)((double)width/3.);
	int half = (int)((double)width/2.);
	int fourth = (int)((double)width/4.);
	int fifth = (int)((double)width/5.);
	int sixth = (int)((double)width/6.);
	int seventh = (int)((double)width/7.);
	int eighth = (int)((double)width/8.);
	int sixteenth = (int)((double)width/16.);
	
	//Sampler parameters
	out_file<<std::setw(width)<<std::left<<
		"Parameters of sampler: "<<std::endl;	
	out_file<<
		std::setw(eighth)<<std::left<<
		"Max/Min dimension: "<<
		std::setw(sixteenth)<<std::left<<
		sampler->max_dim<<
		std::setw(sixteenth)<<std::left<<
		sampler->min_dim<<
		std::setw(fourth)<<std::left<<
		"Length of History: "<<
		std::setw(fourth)<<std::left<<
		sampler->history_length<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Number of Chains: "<<
		std::setw(fourth)<<std::left<<
		sampler->chain_N<<
		std::setw(fourth)<<std::left<<
		"Threads/Stochastic: "<<
		std::setw(fourth)<<std::left<<
		sampler->num_threads<<" / "<<sampler->pool<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Steps per chain: "<<
		std::setw(fourth)<<std::left<<
		sampler->N_steps<<
		std::setw(fourth)<<std::left<<
		"Chain steps/chain swap: "<<
		std::setw(fourth)<<std::left<<
		sampler->swp_freq<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"CPU time - sampler (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu/60.<<
		std::setw(fourth)<<std::left<<
		"CPU time - sampler(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Wall time - sampler (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall/60.<<
		std::setw(fourth)<<std::left<<
		"Wall time - sampler(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"CPU time - auto-corr (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu_ac/60.<<
		std::setw(fourth)<<std::left<<
		"CPU time - auto-corr(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_cpu_ac<<
		std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		"Wall time - auto-corr (min): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall_ac/60.<<
		std::setw(fourth)<<std::left<<
		"Wall time - auto-corr(sec): "<<
		std::setw(fourth)<<std::left<<
		sampler->time_elapsed_wall_ac<<
		std::endl;
	out_file<<
		std::setw(width)<<std::left<<
		"Probabilities of steps (Gaussian, DE, MMALA, FISHER, RJ): "<<std::endl;
	out_file<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][0]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][1]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][2]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][3]<<
		std::setw(fifth)<<std::left<<
		sampler->step_prob[0][4]<<
		std::endl;
	out_file<<std::endl;

	out_file<<std::setw(width)<<"Chain temperature, total number of steps per chain, fraction of accepted swaps"<<std::endl;
	out_file<<
		std::setw(fourth)<<"ID" <<
		std::setw(fourth)<<"Temp" <<
		std::setw(fourth)<<"Step Number" <<
		std::setw(fourth)<<"Accepted Swaps" <<
		std::endl;
	int ts;
	int swpt;
	int swpa;
	for(int i =0; i<sampler->chain_N; i++){
		ts = sampler->fish_accept_ct[i]+sampler->fish_reject_ct[i]+
			sampler->de_accept_ct[i]+sampler->de_reject_ct[i]+
			sampler->mmala_accept_ct[i]+sampler->mmala_reject_ct[i]+
			sampler->RJstep_accept_ct[i]+sampler->RJstep_reject_ct[i]+
			sampler->gauss_accept_ct[i]+sampler->gauss_reject_ct[i];
		swpa = sampler->swap_accept_ct[i];
		swpt = sampler->swap_reject_ct[i] + swpa;
		out_file<<std::setw(fourth)<<i<<
			std::setw(fourth)<<sampler->chain_temps[i]<<
			std::setw(fourth)<<ts<<
			std::setw(fourth)<<(double)swpa/swpt<<
			std::endl;
	}
	out_file<<std::endl;
	double total_step_type;
	out_file<<
		std::setw(width)<<std::left<<
		"Fraction of steps per type : "<<std::endl;
	out_file<<
		std::setw(sixth)<<std::left<<
		"Chain Number"<<
		std::setw(sixth)<<std::left<<
		"Gaussian"<<
		std::setw(sixth)<<std::left<<
		"Diff. Ev."<<
		std::setw(sixth)<<std::left<<
		"MMALA"<<
		std::setw(sixth)<<std::left<<
		"Fisher"<<
		std::setw(sixth)<<std::left<<
		"RJ"<<
		std::endl;
	for (int i =0; i < sampler->chain_N; i++){	
	 	total_step_type= sampler->num_gauss[i]+sampler->num_mmala[i]+
				sampler->num_de[i]+sampler->num_fish[i]+sampler->num_RJstep[i];
		out_file<<
			std::setw(sixth)<<std::left<<
			i<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_gauss[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_de[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_mmala[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_fish[i]/total_step_type<<
			std::setw(sixth)<<std::left<<
			(double)sampler->num_RJstep[i]/total_step_type<<
			std::endl;
		
	}
	out_file<<std::endl;	
	//########################################################
	
	out_file<<
		std::setw(width)<<std::left<<
		"Final width of Gaussian random number per step type: "<<std::endl;
	out_file<<
		std::setw(sixth)<<std::left<<
		"Chain Number"<<
		std::setw(sixth)<<std::left<<
		"Gaussian"<<
		std::setw(sixth)<<std::left<<
		"Diff. Ev."<<
		std::setw(sixth)<<std::left<<
		"MMALA"<<
		std::setw(sixth)<<std::left<<
		"Fisher"<<
		std::setw(sixth)<<std::left<<
		"RJ"<<
		std::endl;
	for (int i =0; i < sampler->chain_N; i++){	
		out_file<<
			std::setw(sixth)<<std::left<<
			i<<
			std::setw(sixth)<<std::left<<
			(double)sampler->randgauss_width[i][0]<<
			std::setw(sixth)<<std::left<<
			(double)sampler->randgauss_width[i][1]<<
			std::setw(sixth)<<std::left<<
			(double)sampler->randgauss_width[i][2]<<
			std::setw(sixth)<<std::left<<
			(double)sampler->randgauss_width[i][3]<<
			std::setw(sixth)<<std::left<<
			(double)sampler->randgauss_width[i][4]<<
			std::endl;
		
	}
	out_file<<std::endl;	

	//#######################################################

	double acc_total=0;
	double rej_total=0;

	out_file<<std::left<<"Fraction of accepted steps for each chain temp (by step): "<<std::endl;
	out_file.width(width);
	
	out_file<<
		std::left<<std::setw(seventh)<<"Chain ID"<<
		std::setw(seventh)<<std::left<<"GAUSS"<<
		std::setw(seventh)<<std::left<<"Diff Ev"<<
		std::setw(seventh)<<std::left<<"MMALA"<<
		std::setw(seventh)<<std::left<<"Fisher"<<
		std::setw(seventh)<<std::left<<"RJ"<<
		std::setw(seventh)<<std::left<<"Total"<<
		std::endl;

	double total ;
	double acc_frac = 0;
	double rej_frac = 0;
	double gtotal ;
	double detotal ;
	double mmtotal ;
	double ftotal ;
	double RJtotal ;
	double gtotal_total =0;
	double detotal_total =0;
	double mmtotal_total =0;
	double ftotal_total =0;
	double RJtotal_total =0;
	double total_total =0;
	double gacc_frac = 0;
	double deacc_frac = 0;
	double mmacc_frac = 0;
	double facc_frac = 0;
	double RJacc_frac = 0;
	double gacc_frac_total = 0;
	double deacc_frac_total = 0;
	double mmacc_frac_total = 0;
	double facc_frac_total = 0;
	double RJacc_frac_total = 0;
	double acc_frac_total = 0;
	for (int i =0; i<sampler->chain_N;i++){
		//####################################
		gtotal = sampler->gauss_accept_ct[i]+sampler->gauss_reject_ct[i];
		gacc_frac = (double)sampler->gauss_accept_ct[i]/gtotal;

		detotal = sampler->de_accept_ct[i]+sampler->de_reject_ct[i];
		deacc_frac = (double)sampler->de_accept_ct[i]/detotal;

		mmtotal = sampler->mmala_accept_ct[i]+sampler->mmala_reject_ct[i];
		mmacc_frac = (double)sampler->mmala_accept_ct[i]/mmtotal;
	
		ftotal = sampler->fish_accept_ct[i]+sampler->fish_reject_ct[i];
		facc_frac = (double)sampler->fish_accept_ct[i]/ftotal;

		RJtotal = sampler->RJstep_accept_ct[i]+sampler->RJstep_reject_ct[i];
		RJacc_frac = (double)sampler->RJstep_accept_ct[i]/RJtotal;

		total = accepted_steps[i]+rejected_steps[i];
		acc_frac = (double)accepted_steps[i]/total;
		
		//####################################
		gtotal_total += gtotal;
		gacc_frac_total+=(double)sampler->gauss_accept_ct[i];

		mmtotal_total += mmtotal;
		mmacc_frac_total+=(double)sampler->mmala_accept_ct[i];

		detotal_total += detotal;
		deacc_frac_total+=(double)sampler->de_accept_ct[i];

		ftotal_total += ftotal;
		facc_frac_total+=(double)sampler->fish_accept_ct[i];

		RJtotal_total += RJtotal;
		RJacc_frac_total+=(double)sampler->RJstep_accept_ct[i];

		total_total += total;
		acc_frac_total+=(double)accepted_steps[i];
		//####################################
		out_file<<std::left<<std::setw(seventh)<<i<<
			std::left<<std::setw(seventh)<<gacc_frac<<
			std::left<<std::setw(seventh)<<deacc_frac<<
			std::left<<std::setw(seventh)<<mmacc_frac<<
			std::left<<std::setw(seventh)<<facc_frac<<
			std::left<<std::setw(seventh)<<RJacc_frac<<
			std::left<<std::setw(seventh)<<acc_frac<<
			std::endl;
	}
	out_file<<
		std::left<<std::setw(seventh)<<"TOTAL: "<<
		std::setw(seventh)<<(double)gacc_frac_total/(gtotal_total)<<
		std::setw(seventh)<<(double)deacc_frac_total/(detotal_total)<<
		std::setw(seventh)<<(double)mmacc_frac_total/(mmtotal_total)<<
		std::setw(seventh)<<(double)facc_frac_total/(ftotal_total)<<
		std::setw(seventh)<<(double)RJacc_frac_total/(RJtotal_total)<<
		std::setw(seventh)<<(double)acc_frac_total/(total_total)<<
		std::endl;

	out_file<<std::endl;	
	
	//Accepted rejected swaps
	out_file<<std::left<<"Number of accepted and rejected swaps for all chain temps: "
			<<std::endl;
	out_file<<std::left<<std::setw(fourth)<<"Accepted: "<<std::setw(fourth)<<accepted_swps
			<<std::setw(fourth)<<"Rejected: "<<std::setw(fourth)<< 
			rejected_swps<<std::endl;
	out_file<<std::left<<std::setw(fourth)<<"Accepted fraction: "<<std::setw(fourth)<<
			accepted_swp_fraction
			<<std::setw(fourth)<<"Rejected fraction: "<<std::setw(fourth)<< 
			rejected_swp_fraction<<std::endl;
	out_file.close();
	free(accepted_steps);free(rejected_steps);
}
/*! \brief Routine that writes metadata and final positions of a sampler to a checkpoint file
 *
 */
void write_checkpoint_file(sampler *sampler, std::string filename)
{
	std::ofstream checkfile;
	checkfile.open(filename);
	checkfile.precision(15);
	checkfile<<sampler->min_dim<<" , "<<sampler->max_dim<<" , "<<sampler->chain_N<<std::endl;//Dim, chain_N
	//Chain temps
	checkfile<<sampler->chain_temps[0];
	for(int i =1 ; i<sampler->chain_N; i++){
		checkfile<<" , "<<sampler->chain_temps[i];
	}
	checkfile<<std::endl;
	//step widths
	for(int i =0 ; i<sampler->chain_N; i++){
		checkfile<<sampler->randgauss_width[i][0];
		for(int j =1 ;j<sampler->types_of_steps; j++){
			checkfile<<" , "<<sampler->randgauss_width[i][j];
		}
		checkfile<<std::endl;
	}
	//final position
	for(int i =0 ; i< sampler->chain_N; i++){
		int pos = sampler->chain_pos[i];
		checkfile<<sampler->output[i][pos][0];
		for(int j =1; j<sampler->max_dim;j++){
			checkfile<<" , "<<sampler->output[i][pos][j];
		}
		for(int j =0; j<sampler->max_dim;j++){
			checkfile<<" , "<<sampler->param_status[i][pos][j];
			//std::cout<<sampler->param_status[i][pos][j]<<std::endl;
		}
		checkfile<<std::endl;
	}
	//history -- only adds if all chains are primed
	bool de_primed= true;
	for(int i =0 ; i<sampler->chain_N; i++){
		if(!sampler->de_primed[i])
			de_primed=false	;
	}
	if(de_primed){
		checkfile<<"true"<<std::endl;
		
		for(int i =0 ; i<sampler->chain_N; i++){
			checkfile<<sampler->history[i][0][0];
			for (int j = 1 ; j<sampler->max_dim; j++){
				checkfile<<" , "<<sampler->history[i][0][j];
			}
			for(int k=1; k<sampler->history_length; k++){
				for (int j = 0 ; j<sampler->max_dim; j++){

					checkfile<<" , "<<sampler->history[i][k][j];
				}
			}
			checkfile<<std::endl;
		}
	}
	else{
		std::cout<<"false"<<std::endl;
	}
	checkfile.close();
}


/*! \brief load temperatures from checkpoint file 
 *
 * Assumed the temps array is already allocated in memory for the correct number of chains
 * 	
 * Just a utility routine to read temperatures from checkpoint file
 *
 * It would be easy to read in the chain number and allocate memory in the function, but I prefer to leave allocation/deallocation up to the client
 */
void load_temps_checkpoint_file(std::string check_file, double *temps, int chain_N)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		std::getline(file_in,line);

		//Second Row -- temps
		std::getline(file_in,line);
		std::stringstream lineStreamtemps(line);
		i=0;
		while(std::getline(lineStreamtemps, item, ',')){
			temps[i] = std::stod(item);
			i++;
		}
	}
	file_in.close();
}

/*! \brief load checkpoint file into sampler struct
 *
 * *NOTE* -- allocate_sampler called in function -- MUST deallocate manually
 *
 * *NOTE* -- sampler->chain_temps allocated internally -- MUST free manually
 */
void load_checkpoint_file(std::string check_file,sampler *sampler)
{
	std::fstream file_in;
	file_in.open(check_file, std::ios::in);
	std::string line;
	std::string item;
	int i;
	if(file_in){
		//First row -- dim, chain_N
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');
		sampler->min_dim = std::stod(item);
		std::getline(lineStream, item, ',');
		sampler->max_dim = std::stod(item);
		if(sampler->min_dim == sampler->max_dim){
			sampler->dimension=sampler->min_dim;
			sampler->RJMCMC = false;
		}
		else{
			sampler->RJMCMC=true;
		}
		std::getline(lineStream, item, ',');
		sampler->chain_N = std::stod(item);

		//Second Row -- temps
		std::getline(file_in,line);
		std::stringstream lineStreamtemps(line);
		sampler->chain_temps = (double *)malloc(sizeof(double)*sampler->chain_N);
		i=0;
		while(std::getline(lineStreamtemps, item, ',')){
			sampler->chain_temps[i] = std::stod(item);
			i++;
		}
		
		//###################################
		//Allocate memory, now we have initial parameters
		allocate_sampler_mem(sampler);
		//###################################

		//third row+chain_N -- step widths
		for(int j =0 ;j<sampler->chain_N; j++){
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreamwidths(line);
			while(std::getline(lineStreamwidths, item, ',')){
				sampler->randgauss_width[j][i] = std::stod(item);
				i++;
			}
		}

		//row 3 +chain_N  to 3+2 chain_N-- initial positions
		for(int j =0 ;j<sampler->chain_N; j++){
			i=0;
			std::getline(file_in,line);
			std::stringstream lineStreampos(line);
			//while(i<sampler->max_dim && std::getline(lineStreampos, item, ',')){
			while(i<sampler->max_dim  ){
				std::getline(lineStreampos, item, ',');
				sampler->output[j][0][i] = std::stod(item);
				i++;
			}
			i=0;
			while(std::getline(lineStreampos, item, ',')){
				sampler->param_status[j][0][i] = std::stod(item);
				i++;
			}
		}
		std::getline(file_in,line);
		std::stringstream lineStreamprimed(line);
		std::getline(lineStreamprimed, item, ',');
		if(item =="true"){
			for(int j =0 ;j<sampler->chain_N; j++){
				std::getline(file_in,line);
				std::stringstream lineStreamhist(line);
				i=0;
				while(std::getline(lineStreamhist,item,',')){	
					int step = i/sampler->max_dim ;
					int pos = i%sampler->max_dim;
					sampler->history[j][step][pos] = std::stod(item);	
					i++;
				}
			}
			for(int j =0 ;j<sampler->chain_N; j++)
				sampler->de_primed[j] =true;
		}
		else{
			for(int j =0 ;j<sampler->chain_N; j++)
				sampler->de_primed[j] =false;

		}
	
		//exit(1);
		
	}
	else{std::cout<<"ERROR -- File "<<check_file<<" not found"<<std::endl; exit(1);}
	file_in.close();
}

void assign_ct_p(sampler *sampler, int step, int chain_index)
{
	if(step ==0) sampler->gauss_accept_ct[chain_index]+=1;
	else if(step ==1) sampler->de_accept_ct[chain_index]+=1;
	else if(step ==2) sampler->mmala_accept_ct[chain_index]+=1;
	else if(step ==3) sampler->fish_accept_ct[chain_index]+=1;
	else if(step ==4) sampler->RJstep_accept_ct[chain_index]+=1;
}
void assign_ct_m(sampler *sampler, int step, int chain_index)
{
	if(step ==0) sampler->gauss_reject_ct[chain_index]+=1;
	else if(step ==1) sampler->de_reject_ct[chain_index]+=1;
	else if(step ==2) sampler->mmala_reject_ct[chain_index]+=1;
	else if(step ==3) sampler->fish_reject_ct[chain_index]+=1;
	else if(step ==4) sampler->RJstep_reject_ct[chain_index]+=1;
}

void assign_initial_pos(sampler *samplerptr,double *initial_pos, int *initial_status,double *seeding_var) 
{
	if(!seeding_var){ 
		for (int j=0;j<samplerptr->chain_N;j++){
			samplerptr->de_primed[j]=false;
			for (int i = 0; i<samplerptr->max_dim; i++)
			{
				//Only doing this last loop because there is sometimes ~5 elements 
				//not initialized on the end of the output, which screw up plotting
				for(int l =0; l<samplerptr->N_steps; l++){
					samplerptr->output[j][l][i] = initial_pos[i];
					samplerptr->param_status[j][l][i] = initial_status[i];
				}
				//std::cout<<initial_pos[i]<<" "<<initial_status[i]<<std::endl;
				
			}
			samplerptr->current_likelihoods[j] =
				 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->max_dim, j)/samplerptr->chain_temps[j];
		}
	}
	//Seed non-zero chains normally with variance as specified
	else{
		int attempts = 0;
		int max_attempts = 10;
		double temp_pos[samplerptr->max_dim];
		int temp_status[samplerptr->max_dim];
		for (int j=0;j<samplerptr->chain_N;j++){
			samplerptr->de_primed[j]=false;
			if(j == 0){
				for (int i = 0; i<samplerptr->max_dim; i++)
				{
				//Only doing this last loop because there is sometimes ~5 elements 
				//not initialized on the end of the output, which screw up plotting
					for(int l =0; l<samplerptr->N_steps; l++){
						samplerptr->output[j][l][i] = initial_pos[i];
						samplerptr->param_status[j][l][i] = initial_status[i];
					}
				}
			}
			else{
				do{
					for(int i =0; i<samplerptr->max_dim; i++){
						if(initial_status[i] ==1){
							temp_pos[i] = gsl_ran_gaussian(samplerptr->rvec[j],seeding_var[i]) + initial_pos[i];
							temp_status[i]=1;
						}
						else{
							temp_pos[i]=0;
							temp_status[i]=0;
						}
					}
					attempts+=1;
					//std::cout<<samplerptr->lp(temp_pos, temp_status,samplerptr->max_dim,j)<<std::endl;
				}while(samplerptr->lp(temp_pos, temp_status,samplerptr->max_dim,j) == limit_inf && attempts<max_attempts);
				attempts =0;
				if(samplerptr->lp(temp_pos, temp_status,samplerptr->max_dim,j) != limit_inf ){
					for(int i =0; i<samplerptr->max_dim;i++){
						for(int l =0; l<samplerptr->N_steps; l++){
							samplerptr->output[j][l][i] = temp_pos[i];
							samplerptr->param_status[j][l][i] = temp_status[i];
						}
					}
				}
				else{
					for(int i =0; i<samplerptr->max_dim;i++){
						for(int l =0; l<samplerptr->N_steps; l++){
							samplerptr->output[j][l][i] = initial_pos[i];
							samplerptr->param_status[j][l][i] = initial_status[i];
						}
					}
			
				}
			}
			
			samplerptr->current_likelihoods[j] =
				 samplerptr->ll(samplerptr->output[j][0],samplerptr->param_status[j][0],samplerptr->max_dim, j)/samplerptr->chain_temps[j];
		}
	}
	if(samplerptr->log_ll || samplerptr->log_lp){
		for(int i = 0 ; i<samplerptr->chain_N; i++){
			samplerptr->ll_lp_output[i][0][0] = samplerptr->current_likelihoods[i];
			samplerptr->ll_lp_output[i][0][1] = 
				samplerptr->lp(samplerptr->output[i][0],
				samplerptr->param_status[i][0],
				samplerptr->max_dim,i);
		}
	}
}

/*! \brief Timescale of the PT dynamics
 *
 * kappa in the the language of arXiv:1501.05823v3
 */
double PT_dynamical_timescale(int t0, /**<Timescale of the dyanmics*/
	int nu, /**<Initial amplitude (number of steps to base dynamics on)*/
	int t/**< current time*/
	)
{
	return (1./nu) * (double)(t0) / (t + t0);
}

/*! \brief updates the temperatures for a sampler such that all acceptance rates are equal
 * 
 * Follows the algorithm outlined in arXiv:1501.05823v3
 *
 * Fixed temperatures for the first and last chain
 *
 * used in MCMC_MH_dynamic_PT_alloc_internal
 *
 * For defined results, this should be used while the sampler is using non-pooling methods
 */
void update_temperatures(sampler *samplerptr,
	int t0,
	int nu,
	int t
	)
{
	//acceptance ratios -- as defined by arXiv:1501.05823v3
	//Need to transform, because I save the total acceptance count per chain
	//not the acceptance count between chains
	double *old_temps = new double[samplerptr->chain_N];

	for(int i =0 ; i<samplerptr->chain_N-1; i++){
		old_temps[i] = samplerptr->chain_temps[i];
	}
	double power;
	double kappa = PT_dynamical_timescale(t0, nu, t);
	//std::cout<<kappa<<std::endl;
	for (int i =1 ; i<samplerptr->chain_N-1; i++){
		power = kappa * (samplerptr->A[i] - samplerptr->A[i+1]);	
		samplerptr->chain_temps[i] = samplerptr->chain_temps[i-1] +
			(old_temps[i] - old_temps[i-1]) * std::exp(power);
		//Replace LL, since I store it already scaled by the old temperature
		samplerptr->current_likelihoods[i] =(old_temps[i] / samplerptr->chain_temps[i]) *samplerptr->current_likelihoods[i]; 
	} 
	delete [] old_temps;
}

/*!  \brief For the dynamic PT sampler, this is the function that starts the full sampler with the max number of chains
 *
 * The output file will be reused, but the positions are set back to zero (copying the current position to position zero)
 *
 * Assumes the output, chain_temps have been allocated in memory for the final number of chains chain_N and steps N_steps
 *
 * Allocates memory for the new sampler sampler_new -> it's the user's responsibility to deallocate with deallocate_sampler_mem
 */
void initiate_full_sampler(sampler *sampler_new, sampler *sampler_old, /**<Dynamic sampler*/
	int chain_N_thermo_ensemble, /**< Number of chains used in the thermodynamic ensemble*/
	int chain_N,/**< Number of chains to use in the static sampler*/
	std::string chain_allocation_scheme /**< Scheme to use to allocate any remaining chains*/
	)
{
	//Check to see if more chains need to be allocated
	bool allocate_chains = true;
	if(sampler_old->chain_N == chain_N){	allocate_chains = false;}
	
	//Allocate new sampler 
	sampler_new->chain_N = chain_N;
	sampler_new->ll = sampler_old->ll;
	sampler_new->lp = sampler_old->lp;
	sampler_new->fish = sampler_old->fish;
	sampler_new->N_steps = sampler_old->N_steps;
	sampler_new->dimension = sampler_old->dimension;
	sampler_new->min_dim = sampler_old->min_dim;
	sampler_new->max_dim = sampler_old->max_dim;
	sampler_new->num_threads = sampler_old->num_threads;
	sampler_new->numThreads = sampler_old->numThreads;
	//sampler_new->chain_temps = sampler_old->chain_temps;
	sampler_new->history_length = sampler_old->history_length;
	sampler_new->history_update = sampler_old->history_update;
	sampler_new->fisher_update_number = sampler_old->fisher_update_number;
	sampler_new->pool = sampler_old->pool;
	sampler_new->swp_freq = sampler_old->swp_freq;

	//This is a reference copy, not a value.
	sampler_new->output = sampler_old->output;
	//sampler_new->output = allocate_3D_array(chain_N, sampler_old->N_steps, sampler_old->dimension);

	sampler_new->fisher_exist = sampler_old->fisher_exist;
	sampler_new->progress = sampler_old->progress;

	//TODO
	//Need to populate the temps first, needed for allocation...
	for(int i = 0 ; i<sampler_old->chain_N; i++){
		sampler_new->chain_temps[i] = sampler_old->chain_temps[i];
	}
	for(int i = sampler_old->chain_N ; i<sampler_new->chain_N; i++){
		if(chain_allocation_scheme=="cold"){
			sampler_new->chain_temps[i] = sampler_old->chain_temps[0];
		}
		else if(chain_allocation_scheme=="double"){
			sampler_new->chain_temps[i] = sampler_old->chain_temps[i%sampler_old->chain_N];
		}
		else if(chain_allocation_scheme=="refine"){
			double prev_temp = sampler_old->chain_temps[i%sampler_old->chain_N];
			double next_temp = sampler_old->chain_temps[(i+1)%sampler_old->chain_N];
			sampler_new->chain_temps[i] = std::sqrt(prev_temp*next_temp);
		}
		else if(chain_allocation_scheme=="half_ensemble"){
			sampler_new->chain_temps[i] = 10;	
		}
	}
	allocate_sampler_mem(sampler_new);


	
	//Copy over pertinent sampler data: histories, step_widths, de_primed, current_hist_pos, copy current pos to first pos, current LL
	//Only for chains 0 - chain_N_thermo_ensemble
	for(int i =0; i<sampler_old->chain_N; i++){
		transfer_chain(sampler_new, sampler_old, i, i, false);
	}

	//Allocate extra chains if needed
	//Only for chains chain_N_thermo_ensemble - chain_N
	if(allocate_chains){
		if(chain_allocation_scheme =="cold"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, 0,true );
			}		
		}
		if(chain_allocation_scheme =="double"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, i%sampler_old->chain_N, true);
			}
		}
		//Continuously refine -- needs correction -- right now, its just adding the ``refinement'' temperatures repeatedly
		if(chain_allocation_scheme =="refine"){	
			for(int i =sampler_old->chain_N; i<chain_N; i++){
				transfer_chain(sampler_new, sampler_old, i, i%sampler_old->chain_N,true);
				double prev_temp = sampler_old->chain_temps[i%sampler_old->chain_N];
				double next_temp = sampler_old->chain_temps[(i+1)%sampler_old->chain_N];
				sampler_new->chain_temps[i] = std::sqrt(prev_temp*next_temp);
			}
		}
		if(chain_allocation_scheme =="half_ensemble"){	
			int i = sampler_old->chain_N;
			bool even_odd = true;
			int j = 2;
			while(i	< chain_N){
			//{
				transfer_chain(sampler_new, sampler_old, i, 0, true);
				i++;
				if(i==chain_N) break;
				while(j<sampler_old->chain_N){
					transfer_chain(sampler_new, sampler_old, i, j, true);
					j+=2;
					i++;
					if(i==chain_N) break;
				}
				if (even_odd){
					j = 1;
					even_odd = false;
				}
				else{
					j = 2;
					even_odd = true;
				}
			}
		}
	}


}
/*! \brief Utility to write out the parameters and status of a sampler to a file
 */
void write_output_file(std::string file, int step_num, int max_dimension, double **output, int **status)
{
	std::ofstream out_file;
	out_file.open(file);
	out_file.precision(15);
	for(int i =0; i<step_num; i++){
		for(int j=0; j<max_dimension;j++){
			out_file<<output[i][j]<<" , ";
		}
		for(int j = 0 ; j<max_dimension; j++){
			if(j==max_dimension-1)
				out_file<<status[i][j]<<std::endl;
			else
				out_file<<status[i][j]<<" , ";
		}
	}
	out_file.close();

}
