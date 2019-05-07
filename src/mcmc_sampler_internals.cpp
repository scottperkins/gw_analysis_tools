#include "mcmc_sampler_internals.h"
#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <eigen3/Eigen/Eigen>
#include "util.h"
#include <limits>

double limit_inf = std::numeric_limits<double>::infinity();
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

	
	if (alpha<sampler->prob_boundaries[chain_number][0])
	{
		gaussian_step(sampler, current_param, proposed_param);
	}
	else if (alpha<sampler->prob_boundaries[chain_number][1])
	{
		diff_ev_step(sampler, current_param, proposed_param, chain_number);
	}
	else if (alpha<sampler->prob_boundaries[chain_number][2])
	{
		mmala_step(sampler, current_param, proposed_param);
	}
	else 
	{
		fisher_step(sampler, current_param, proposed_param, chain_number);
	}
	
	double current_lp = sampler->lp(current_param, sampler->dimension);
	double proposed_lp = sampler->lp(proposed_param, sampler->dimension);
	double MH_ratio;

	if(current_lp == limit_inf || proposed_lp == limit_inf){MH_ratio = 0.;}
	else{
		//Calculate log_likelihood and log prior
		double current_ll = sampler->ll(current_param, sampler->dimension);
		current_ll = (current_ll )/sampler->chain_temps[chain_number];
		double current_lp = sampler->lp(current_param, sampler->dimension);
		double proposed_ll = sampler->ll(proposed_param, sampler->dimension);
		proposed_ll = (proposed_ll )/sampler->chain_temps[chain_number];
		//std::cout<<"mcmc update current pos: "<<current_param[0]<<" mcmc update proposed pos: "<<proposed_param[0]<<std::endl;
		//std::cout<<"mcmc update current pos: "<<current_param[1]<<" mcmc update proposed pos: "<<proposed_param[1]<<std::endl;
		//std::cout<<"mcmc update current ll: "<<current_ll<<" mcmc update proposed ll: "<<proposed_ll<<std::endl;
		//std::cout<<"Chain number: "<<chain_number<<std::endl;

		//Calculate MH ratio
		MH_ratio = std::exp(-current_ll+proposed_ll-current_lp + proposed_lp);
	}

	int i;
	if (MH_ratio>1.)
	{
		for ( i=0;i<sampler->dimension; i ++)
		{
			next_param[i] = proposed_param[i];
		}
		return 1;
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
	//double alpha = .0005;
	for (i=0;i<sampler->dimension;i++){
		proposed_param[i] = gsl_ran_gaussian(sampler->r, alpha)+current_param[i];
	}
}

/*!\brief Fisher informed gaussian step
 */
void fisher_step(sampler *sampler, /**< Sampler struct*/
		double *current_param, /**< current position in parameter space*/
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int chain_index
		)
{
	//Fisher calculation
	double **fisher=(double **)malloc(sizeof(double*)*sampler->dimension);	
	for (int i =0; i<sampler->dimension;i++){
		fisher[i] = (double*)malloc(sizeof(double)*sampler->dimension);
	}
	sampler->fish(current_param, sampler->dimension, fisher);

	//Convert to 1D array for Eigen
	double *oneDfisher=(double *)malloc(sizeof(double)*sampler->dimension*sampler->dimension);
	for (int i =0; i<sampler->dimension;i++){
		for (int j = 0; j<sampler->dimension; j++){
			oneDfisher[sampler->dimension*i+j] = fisher[i][j];///
					//sampler->chain_temps[chain_index];
			//std::cout<<oneDfisher[sampler->dimension*i+j]<<std::endl;
		}
		
	}
	
	//Find eigen vectors and eigen values
	Eigen::Map<Eigen::MatrixXd> m(oneDfisher,sampler->dimension,sampler->dimension);
 	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m);
 	Eigen::MatrixXd eigen_vecs = eigensolver.eigenvectors();
 	Eigen::VectorXd eigen_vals = eigensolver.eigenvalues();
	
	//beta determines direction to step in eigen directions
	int beta = (int)(sampler->dimension)*(gsl_rng_uniform(sampler->r));
	//std::cout<<"Val: "<<eigen_vals(beta)<<std::endl;
	//std::cout<<"Direction: "<<beta<<std::endl;
	
	//##############################################################
	//Numerical matrix inversion can be tricky - catch nans here and replace
	//with gaussian step just to not have the program crash because of 
	//one position in parameter space
	//
	//To see impact of nans, variable is stored in sampler->nan_counter
	//##############################################################
	int nansum = 0;
	for(int i =0; i<sampler->dimension; i++)
		nansum+= std::isnan(eigen_vecs.col(beta)(i));	
	nansum+= std::isnan(eigen_vals(beta));
	if(nansum){
		gaussian_step(sampler,current_param, proposed_param);	
		sampler->nan_counter+=1;
	}
	else{
		double alpha = gsl_ran_gaussian(sampler->r, 1);

		double scaling = 0.0;
		for(int i =0; i< sampler->dimension;i++)
		{
			if(abs(eigen_vals(beta))<100){scaling = 100.;}
			else if(abs(eigen_vals(beta))>1000){scaling = 1000.;}
			else{scaling = sampler->chain_temps[chain_index]*abs(eigen_vals(beta));}
			proposed_param[i] = current_param[i] +
				alpha/sqrt(scaling) *eigen_vecs.col(beta)(i);
		}


	}
	for (int i =0; i<sampler->dimension;i++){
		free(fisher[i]);
	}
	free(fisher);
	free(oneDfisher);
}

/*!\brief MMALA informed step -- Currently not supported
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
		double *proposed_param, /**< [out] Proposed position in parameter space*/
		int chain_id
		)
{
	int i = (int)(sampler->history_length)*(gsl_rng_uniform(sampler->r));
	int j;
	do{
		j=(int)(sampler->history_length)*(gsl_rng_uniform(sampler->r));	
	}while(j==i);
		
	double alpha = 1.;
	double beta = gsl_rng_uniform(sampler->r);
	if(beta<.9)
		alpha=gsl_ran_gaussian(sampler->r,1.);
	for (int k = 0; k<sampler->dimension; k++)
	{
		proposed_param[k] = current_param[k] + alpha*
			(sampler->history[chain_id][i][k]-sampler->history[chain_id][j][k]);
	}
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
	for (int i =0; i < sampler->chain_N-1; i ++)
	{
		int success = single_chain_swap(sampler, output[i][step_num],output[i+1][step_num], i, i+1);
		if(success==1)
			*swp_accepted+=1;
		else
			*swp_rejected+=1;
	
	}
}

/*! \brief subroutine to actually swap two chains
 */
int single_chain_swap(sampler *sampler, /**< sampler structure*/
			double *chain1, /**< parameter position of chain that could be changed*/
			double *chain2, /**< chain that is not swapped, but provides parameters to be swapped by the other chain*/
			int T1_index,	/**<number of chain swappe in chain_temps*/
			int T2_index	/**<number of chain swapper in chain_temps*/
			)
{
	double ll1 =  sampler->ll(chain1, sampler->dimension);
	double ll2 =  sampler->ll(chain2, sampler->dimension);
	double T1 = sampler->chain_temps[T1_index];
	double T2 = sampler->chain_temps[T2_index];
	//std::cout<<"ll for chain 1 (swapping): "<<ll1<<std::endl;
	//std::cout<<"ll for chain 2 (swapping): "<<ll2<<std::endl;
	//std::cout<<"Position chain 1 (swapping): "<<chain1[0]<<std::endl;
	//std::cout<<"Position chain 2 (swapping): "<<chain2[0]<<std::endl;
	double MH_ratio = std::exp(ll1/T2 + ll2/T1 - ll1/T1 - ll2/T2);
	if(MH_ratio>1.)
	{
		double temp[sampler->dimension];
		for(int i =0; i < sampler->dimension;i++)
		{
			temp[i] = chain1[i];
			chain1[i] = chain2[i];
			chain2[i]=temp[i];
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
			double temp[sampler->dimension];
			for(int i =0; i < sampler->dimension;i++)
			{
				temp[i] = chain1[i];
				chain1[i] = chain2[i];
				chain2[i]=temp[i];
			}
			return 1;
		}
			
	}
}

void assign_probabilities(sampler *sampler, int chain_index)
{

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
		sampler->step_prob[chain_index][0]=.3;
		sampler->step_prob[chain_index][1]=0;
		sampler->step_prob[chain_index][2]=.0;
		sampler->step_prob[chain_index][3]=.7;

	}
	//No fisher, but de ready
	else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
	{
		//testing
		//sampler->step_prob[chain_index][0]=0;
		//sampler->step_prob[chain_index][1]=1;
		
		sampler->step_prob[chain_index][0]=.3;
		sampler->step_prob[chain_index][1]=.7;
		sampler->step_prob[chain_index][2]=.0;
		sampler->step_prob[chain_index][3]=.0;

	}
	//all methods available
	else
	{
		//sampler->step_prob[chain_index][0]=.2;
		//sampler->step_prob[chain_index][1]=.3;
		//sampler->step_prob[chain_index][2]=.2;
		//sampler->step_prob[chain_index][3]=.3;
		//Testing
		sampler->step_prob[chain_index][0]=.2;
		sampler->step_prob[chain_index][1]=.4;
		sampler->step_prob[chain_index][2]=.0;
		sampler->step_prob[chain_index][3]=.4;

	}
	//Split probabilities into boundaries for if-else loop
	sampler->prob_boundaries[chain_index][0] = sampler->step_prob[chain_index][0];
	sampler->prob_boundaries[chain_index][1] = sampler->step_prob[chain_index][1]+sampler->prob_boundaries[chain_index][0];
	sampler->prob_boundaries[chain_index][2] = sampler->step_prob[chain_index][2]+sampler->prob_boundaries[chain_index][1];
	sampler->prob_boundaries[chain_index][3] = sampler->step_prob[chain_index][3]+sampler->prob_boundaries[chain_index][2];
}	

void allocate_sampler_mem(sampler *sampler)
{
	int i;
	sampler->step_prob = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->prob_boundaries = (double **)malloc(sizeof(double *) * sampler->chain_N);
	sampler->de_primed = (bool *)malloc(sizeof(bool ) * sampler->chain_N);
	sampler->current_hist_pos = (int *)malloc(sizeof(int ) * sampler->chain_N);
	for (i =0; i<sampler->chain_N; i++)
	{
		sampler->step_prob[i] = (double *)malloc(sizeof(double)*4);
		sampler->prob_boundaries[i] = (double *)malloc(sizeof(double)*4);
		sampler->de_primed[i] = false;
		sampler->current_hist_pos[i] = 0;

	}		
	sampler->history = allocate_3D_array(sampler->chain_N, 
				sampler->history_length, sampler->dimension);
}
void deallocate_sampler_mem(sampler *sampler)
{
	int i;
	for (i =0; i<sampler->chain_N; i++)
	{
		free(sampler->step_prob[i]);
		free(sampler->prob_boundaries[i]); 

	}		
	free(sampler->step_prob); 
	free(sampler->prob_boundaries); 
	free(sampler->de_primed);
	free(sampler->current_hist_pos);
	deallocate_3D_array(sampler->history,sampler->chain_N, 
				sampler->history_length, sampler->dimension);
 
	
}

void update_history(sampler *sampler, double *new_params, int chain_index)
{
	if(sampler->current_hist_pos[chain_index] < sampler->history_length-1)
	{
		sampler->current_hist_pos[chain_index]+=1;
	}
	else
	{
		sampler->current_hist_pos[chain_index] = 0;
	}
	for (int i =0; i < sampler->dimension; i++)
	{
		sampler->history[chain_index][sampler->current_hist_pos[chain_index]][i] =
			new_params[i];
	}
		

}


//Calculate the autocorrelation of a chain - if the chain is >100,000,
//the program will use a box-search method to help with computation time
double auto_correlation(double *arr, int length , double tolerance){
	
	//if the chain is short enough, its easier to just calculate it in serial
	if(length<=100000){return auto_correlation_serial(arr,length);}
	
	double sum =0;
	int k;	
	for(k =0 ; k< length; k++){
		sum+= arr[k];
	}
	double ave = sum/length;
	double gamma_0_sum = 0;

	for (k=0;k<length;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	double step_multiplier = 1. + .01*tolerance;
	double error_tol = .01*tolerance;	/*error tolerance to find stopping point*/
	
	double rho = 1;
	int h= (int)(0.1*length);
	int direction = 1;	/*Variable to track which direction h should change each iteration*/
	double gamma_sum, gamma;
	int counter = 0;

	while(rho>.01+error_tol || rho < .01-error_tol || rho <0.){
		if(counter%1000==0)std::cout<<"Rho: "<<rho<<std::endl;
		
		if(counter%10000 == 0)h = 1.12*h;
		
		/* Pick new h based on direction */
		if (direction > 0&& h< (int)(length/step_multiplier)){
			h = (int)(step_multiplier*h);
		}
		else{
			h = (int)(h/step_multiplier);
		}
	
		/*calculate new Gamma and rho*/
		gamma_sum=0;
	
		for(k=0;k<(length-h);k++){
			gamma_sum += (arr[k+h] - ave)*(arr[k]-ave);
		}

		gamma = gamma_sum/(length-h);
		rho = gamma/gamma_0;
		/*Update direction for next iteration*/
		if(rho - (.01+error_tol)>0){
			direction = 1;
		}
		else if ((0.01 - error_tol) - rho>0){
			direction = -1;
		} 
		counter++;
	}	
	printf("Loops Required %i, rho: %f \n",counter,rho);
	return h;
}	

//Serial version for short chains
double auto_correlation_serial(double *arr, int length  ){

	double sum =0;
	int k;	
	for(k =0 ; k< length; k++){
		sum+= arr[k];
	}
	double ave = sum/length;
	double gamma_0_sum = 0;
	for (k=0;k<length;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	double rho = 1;
	double gamma_sum, gamma;
	int counter = 0;
	int h = 1;
	while(rho>.01){	
		h++;
		gamma_sum=0;
		for(k=0;k<(length-h);k++){
			gamma_sum += (arr[k+h] - ave)*(arr[k]-ave);
		}

		gamma = gamma_sum/(length-h);
		rho = gamma/gamma_0;
	}	
	return h;
}

