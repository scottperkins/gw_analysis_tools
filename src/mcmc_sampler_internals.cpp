#include "mcmc_sampler_internals.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <eigen3/Eigen/Eigen>
#include "util.h"
#include <limits>
#include <iomanip>

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

	int step;	
	if (alpha<sampler->prob_boundaries[chain_number][0])
	{
		gaussian_step(sampler, current_param, proposed_param);
		sampler->num_gauss+=1;
		step = 0;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][1])
	{
		diff_ev_step(sampler, current_param, proposed_param, chain_number);
		sampler->num_de+=1;
		step = 1;
	}
	else if (alpha<sampler->prob_boundaries[chain_number][2])
	{
		mmala_step(sampler, current_param, proposed_param);
		sampler->num_mmala+=1;
		step= 2;
	}
	else 
	{
		fisher_step(sampler, current_param, proposed_param, chain_number);
		step = 3;
	}
	
	double current_lp = sampler->lp(current_param, sampler->dimension);
	double proposed_lp = sampler->lp(proposed_param, sampler->dimension);
	double MH_ratio;
	double power;

	if(current_lp == limit_inf || proposed_lp == limit_inf){
		MH_ratio =-1e20;
		//std::cout<<std::exp(proposed_param[0])/MPC_SEC<<std::endl;		
		//std::cout<<proposed_param[1]<<std::endl;		
		//std::cout<<proposed_param[2]<<std::endl;		
		//std::cout<<std::exp(proposed_param[3])/MSOL_SEC<<std::endl;		
		//std::cout<<proposed_param[4]<<std::endl;		
		//std::cout<<proposed_param[5]<<std::endl;		
		//std::cout<<proposed_param[6]<<std::endl;		
		//std::cout<<"YIKES"<<std::endl;
	}
	else{
		//Calculate log_likelihood and log prior
		double current_ll = sampler->ll(current_param, sampler->dimension);
		current_ll = (current_ll )/sampler->chain_temps[chain_number];
		//double current_lp = sampler->lp(current_param, sampler->dimension);
		double proposed_ll = sampler->ll(proposed_param, sampler->dimension);
		proposed_ll = (proposed_ll )/sampler->chain_temps[chain_number];
		//std::cout<<"mcmc update current pos: "<<current_param[0]<<" mcmc update proposed pos: "<<proposed_param[0]<<std::endl;
		//std::cout<<"mcmc update current pos: "<<current_param[1]<<" mcmc update proposed pos: "<<proposed_param[1]<<std::endl;
		//std::cout<<"mcmc update current ll: "<<current_ll<<" mcmc update proposed ll: "<<proposed_ll<<std::endl;
		//std::cout<<"Chain number: "<<chain_number<<std::endl;

		//Calculate MH ratio
		//power = -current_ll+proposed_ll-current_lp + proposed_lp;
		MH_ratio = -current_ll+proposed_ll-current_lp + proposed_lp;
		//if(power>0.){MH_ratio=1.1;}
		//else{MH_ratio = std::exp(power);}
		//std::cout<<power<<" "<<MH_ratio<<std::endl;
	}

	int i;
	//Random number to determine step acceptance
	double beta = log(gsl_rng_uniform(sampler->r));
	//std::cout<<beta<<" "<<MH_ratio<<std::endl;
	if(MH_ratio< beta){
		for ( i=0;i<sampler->dimension; i ++)
		{
			next_param[i] = current_param[i];
		}
		assign_ct_m(sampler, step,chain_number);
		return -1;
	}	
	else
	{
		for ( i=0;i<sampler->dimension; i ++)
		{
			next_param[i] = proposed_param[i];
		}
		assign_ct_p(sampler, step, chain_number);
		return 1;
	}		
	
	//if (MH_ratio>1.)
	//{
	//	for ( i=0;i<sampler->dimension; i ++)
	//	{
	//		next_param[i] = proposed_param[i];
	//	}
	//	assign_ct_p(sampler, step,chain_number);
	//	return 1;
	//}	
	//else 
	//{
	//	//Random number to determine step acceptance
	//	double beta = gsl_rng_uniform(sampler->r);
	//	if(MH_ratio< beta){
	//		for ( i=0;i<sampler->dimension; i ++)
	//		{
	//			next_param[i] = current_param[i];
	//		}
	//		assign_ct_m(sampler, step,chain_number);
	//		return -1;
	//	}	
	//	else
	//	{
	//		for ( i=0;i<sampler->dimension; i ++)
	//		{
	//			next_param[i] = proposed_param[i];
	//		}
	//		assign_ct_p(sampler, step, chain_number);
	//		return 1;
	//	}		
	//
	//}
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
		proposed_param[i] = gsl_ran_gaussian(sampler->r, 0.05*alpha)+current_param[i];
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
	sampler->num_fish+=1;
	if(sampler->fisher_update_ct[chain_index]==sampler->fisher_update_number)
		update_fisher(sampler, current_param, chain_index);	

	sampler->fisher_update_ct[chain_index] += 1;
	
	//beta determines direction to step in eigen directions
	int beta = (int)((sampler->dimension)*(gsl_rng_uniform(sampler->r)));
	
	double alpha = gsl_ran_gaussian(sampler->r, .1);

	double scaling = 0.0;
	if(abs(sampler->fisher_vals[chain_index][beta])<10){scaling = 10.;}
	else if(abs(sampler->fisher_vals[chain_index][beta])>1000){scaling = 1000.;}

	else{scaling = abs(sampler->fisher_vals[chain_index][beta])/
				sampler->chain_temps[chain_index];}
	for(int i =0; i< sampler->dimension;i++)
	{
		//std::cout<<sampler->fisher_vecs[chain_index][beta][i]<<std::endl;
		proposed_param[i] = current_param[i] +
			alpha/sqrt(scaling) *sampler->fisher_vecs[chain_index][beta][i];
	}
	//std::cout<<alpha/sqrt(scaling)<<std::endl;
	//std::cout<<std::endl;

}


void update_fisher(sampler *sampler, double *current_param, int chain_index)
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
	
	//##############################################################
	//Numerical matrix inversion can be tricky - catch nans here and replace
	//with gaussian step just to not have the program crash because of 
	//one position in parameter space
	//
	//To see impact of nans, variable is stored in sampler->nan_counter
	//##############################################################
	int nansum = 0;
	for(int j = 0 ; j<sampler->dimension; j++){
		for(int i =0; i<sampler->dimension; i++)
			nansum+= std::isnan(eigen_vecs.col(j)(i));	
		nansum+= std::isnan(eigen_vals(j));
	}
	if(!nansum){
		for (int i =0; i < sampler-> dimension; i++)
		{
			for(int j = 0; j<sampler->dimension; j++)
			{
				sampler->fisher_vecs[chain_index][i][j] = eigen_vecs.col(i)(j);
			}
			sampler->fisher_vals[chain_index][i]=eigen_vals[i];
			//std::cout<<sampler->fisher_vals[chain_index][i]<<std::endl;
		}
		sampler->fisher_update_ct[chain_index]=0;
	}
	else{ sampler->fisher_update_ct[chain_index]=sampler->fisher_update_number-1;}

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
	int i = (int)((sampler->history_length-1)*(gsl_rng_uniform(sampler->r)));
	int j;
	do{
		j=(int)((sampler->history_length-1)*(gsl_rng_uniform(sampler->r)));	
	}while(j==i);
		
	double alpha = .1;
	double beta = gsl_rng_uniform(sampler->r);
	if(beta<.9)
		alpha=gsl_ran_gaussian(sampler->r,.5);
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
	double pow = ll1/T2 + ll2/T1 - ll1/T1 - ll2/T2;
	double MH_ratio;
	if(pow>1){MH_ratio = 1.1;}
	else{MH_ratio = std::exp(ll1/T2 + ll2/T1 - ll1/T1 - ll2/T2);}
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
		sampler->step_prob[chain_index][0]=.1;
		sampler->step_prob[chain_index][1]=0;
		sampler->step_prob[chain_index][2]=.0;
		sampler->step_prob[chain_index][3]=.9;

	}
	//No fisher, but de ready
	else if (!sampler->fisher_exist && sampler->de_primed[chain_index])
	{
		//testing
		//sampler->step_prob[chain_index][0]=0;
		//sampler->step_prob[chain_index][1]=1;
		
		sampler->step_prob[chain_index][0]=.2;
		sampler->step_prob[chain_index][1]=.8;
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
		sampler->step_prob[chain_index][0]=.1;
		sampler->step_prob[chain_index][1]=.4;
		sampler->step_prob[chain_index][2]=.0;
		sampler->step_prob[chain_index][3]=.5;

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
	sampler->fisher_update_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->swap_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_accept_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	sampler->step_reject_ct = (int *)malloc(sizeof(int) * sampler->chain_N);
	for (i =0; i<sampler->chain_N; i++)
	{
		sampler->step_prob[i] = (double *)malloc(sizeof(double)*4);
		sampler->prob_boundaries[i] = (double *)malloc(sizeof(double)*4);
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
		sampler->fisher_update_ct[i] = sampler->fisher_update_number;
		sampler->chain_pos[i]=0;
		sampler->waiting[i]=true;
		sampler->swap_accept_ct[i]=0;
		sampler->swap_reject_ct[i]=0;
		sampler->step_accept_ct[i]=0;
		sampler->step_reject_ct[i]=0;

	}		
	sampler->history = allocate_3D_array(sampler->chain_N, 
				sampler->history_length, sampler->dimension);
	sampler->fisher_vecs = allocate_3D_array(sampler->chain_N, 
				sampler->dimension, sampler->dimension);
	sampler->fisher_vals = allocate_2D_array(sampler->chain_N, sampler->dimension);
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
	free(sampler->fish_accept_ct);
	free(sampler->fish_reject_ct);
	free(sampler->de_accept_ct);
	free(sampler->de_reject_ct);
	free(sampler->gauss_accept_ct);
	free(sampler->gauss_reject_ct);
	free(sampler->mmala_accept_ct);
	free(sampler->mmala_reject_ct);
	free(sampler->chain_pos);
	free(sampler->waiting);
	free(sampler->swap_accept_ct);
	free(sampler->swap_reject_ct);
	free(sampler->step_accept_ct);
	free(sampler->step_reject_ct);
	deallocate_3D_array(sampler->history,sampler->chain_N, 
				sampler->history_length, sampler->dimension);
	deallocate_3D_array(sampler->fisher_vecs, sampler->chain_N, sampler->dimension, sampler->dimension);
	deallocate_2D_array(sampler->fisher_vals, sampler->chain_N, sampler->dimension);
 
	free(sampler->fisher_update_ct);
	gsl_rng_free(sampler->r);
	
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

/*! \brief Function that computes the autocorrelation length on an array of data at set intervals to help determine convergence
 * 
 */
void auto_corr_intervals(double *data, /**<Input data */
			int length, /**< length of input data*/
			double *output, /**<[out] array that stores the auto-corr lengths -- array[num_segments]*/
			int num_segments, /**< number of segements to compute the auto-corr length*/
			double accuracy /**< longer chains are computed numerically, this specifies the tolerance*/
			)
{
	double stepsize = (double)length/num_segments;
	int lengths[num_segments];
	for (int i =0; i<num_segments;i++)
		lengths[i]=(int)(stepsize*(1. + i));
	double *temp = (double *)malloc(sizeof(double)*length);		
	for(int l =0; l<num_segments; l++){
		for(int j =0; j< lengths[l]; j++){
			temp[j] = data[j];
		}
		//output[l]=auto_correlation(data,lengths[l], accuracy);
		output[l]=auto_correlation_serial(data,lengths[l]);
	}
	free(temp);
	
}


void write_stat_file(sampler *sampler, 
		std::string filename, 
		int *accepted_steps,
		int *rejected_steps,
		int accepted_swps, 
		int rejected_swps
		)
{
	double total_swps= accepted_swps + rejected_swps;
	double accepted_swp_fraction = (double)accepted_swps/(total_swps);		
	double rejected_swp_fraction = (double)rejected_swps/(total_swps);		
	
	std::ofstream out_file;
	out_file.open(filename);	
	//File variables
	int width= 80;
	int third = (int)((double)width/3.);
	int half = (int)((double)width/2.);
	int fourth = (int)((double)width/4.);
	int fifth = (int)((double)width/5.);
	int sixth = (int)((double)width/6.);
	
	//Sampler parameters
	out_file<<std::setw(width)<<std::left<<
		"Parameters of sampler: "<<std::endl;	
	out_file<<
		std::setw(fourth)<<std::left<<
		"Dimension: "<<
		std::setw(fourth)<<std::left<<
		sampler->dimension<<
		std::setw(fourth)<<std::left<<
		"Length of History: "<<
		std::setw(fourth)<<std::left<<
		sampler->history_length<<
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
		"Probabilities of steps (Gaussian, DE, MMALA, FISHER): "<<std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		sampler->step_prob[0][0]<<
		std::setw(fourth)<<std::left<<
		sampler->step_prob[0][1]<<
		std::setw(fourth)<<std::left<<
		sampler->step_prob[0][2]<<
		std::setw(fourth)<<std::left<<
		sampler->step_prob[0][3]<<
		std::endl;
	out_file<<std::endl;

	out_file<<
		std::setw(width)<<std::left<<
		"Number of steps per type (Gaussian, DE, MMALA, FISHER): "<<std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		sampler->num_gauss<<
		std::setw(fourth)<<std::left<<
		sampler->num_de<<
		std::setw(fourth)<<std::left<<
		sampler->num_mmala<<
		std::setw(fourth)<<std::left<<
		sampler->num_fish<<
		std::endl;
	double total_step_type = sampler->num_gauss+sampler->num_mmala+sampler->num_de+sampler->num_fish;
	out_file<<
		std::setw(width)<<std::left<<
		"Fraction of steps per type (Gaussian, DE, MMALA, FISHER): "<<std::endl;
	out_file<<
		std::setw(fourth)<<std::left<<
		(double)sampler->num_gauss/total_step_type<<
		std::setw(fourth)<<std::left<<
		(double)sampler->num_de/total_step_type<<
		std::setw(fourth)<<std::left<<
		(double)sampler->num_mmala/total_step_type<<
		std::setw(fourth)<<std::left<<
		(double)sampler->num_fish/total_step_type<<
		std::endl;

	out_file<<std::endl;	

	double acc_total=0;
	double rej_total=0;
	//Accepted rejected steps
	//out_file.width(width);
	//out_file<<std::left<<"Number of accepted and rejected steps for each chain temp: "<<std::endl;
	//out_file.width(width);
	//
	//out_file<<std::left<<std::setw(third)<<"Chain Temp"<<std::setw(third)<<"accepted"<<std::left<<std::setw(third)<<"rejected"<<std::endl;
	//for (int i =0; i<sampler->chain_N;i++){
	//	acc_total += accepted_steps[i];
	//	rej_total += rejected_steps[i];
	//	out_file<<std::left<<std::setw(third)<<sampler->chain_temps[i]<<
	//		std::setw(third)<<accepted_steps[i]<<
	//		std::left<<std::setw(third)<<rejected_steps[i]<<std::endl;
	//}
	//out_file<<
	//	std::left<<std::setw(third)<<"TOTAL: "<<
	//	std::setw(third)<<acc_total<<
	//	std::left<<std::setw(third)<<rej_total<<
	//	std::endl;


	//out_file<<std::endl;	

	//out_file<<std::left<<"Fraction of accepted and rejected steps for each chain temp: "<<std::endl;
	//out_file.width(width);
	//
	//out_file<<std::left<<std::setw(third)<<"Chain Temp"<<std::setw(third)<<"accepted"<<std::left<<std::setw(third)<<"rejected"<<std::endl;
	//double total ;
	//double acc_frac = 0;
	//double rej_frac = 0;
	//for (int i =0; i<sampler->chain_N;i++){
	//	total = accepted_steps[i]+rejected_steps[i];
	//	acc_frac = (double)accepted_steps[i]/total;
	//	rej_frac = (double)rejected_steps[i]/total;
	//	out_file<<std::left<<std::setw(third)<<sampler->chain_temps[i]<<
	//		std::setw(third)<<acc_frac<<
	//		std::left<<std::setw(third)<<rej_frac<<std::endl;
	//}
	//out_file<<
	//	std::left<<std::setw(third)<<"TOTAL: "<<
	//	std::setw(third)<<(double)acc_total/(acc_total+rej_total)<<
	//	std::left<<std::setw(third)<<(double)rej_total/(acc_total+rej_total)<<
	//	std::endl;

	//out_file<<std::endl;	

	out_file<<std::left<<"Fraction of accepted steps for each chain temp (by step): "<<std::endl;
	out_file.width(width);
	
	out_file<<
		std::left<<std::setw(sixth)<<"Chain Temp"<<
		std::setw(sixth)<<std::left<<"GAUSS"<<
		std::setw(sixth)<<std::left<<"Diff Ev"<<
		std::setw(sixth)<<std::left<<"MMALA"<<
		std::setw(sixth)<<std::left<<"Fisher"<<
		std::setw(sixth)<<std::left<<"Total"<<
		std::endl;

	double total ;
	double acc_frac = 0;
	double rej_frac = 0;
	double gtotal ;
	double detotal ;
	double mmtotal ;
	double ftotal ;
	double gtotal_total =0;
	double detotal_total =0;
	double mmtotal_total =0;
	double ftotal_total =0;
	double total_total =0;
	double gacc_frac = 0;
	double deacc_frac = 0;
	double mmacc_frac = 0;
	double facc_frac = 0;
	double gacc_frac_total = 0;
	double deacc_frac_total = 0;
	double mmacc_frac_total = 0;
	double facc_frac_total = 0;
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

		total_total += total;
		acc_frac_total+=(double)accepted_steps[i];
		//####################################
		out_file<<std::left<<std::setw(sixth)<<sampler->chain_temps[i]<<
			std::left<<std::setw(sixth)<<gacc_frac<<
			std::left<<std::setw(sixth)<<deacc_frac<<
			std::left<<std::setw(sixth)<<mmacc_frac<<
			std::left<<std::setw(sixth)<<facc_frac<<
			std::left<<std::setw(sixth)<<acc_frac<<
			std::endl;
	}
	out_file<<
		std::left<<std::setw(sixth)<<"TOTAL: "<<
		std::setw(sixth)<<(double)gacc_frac_total/(gtotal_total)<<
		std::setw(sixth)<<(double)deacc_frac_total/(detotal_total)<<
		std::setw(sixth)<<(double)mmacc_frac_total/(mmtotal_total)<<
		std::setw(sixth)<<(double)facc_frac_total/(ftotal_total)<<
		std::setw(sixth)<<(double)acc_frac_total/(total_total)<<
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
}

void assign_ct_p(sampler *sampler, int step, int chain_index)
{
	if(step ==0) sampler->gauss_accept_ct[chain_index]+=1;
	else if(step ==1) sampler->de_accept_ct[chain_index]+=1;
	else if(step ==2) sampler->mmala_accept_ct[chain_index]+=1;
	else if(step ==3) sampler->fish_accept_ct[chain_index]+=1;
}
void assign_ct_m(sampler *sampler, int step, int chain_index)
{
	if(step ==0) sampler->gauss_reject_ct[chain_index]+=1;
	else if(step ==1) sampler->de_reject_ct[chain_index]+=1;
	else if(step ==2) sampler->mmala_reject_ct[chain_index]+=1;
	else if(step ==3) sampler->fish_reject_ct[chain_index]+=1;
}
