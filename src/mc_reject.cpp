#include "mc_reject.h"
#include "util.h"
#include "threadPool.h"
#include <gsl/gsl_rng.h>
#include <functional>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! \file 
 *
 * Monte Carlo Rejection sampling 
 *
 * Distribution should be normalized to 1
 */


/*! \brief Class constructor for Monte Carlo Rejection Sampler
 *
 * Assumes no parallel threading
 *
 * ``probability'' is the function that should return the probability of a given sample
 *
 * ``draw_parameters'' -- NULL, if default should be used. Otherwise, this should propose another parameter set
 *
 * param_ranges_max/min are only used if the default draw_parameter is used, otherwise can be NULL
 */
mcr_sampler::mcr_sampler(void (*probability)(double *prob, /**< [out] Probability*/
		void *parameters, /**<Parameters*/
		int d, /**<Dimension*/
		int threadid /**<Thread ID*/),
	void (*draw_parameters)(void *prop_parameters, /**<[out] Proposed Parameters*/
		int d, /**<Dimension*/
		int threadid, /**< Thread ID*/
		gsl_rng *r /**< GSL random number seed*/), 
	int dimension, /**< Dimension*/
	double *param_ranges_max, /**< Max allowed parameters -- Shape [dimension]*/
	double * param_ranges_min/**< Min allowed parameters -- Shape [dimension]*/
	)
{
	this->init(probability, draw_parameters, dimension, param_ranges_max,param_ranges_min, 1, false);
}
/*! \brief Class constructor for Monte Carlo Rejection Sampler
 *
 * Allows parallel threading 
 *
 * ``probability'' is the function that should return the probability of a given sample
 *
 * ``draw_parameters'' -- NULL, if default should be used. Otherwise, this should propose another parameter set
 *
 * param_ranges_max/min are only used if the default draw_parameter is used, otherwise can be NULL
 */
mcr_sampler::mcr_sampler(void (*probability)(double *prob, /**< [out] Probability*/
		void *parameters, /**<Parameters*/
		int d, /**<Dimension*/
		int threadid /**<Thread ID*/),
	void (*draw_parameters)(void *prop_parameters, /**<[out] Proposed Parameters*/
		int d, /**<Dimension*/
		int threadid, /**< Thread ID*/
		gsl_rng *r /**< GSL random number seed*/), 
	int dimension, /**< Dimension*/
	double *param_ranges_max, /**< Max allowed parameters -- Shape [dimension]*/
	double * param_ranges_min,/**< Min allowed parameters -- Shape [dimension]*/
	int thread_num, /**<Thread number to use*/
	bool thread_safe /**< Bool thread safe -- true if parallel threading should be used*/)
{
	this->init(probability, draw_parameters, dimension, param_ranges_max,param_ranges_min, thread_num, thread_safe);
}
void mcr_sampler::init(void (*probability)(double *prob, 
		void *parameters, 
		int d, 
		int threadid), 
	void (*draw_parameters)(void *prop_parameters, 
		int d, 
		int threadid, 
		gsl_rng *r), 
	int dimension, 
	double *param_ranges_max, 
	double * param_ranges_min,
	int thread_num, 
	bool thread_safe)
{
	this->thread_num = thread_num;
	this->thread_safe=thread_safe;
	this->dimension = dimension;
	this->param_ranges_max = param_ranges_max;
	this->param_ranges_min = param_ranges_min;
	this->pfn = probability;
	if(draw_parameters){
		this->draw_param = draw_parameters;
	}
	else{
		this->draw_param=std::bind(&mcr_sampler::draw_param_standard, this, 
			std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
			std::placeholders::_4);
	}

	gsl_rng_env_setup();
	this->T = gsl_rng_default;
	this->rvec = new gsl_rng*[thread_num];
	for(int i = 0; i<thread_num ; i++){
		this->rvec[i] = gsl_rng_alloc(this->T);
		gsl_rng_set(this->rvec[i],i+1);
	}

}
/*! \brief Destructor function for mcr_sampler
 *
 * Deallocates the gsl_rng memory
 */
mcr_sampler::~mcr_sampler()
{
	for(int i = 0 ; i< this->thread_num; i++){
		gsl_rng_free(this->rvec[i]);
	}
	delete [] rvec;
}


/*! \brief Main driver -- Samples from distribution that has been associated with a sampler object
 *
 */
void mcr_sampler::sample_distribution(int N_samples, /*<< Number of samples to produce*/
	void **output/**< [out] Output array shape [N_samples][dimension] Caste to void ** */
	)
{
	if(!this->thread_safe){
		int i =0 ;
		while(i<N_samples){
			this->sample(output[i],0);						
			i++;
		}
	}
	else{
		int i =0 ;
		this->parallel_job_num = 10*this->thread_num;
		mcr_job *jobs = new mcr_job[this->parallel_job_num];	
		int job_length = N_samples/this->parallel_job_num;
		int j = 0 ;
		//while(i < N_samples){
		while(j < this->parallel_job_num){
			mcr_job temp;
			temp.sampler = this;
			temp.starting_index = i;
			temp.output = output;
			if(j<this->parallel_job_num-1){
				temp.job_length = job_length;
				i+=job_length;
			}
			else{
				temp.job_length = N_samples-i;
				i = N_samples;
			}
			jobs[j] = temp;
			j++;
		}
		threadPool<mcr_job,default_comp<mcr_job>> pool(this->thread_num, parallel_sample);
		for(int k = 0 ; k<this->parallel_job_num; k++){
			pool.enqueue(jobs[k]);
		}
	}
}
/*! \brief Standard option to draw samples randomly 
 *
 * Assumes the parameter vector is of type double * 
 *
 * Uses the ranges provided by the user to randomly draw samples uniformly
 */
void mcr_sampler::draw_param_standard(void *prop_params, /**<[out] Proposed Param*/
	int d, /**< Dimension of parameter space*/
	int threadid, /**< Thread ID for parallel execution*/
	gsl_rng *r /**< GSL random number seed*/
	)
{
	double *tmp_params = (double *) prop_params;
	for(int i = 0 ; i<this->dimension; i++){
		tmp_params[i] = gsl_rng_uniform(r) *
			(this->param_ranges_max[i]-this->param_ranges_min[i]) + 
			this->param_ranges_min[i];
	}

}

/*! \brief Internal routine to return one accepted sample
 *
 * Should not be used directly by user -- couldn't make private though
 */
void mcr_sampler::sample(void *output, /**< [out] Output, accepted parameter set*/
	int threadid /**< Thread id -- not currently used, but included for future proofing*/
	)
{
	bool successful_draw = false;
	double r, prob;
	while(!successful_draw ){
		this->draw_param(output, this->dimension, threadid, this->rvec[threadid]);
		this->pfn(&prob, output, this->dimension, threadid);	
		r = gsl_rng_uniform(this->rvec[threadid]);
		if(r<prob){
			successful_draw=true;
			std::cout<<"Success"<<std::endl;
		}
	}

}
/*! \brief internal routine to draw N samples
 *
 * Should not be called by the user
 *
 * Just made this separate from the class for now for convenience
 */
void parallel_sample(int threadid, /*<<thread id*/
	mcr_job job/**< mcr job -- see mc_reject.h */
	)
{
	for(int i = job.starting_index ; i<job.starting_index+job.job_length; i++){
		job.sampler->sample(job.output[i], threadid);
	}

}

