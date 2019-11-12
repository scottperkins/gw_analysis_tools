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
 */
mcr_sampler::mcr_sampler(void (*probability)(double *prob, void *parameters, int d, int threadid), 
	void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r), 
	int dimension, 
	double *param_ranges_max, 
	double *param_ranges_min)
{
	this->init(probability, draw_parameters, dimension, param_ranges_max,param_ranges_min, 1, false);
}
/*! \brief Class constructor for Monte Carlo Rejection Sampler
 *
 * Allows parallel threading 
 */
mcr_sampler::mcr_sampler(void (*probability)(double *prob, void *parameters, int d, int threadid), 
	void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r), 
	int dimension, 
	double *param_ranges_max, 
	double * param_ranges_min,
	int thread_num, 
	bool thread_safe)
{
	this->init(probability, draw_parameters, dimension, param_ranges_max,param_ranges_min, thread_num, thread_safe);
}
void mcr_sampler::init(void (*probability)(double *prob, void *parameters, int d, int threadid), 
	void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r), 
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
	//this->pfn = [&probability](double *prob_out, void * param, int d, int tid){ probability(prob_out, param, d, tid);};
	this->pfn = probability;
	if(draw_parameters){
		this->draw_param = [&draw_parameters](void * prop_param, int d, int tid, gsl_rng * r){ draw_parameters( prop_param, d, tid, r);};
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


void mcr_sampler::sample_distribution(int N_samples, 
	void **output)
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
void mcr_sampler::draw_param_standard(void *prop_params, int d, int threadid, gsl_rng *r)
{
	double *tmp_params = (double *) prop_params;
	for(int i = 0 ; i<this->dimension; i++){
		tmp_params[i] = gsl_rng_uniform(r) *
			(this->param_ranges_max[i]-this->param_ranges_min[i]) + 
			this->param_ranges_min[i];
	}

}

void mcr_sampler::sample(void *output, int threadid)
{
	bool successful_draw = false;
	double r, prob;
	while(!successful_draw ){
		this->draw_param(output, this->dimension, threadid, this->rvec[threadid]);
		this->pfn(&prob, output, this->dimension, threadid);	
		r = gsl_rng_uniform(this->rvec[threadid]);
		if(r<prob){
			successful_draw=true;
		}
	}

}
void parallel_sample(int threadid, mcr_job job)
{
	for(int i = job.starting_index ; i<job.starting_index+job.job_length; i++){
		job.sampler->sample(job.output[i], threadid);
	}

}

