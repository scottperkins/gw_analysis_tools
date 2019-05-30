#ifndef MCMC_SAMPLER_INTERNALS_H
#define MCMC_SAMPLER_INTERNALS_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <functional>

/*! \file
 * Internal functions of the generic MCMC sampler (nothing specific to GW)
 */

//typedef double (*log_prior_thread_safe)(double *param, int dimension);
//typedef double (*log_likelihood_thread_safe)(double *param, int dimension);
//typedef void (*fisher_thread_safe)(double *param, int dimension,double **fisher);

//typedef double (*log_prior)(double *param, int dimension, int thread_id);
//typedef double (*log_likelihood)(double *param, int dimension, int thread_id);
//typedef void (*fisher)(double *param, int dimension,double **fisher, int thread_id);
//std::function<double(double*,int, int)> log_prior;
//std::function<double(double*,int, int)> log_likelihood;
//std::function<void(double*,int,double**,int)> fisher;
/*! Structure storing everything that defines an instance of the sampler
 */
struct sampler
{
	double **step_prob;
	double **prob_boundaries;
	double *chain_temps;
	bool *waiting;
	int *chain_pos;
	double swp_freq;
	int chain_N;
	int numThreads;
	int N_steps;
	int dimension;
	bool fisher_exist;
	bool *de_primed;
	int *priority;
	double ***output;
	bool pool;
	int progress=0;

	int history_length;
	int history_update;
	int *current_hist_pos;
	double ***history;
	double *current_likelihoods;

	double ***fisher_vecs;
	double **fisher_vals;
	int *fisher_update_ct;
	int fisher_update_number;

	//log_prior lp;
	//log_likelihood ll;
	//fisher fish;
	std::function<double(double*,int, int)> lp;
	std::function<double(double*,int, int)> ll;
	std::function<void(double*,int,double**,int)> fish;
 
	gsl_rng ** rvec;

	int *nan_counter;
	int *num_gauss ;
	int *num_fish ;
	int *num_de ;
	int *num_mmala ;

	double time_elapsed_cpu;
	double time_elapsed_wall;
	double time_elapsed_cpu_ac;
	double time_elapsed_wall_ac;

	int *fish_accept_ct;
	int *fish_reject_ct;
	int *de_accept_ct;
	int *de_reject_ct;
	int *gauss_accept_ct;
	int *gauss_reject_ct;
	int *mmala_accept_ct;
	int *mmala_reject_ct;

	int *swap_accept_ct;
	int *swap_reject_ct;
	int *step_accept_ct;
	int *step_reject_ct;
};

int mcmc_step(sampler *sampler, double *current_param,double *next_param, int chain_number);

void gaussian_step(sampler *sampler, double *current_param,double *proposed_param, int chain_id);

void fisher_step(sampler *sampler,double *current_param, double *proposed_param, int chain_index);

void update_fisher(sampler *sampler, double *current_param, int chain_index);

void mmala_step(sampler *sampler,double *current_param, double *proposed_param);

void diff_ev_step(sampler *sampler,double *current_param, double *proposed_param, int chain_id);

void chain_swap(sampler *sampler, double ***output, int step_num,int *swp_accepted, int *swp_rejected);

int single_chain_swap(sampler *sampler, double *chain1, double *chain2,int T1_index, int T2_index);

void assign_probabilities(sampler *sampler, int chain_index);

void allocate_sampler_mem(sampler *sampler);

void deallocate_sampler_mem(sampler *sampler);

void update_history(sampler *sampler, double *new_params, int chain_index);

double auto_correlation(double *arr, int length, double tolerance);

double auto_correlation_serial(double *arr, int length);
void auto_corr_intervals(double *data, int length,double *output, int num_segments,  double accuracy);
void write_stat_file(sampler *sampler, std::string filename, int *accepted_steps, int *rejected_steps,int accepted_swps, int rejected_swps);

void assign_ct_p(sampler *sampler, int step, int chain_index);
void assign_ct_m(sampler *sampler, int step, int chain_index);
#endif
