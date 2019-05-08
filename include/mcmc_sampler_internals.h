#ifndef MCMC_SAMPLER_INTERNALS_H
#define MCMC_SAMPLER_INTERNALS_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <string>

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
	double **step_prob;
	double **prob_boundaries;
	double *chain_temps;
	double swp_freq;
	int chain_N;
	int N_steps;
	int dimension;
	bool fisher_exist;
	bool *de_primed;
	int history_length;
	int *current_hist_pos;
	double ***history;
	log_prior lp;
	log_likelihood ll;
	fisher fish;
	gsl_rng * r;
	int nan_counter=0;
	int num_gauss =0;
	int num_fish = 0;
	int num_de = 0;
	int num_mmala = 0;
	double time_elapsed_cpu;
	double time_elapsed_wall;
	double time_elapsed_cpu_ac;
	double time_elapsed_wall_ac;
};

int mcmc_step(sampler *sampler, double *current_param,double *next_param, int chain_number);

void gaussian_step(sampler *sampler, double *current_param,double *proposed_param);

void fisher_step(sampler *sampler,double *current_param, double *proposed_param, int chain_index);

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
#endif
