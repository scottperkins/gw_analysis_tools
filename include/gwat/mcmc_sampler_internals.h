#ifndef MCMC_SAMPLER_INTERNALS_H
#define MCMC_SAMPLER_INTERNALS_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <functional>
#include <limits>
#include <iomanip>

/*! \file
 * Internal functions of the generic MCMC sampler (nothing specific to GW)
 */

//typedef double (*log_prior_thread_safe)(double *param, int dimension);
//typedef double (*log_likelihood_thread_safe)(double *param, int dimension);
//typedef void (*fisher_thread_safe)(double *param, int dimension,double **fisher);

const double limit_inf = -std::numeric_limits<double>::infinity();


class mcmc_data_interface{
public:
	int min_dim;
	int chain_id;
	int max_dim;
	int chain_number;
	double RJ_step_width=1;
	bool burn_phase=false;
	~mcmc_data_interface(){};
};

/*! Class storing everything that defines an instance of the sampler
 */
class sampler
{
public:
	mcmc_data_interface **interfaces;
	bool tune=true;
	int types_of_steps = 5;
	double **step_prob;
	double **prob_boundaries;
	double *chain_temps;
	//########
	double **chain_neighborhoods;
	int *chain_neighbors;
	int **chain_neighborhoods_ids;
	int *chain_neighbors_ids;
	/* Chain radius controls the temperature difference allowed for swap proposals whether the chains are isolated or not*/
	int chain_radius=2;
	/*Restricting the swapping means that swaps can be proposed for chains with similar temperatures, in any of the ensembles*/
	bool restrict_swapping=true;
	/* Isolating the ensemble means the chains can only swap with those inside their ensemble*/
	bool isolate_ensembles=false;
	double swap_rate=1./2.;
	bool burn_phase=false;
	/* If true, the ensembles are only allowed to PT swap with those in their ensemble, but each ensemble is allowed to propose steps with the ensemble-approach. The frequency with which two ensembles are walked forward with this proposal is 1/ensemble_rate*/
	bool ensemble_proposal=true;
	double ensemble_proposal_rate = 0.25;
	double ensemble_prop_a = 1.;
	//########
	bool *waiting;
	int *chain_pos;
	double swp_freq;
	int chain_N;
	int numThreads;
	int N_steps;
	int dimension;
	int min_dim;
	int max_dim;
	bool fisher_exist;
	bool *de_primed;
	int *priority;
	bool *ref_chain_status;
	int *ref_chain_ids;
	int ref_chain_num;

	//Sets the cold chains at higher priority to push them up the queue for
	//swapping and stepping
	bool prioritize_cold_chains = false;

	double ***output;

	bool pool;
	int progress=0;
	bool show_progress;
	int num_threads;

	//NOTE: currently, update the history every step until length is 
	//reached, then the history is updated every 20th step, always only
	//keeping history of length history_length (overwrites the list as 
	//it walks forward when it reaches the end)
	int history_length=500;
	//int history_length=5;
	int history_update=5;
	int *current_hist_pos;
	double ***history;
	int ***history_status;
	double *current_likelihoods;


	int *check_stepsize_freq;
	double *max_target_accept_ratio;
	double *min_target_accept_ratio;
	int *gauss_last_accept_ct;
	int *gauss_last_reject_ct;
	int **gauss_last_accept_ct_per_dim;
	int **gauss_last_reject_ct_per_dim;
	int *de_last_accept_ct;
	int *de_last_reject_ct;
	int *fish_last_accept_ct;
	int *fish_last_reject_ct;
	int *RJstep_last_accept_ct;
	int *RJstep_last_reject_ct;
	int **randgauss_width_number;
	double ***randgauss_width;

	double ***fisher_vecs;
	double **fisher_vals;
	double ***fisher_vecs_prev;
	double **fisher_vals_prev;
	double ***fisher_matrix;
	double ***fisher_matrix_prev;
	int *fisher_update_ct;
	double *prop_MH_factor;
	//Number of steps to take with the fisher before updating the fisher 
	//to a new value 
	//NOTE: if this is too low, detailed balance isn't maintained without 
	//accounting for the changing fisher (doesn't cancel in MH ratio)
	//but if the number is high enough, detailed balance is approximately 
	//kept without calculating second fisher
	int fisher_update_number=200;
	//int fisher_update_number=500000;
	//int fisher_update_number=2;

	//log_prior lp;
	//log_likelihood ll;
	//fisher fish;
	std::function<double(double*,int *, mcmc_data_interface *, void *)> lp;
	std::function<double(double*,int *,mcmc_data_interface *, void *)> ll;
	std::function<void(double*,int *,double**,mcmc_data_interface *, void *)> fish;
	
	void ** user_parameters=NULL;
	bool local_param_allocation=false;
 
	gsl_rng ** rvec;

	int *nan_counter;
	int *num_gauss ;
	int *num_fish ;
	int *num_de ;
	int *num_mmala ;
	int *num_RJstep ;

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
	int **gauss_accept_ct_per_dim;
	int **gauss_reject_ct_per_dim;
	int *mmala_accept_ct;
	int *mmala_reject_ct;
	int *RJstep_accept_ct;
	int *RJstep_reject_ct;

	int *swap_accept_ct;
	int *swap_reject_ct;
	int *step_accept_ct;
	int *step_reject_ct;

	//Pointer for testing -- stores log likelihood and log prior for output
	//Quite a bit of memory (order chain_N * N_steps * 2), so probably not good to run every single time. Just for trouble
	//shooting.
	double ***ll_lp_output;
	bool log_ll=false;
	bool log_lp=false;



	//Parameters for dynamic PT allocation
	int *A;
	bool PT_alloc=false;
	bool linear_swapping=true;

	//RJPTMCMC Parameterts
	int ***param_status;
	bool RJMCMC=false;
	std::function<void(double*,double*,int *,int *,mcmc_data_interface *,void *)> rj;
	bool update_RJ_width=true;
	
};


int mcmc_step(sampler *sampler, double *current_param,double *next_param, int *current_status, int *next_status, int chain_number);

void gaussian_step(sampler *sampler, double *current_param,double *proposed_param, int *current_status, int *proposed_status, int chain_id, int *selected_dimension);

void fisher_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status, int chain_index);

void update_fisher(sampler *sampler, double *current_param, int *param_status,int chain_index);

void mmala_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status,int chain_index);

void diff_ev_step(sampler *sampler,double *current_param, double *proposed_param,int *current_status, int *proposed_status, int chain_id);

void RJ_smooth_history(sampler *sampler, double *current_param, int *current_param_status,int base_history_id, double *eff_history_coord, int *eff_history_status, int chain_id);

void RJ_step(sampler *sampler, 
	double *current_param, 
	double *proposed_param, 
	int *current_status, 
	int *proposed_status, 
	int chain_number);

void chain_swap(sampler *sampler, double ***output, int ***param_status,int step_num,int *swp_accepted, int *swp_rejected);

int single_chain_swap(sampler *sampler, double *chain1, double *chain2,int *chain1_status, int *chain2_status, int T1_index, int T2_index);

void assign_probabilities(sampler *sampler, int chain_index);

void transfer_chain(sampler *samplerptr_dest,sampler *samplerptr_source, int id_destination, int id_source, bool transfer_output);

bool check_sampler_status(sampler *samplerptr);

void update_step_widths(sampler *samplerptr, int chain_id);

void allocate_sampler_mem(sampler *sampler);

void deallocate_sampler_mem(sampler *sampler);

void update_history(sampler *sampler, double *new_params, int *new_param_status,int chain_index);


//void auto_correlation_spectral(double *chain, int length, double *autocorr);
//
//double auto_correlation(double *arr, int length, double tolerance);
//
//double auto_correlation_serial(double *arr, int length);
//
//double auto_correlation_grid_search(double *arr, int length, int box_num=10, int final_length= 50, double target_length=.01);
//
//double auto_correlation_internal(double *arr, int length, int lag, double ave);
//
//void auto_corr_intervals(double *data, int length,double *output, int num_segments,  double accuracy);
//
//void write_auto_corr_file_from_data(std::string autocorr_filename, double **output,int intervals, int dimension, int N_steps);
//
//void write_auto_corr_file_from_data_file(std::string autocorr_filename, std::string output_file,int intervals, int dimension, int N_steps);

//void write_stat_file(sampler *sampler, std::string filename, int *accepted_steps, int *rejected_steps,int accepted_swps, int rejected_swps);
void write_stat_file(sampler *sampler, std::string filename);

void write_checkpoint_file(sampler *sampler, std::string filename);

int chain_number_from_checkpoint_file(std::string check_file, int *chain_N);
int dimension_from_checkpoint_file(std::string check_file, int *min_dimension, int *max_dimension);

void load_checkpoint_file(std::string check_file, sampler *sampler);

void load_temps_checkpoint_file(std::string check_file, double *temps, int chain_N);

void assign_ct_p(sampler *sampler, int step, int chain_index, int gauss_dim);
void assign_ct_m(sampler *sampler, int step, int chain_index, int gauss_dim);

void assign_initial_pos(sampler *samplerptr,double *initial_pos, int *initial_status,double *seeding_var) ;

double PT_dynamical_timescale(int t0, int nu, int t);

void update_temp_neighborhoods(sampler *samplerptr);
void update_temperatures_full_ensemble(sampler *samplerptr,
	int t0,
	int nu,
	int t
	);
void update_temperatures(sampler *samplerptr,
	int t0,
	int nu,
	int t
	);
void initiate_full_sampler(sampler *sampler_new, sampler *sampler_old, 
	int chain_N_thermo_ensemble, 
	int chain_N,
	std::string chain_allocation_scheme,
	std::string checkpoint_file_start
	);
void copy_base_checkpoint_properties(std::string check_file,sampler *samplerptr);
void write_output_file(std::string file, int step_num, int max_dimension, double ***output, int ***status, int chain_N,double *temps,bool RJ);

void reduce_output(int step_num, int max_dimension, double ***output_old, int ***status_old,double **output_new, int **status_new,int chain_N,double *temps,bool RJ);

int count_cold_chains(double *temps, int chain_N);
#endif
