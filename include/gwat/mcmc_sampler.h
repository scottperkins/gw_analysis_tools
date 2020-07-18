#ifndef MCMC_SAMPLER_H
#define MCMC_SAMPLER_H
#include <iostream>
#include <functional>
#include <math.h>
/*! \file 
 * Header file for mcmc_sampler
 */
#include "mcmc_sampler_internals.h"
#include "util.h"



struct dump_file_struct{
	std::string filename;
	bool trimmed;
	bool cold_only;
	int *file_trim_lengths=NULL;
};
/*! \brief Class that contains all output information from sampler
 *
 * Destructor takes care of all internal memory allocation
 *
 * It's assumed that the chain distribution is fixed, ie when you append output to the internal output, chain 7 is still chain 7
 *
 * It's assumed (for autocorrelation calculations) that all the cold chains have the same total length. If using GWAT MCMC routines, this is the case. If it's not the case, the ac must be calculated manually.
 */
class mcmc_sampler_output
{
public:
	mcmc_sampler_output( int chain_N, int dim);
	~mcmc_sampler_output();
	void populate_chain_temperatures(double *temperatures);
	void update_cold_chain_list();
	void populate_initial_output(double ***new_output,double ***new_logL_logP,int *chain_positions);
	void append_to_output(double ***new_output,double ***new_logL_logP, int *chain_positions);
	void dealloc_output();
	void dealloc_logL_logP();
	void calc_ac_vals( bool trim);
	void count_indep_samples(bool trim);
	int create_data_dump(bool cold_only,bool trim,std::string filename);
	int append_to_data_dump(std::string filename);
	int write_flat_thin_output(std::string filename, bool use_stored_ac,bool trim);
	void set_trim(int trim);

	int chunk_steps = 1000;
	int chain_number;
	double *chain_temperatures=NULL;
	int *cold_chain_ids=NULL;
	int cold_chain_number;
	double ***output=NULL;
	double ***logL_logP=NULL;
	int *chain_lengths=NULL;
	int dimension;
	int **ac_vals=NULL;
	int cold_chain_number_ac_alloc;
	double target_correlation = 0.01;
	int threads = 4;
	int indep_samples=0;
	int *max_acs=NULL;
	int *trim_lengths=NULL;
private:
	int *file_trim_lengths =NULL;
	bool trimmed_file=false;
	std::vector<dump_file_struct *> dump_files;
	std::vector<std::string> dump_file_names;
};


void mcmc_step_threaded(int j);
void mcmc_swap_threaded(int i, int j);

void continue_RJPTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, int *status, mcmc_data_interface *interface, void * parameters),	
	double (*log_likelihood)(double *param, int *status, mcmc_data_interface *interface, void * parameters),
	void (*fisher)(double *param, int *status, double **fisher, mcmc_data_interface *interface, void * parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, mcmc_data_interface *interface,  void * parameters),
	void **user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal_driver(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension, 	
	int N_steps,	
	int chain_N,
	int max_chain_N_thermo_ensemble,
	int swp_freq,	
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int* ,mcmc_data_interface *,void *)> log_prior,
	std::function<double(double*,int*,mcmc_data_interface *,void *)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *,void *)>fisher,
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void continue_RJPTMCMC_MH_internal(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int N_steps,
	int swp_freq,
	std::function<double(double*,int *, mcmc_data_interface *, void *)> log_prior,
	std::function<double(double*,int*,mcmc_data_interface *,void *)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *,void *)>fisher,
	std::function<void(double*,double*, int*,int*,mcmc_data_interface *, void *)> RJ_proposal,
	void ** user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string end_checkpoint_file
	);
void RJPTMCMC_MH_internal(double ***output, 
	int ***parameter_status, 	
	int max_dimension, 	
	int min_dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	int *initial_status, 	
	double *seeding_var, 	
	double *chain_temps,	
	int swp_freq,	
	std::function<double(double*,int *, mcmc_data_interface *, void *)> log_prior,
	std::function<double(double*,int*,mcmc_data_interface *,void *)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *, void *)>fisher,
	std::function<void(double*,double*, int*,int*,mcmc_data_interface *,void *)> RJ_proposal,
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool update_RJ_width, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void RJPTMCMC_MH(double ***output, 
	int ***parameter_status, 
	int max_dimension, 
	int min_dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	int *initial_status, 	
	double *seeding_var, 	
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, int *status, mcmc_data_interface *, void * parameters),	
	double (*log_likelihood)(double *param, int *status, mcmc_data_interface *, void * parameters),
	void (*fisher)(double *param, int *status, double **fisher, mcmc_data_interface *, void *parameters),
	void(*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, mcmc_data_interface *, void *parameters),
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);

void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated(std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output, 
	int N_steps,	
	int max_chain_N_thermo_ensemble,	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, mcmc_data_interface *interface, void * parameters),	
	double (*log_likelihood)(double *param, mcmc_data_interface *interface, void *parameters),
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface, void *parameters),
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output, 
	int N_steps,	
	int max_chain_N_thermo_ensemble,	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int *,mcmc_data_interface *, void *)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *, void*)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface*, void*)>fisher,
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated(mcmc_sampler_output *sampler_output,
	double **output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, mcmc_data_interface *interface, void * parameters),	
	double (*log_likelihood)(double *param, mcmc_data_interface *interface, void *parameters),
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface, void *parameters),
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_internal(mcmc_sampler_output *sampler_output,
	double **output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int *,mcmc_data_interface *interface, void *)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *interface, void*)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *interface, void*)>fisher,
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_file
	);
void dynamic_temperature_internal(sampler *samplerptr, 
	int N_steps, 
	double nu, 
	int t0,
	int swp_freq, 
	int max_chain_N_thermo_ensemble, 
	bool dynamic_chain_number,
	bool show_prog);
void continue_PTMCMC_MH_dynamic_PT_alloc_internal(std::string checkpoint_file_start,
	double ***output, 
	int N_steps,	
	int max_chain_N_thermo_ensemble,	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int *,mcmc_data_interface *interface ,void *)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *interface , void *)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *interface, void *)>fisher,
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool dynamic_chain_number,
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_file,
	bool burn_phase
	);
void continue_PTMCMC_MH_dynamic_PT_alloc(std::string checkpoint_file_start,
	double ***output, 
	int N_steps,	
	int max_chain_N_thermo_ensemble,	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, mcmc_data_interface *interface, void *parameters),	
	double (*log_likelihood)(double *param,  mcmc_data_interface *interface, void *parameters),
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface , void *parameters),
	void ** user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_file
	);
void PTMCMC_MH_dynamic_PT_alloc_internal(double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	std::function<double(double*,int *,mcmc_data_interface *interface, void *)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *interface, void *)> log_likelihood,
	std::function<void(double*,int*,double**,mcmc_data_interface *interface, void *)>fisher,
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	bool dynamic_chain_number,
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_file,
	bool burn_phase 
	);
void PTMCMC_MH_dynamic_PT_alloc(double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	int max_chain_N_thermo_ensemble,	
	double *initial_pos, 	
	double *seeding_var, 	
	double *chain_temps,
	int swp_freq,	
	int t0,
	int nu,
	std::string chain_distribution_scheme, 
	double (*log_prior)(double *param, mcmc_data_interface *interface , void *parameters),	
	double (*log_likelihood)(double *param, mcmc_data_interface *interface, void *parameters),
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface, void * parameters),
	void **user_parameters,
	int numThreads, 
	bool pool, 
	bool show_prog, 
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_file
	);
void continue_PTMCMC_MH(std::string start_checkpoint_file,
	double ***output,
	int N_steps,
	int swp_freq,
	double (*log_prior)(double *param, mcmc_data_interface *interface, void *parameters),
	double (*log_likelihood)(double *param, mcmc_data_interface *interface , void *parameters),
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface, void *parameters),
	void **user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string end_checkpoint_file,
	bool tune
	);
void PTMCMC_MH_loop(sampler *sampler);

void PTMCMC_MH_step_incremental(sampler *sampler, int increment);

void PTMCMC_MH(	double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	double *seeding_var,
	double *chain_temps,	
	int swp_freq,	
	double (*log_prior)(double *param, mcmc_data_interface *interface , void *parameters),	
	double (*log_likelihood)(double *param, mcmc_data_interface *interface, void *parameters),	
	void (*fisher)(double *param,  double **fisher, mcmc_data_interface *interface, void *parameters),
	void **user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_filename
	);
void PTMCMC_MH_internal(	double ***output, 
	int dimension, 	
	int N_steps,	
	int chain_N,	
	double *initial_pos, 	
	double *seeding_var,
	double *chain_temps,	
	int swp_freq,	
	std::function<double(double*,int *,mcmc_data_interface *, void*)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *,void *)> log_likelihood,
	std::function<void(double*,int *,double**,mcmc_data_interface *,void *)>fisher,
	void **user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string checkpoint_filename,
	bool tune,
	bool burn_phase
	);
void continue_PTMCMC_MH_internal(sampler *samplerptr, 
	std::string start_checkpoint_file,
	double ***output,
	int N_steps,
	int swp_freq,
	std::function<double(double*,int *,mcmc_data_interface *,void*)> log_prior,
	std::function<double(double*,int *,mcmc_data_interface *,void *)> log_likelihood,
	std::function<void(double*,int *,double**,mcmc_data_interface *, void*)>fisher,
	void **user_parameters,
	int numThreads,
	bool pool,
	bool show_prog,
	std::string statistics_filename,
	std::string chain_filename,
	std::string end_checkpoint_file,
	bool tune,
	bool burn_phase
	);
#endif
