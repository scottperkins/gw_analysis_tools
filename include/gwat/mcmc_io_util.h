#ifndef MCMC_IO_UTIL_H
#define MCMC_IO_UTIL_H

#include <iostream>
#include <vector>
#include <functional>
#include <math.h>
/*! \file 
 * Header file for mcmc_io_util
 */

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
	mcmc_sampler_output( int chain_N, int dim,int nested_model_N=0);
	~mcmc_sampler_output();
	void populate_chain_temperatures(double *temperatures);
	void update_cold_chain_list();

	void populate_initial_output(
		double ***new_output,	
		int ***new_status,
		int **new_model_status,
		double ***new_logL_logP,
		int *chain_positions);
	void append_to_output(
		double ***new_output,
		int ***new_status,
		int **new_model_status,
		double ***new_logL_logP, 
		int *chain_positions);

	void calc_ac_vals( bool trim);
	void count_indep_samples(bool trim);

	int create_data_dump(bool cold_only,bool trim,std::string filename);
	int append_to_data_dump(std::string filename);
	int write_flat_thin_output(std::string filename, bool use_stored_ac,bool trim);
	void set_trim(int trim);

	void append_integrated_likelihoods(double *integrated_likelihoods, int * integrated_likelihoods_terms, int ensemble_size_new);
	void calculate_evidence();


	void dealloc_output();
	void dealloc_status();
	void dealloc_model_status();
	void dealloc_logL_logP();
	void dealloc_integrated_likelihoods();

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
	int *integrated_likelihoods_terms=NULL;
	double *integrated_likelihoods=NULL;
	double evidence=0;
	double evidence_error=0;
	int ensemble_size;
	bool calculated_evidence=false;

	bool RJ=false;
	int ***status=NULL;
	int **model_status=NULL;
	int nested_model_number = 0;
private:
	int *file_trim_lengths =NULL;
	bool trimmed_file=false;
	std::vector<dump_file_struct *> dump_files;
	std::vector<std::string> dump_file_names;
};


#endif
