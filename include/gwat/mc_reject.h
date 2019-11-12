#ifndef MC_REJECT_H
#define MC_REJECT_H
#include <functional>
#include <gsl/gsl_rng.h>


class mcr_sampler
{
public:
	//###########################
	bool thread_safe;
	int thread_num;
	int dimension;
	int parallel_job_num;
	double *param_ranges_max;
	double *param_ranges_min;
	std::function<void(double *, void *, int , int)> pfn;
	std::function<void( void *, int , int, gsl_rng *)> draw_param;
	//###########################
	
	//###########################
	double **samples;
	//###########################
	
	mcr_sampler(void (*probability)(double *prob, void *parameters, int d, int threadid), 
		void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r),
		int dimension, 
		double *param_ranges_max,
		double *param_ranges_min);

	mcr_sampler(void (*probability)(double *prob, void *parameters, int d, int threadid), 
		void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r),
		int dimension, 
		double *param_ranges_max, 
		double *param_ranges_min, 
		int thread_num, 
		bool thread_safe);

	~mcr_sampler();
	
	void sample_distribution(int N_samples, void **output);
	void sample(void *output, int threadid);

private:
	const gsl_rng_type *T;
	void draw_param_standard(void *prop_params, int d, int threadid, gsl_rng *r);
	gsl_rng **rvec;
	void init(void (*probability)(double *prob, void *parameters, int d, int threadid), 
		void (*draw_parameters)(void *prop_parameters, int d, int threadid, gsl_rng *r),
		int dimension, 
		double *param_ranges_max, 
		double *param_ranges_min, 
		int thread_num, 
		bool thread_safe);
};
struct mcr_job
{
	mcr_sampler *sampler;
	int starting_index=0;
	int job_length=0;
	void **output;
};
void parallel_sample(int threadid,mcr_job job);

#endif
