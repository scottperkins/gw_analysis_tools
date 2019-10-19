#ifndef FISHER_H
#define FISHER_H
#include "util.h"

/*! \file 
 */
using namespace std;


void fisher_numerical(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	double **output,
	int dimension, 
	//double *parameters,
	gen_params_base<double> *parameters,
	int order,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL
	);
void calculate_fisher_elements(double *frequency, int length, int dimension, std::complex<double> **response_deriv, double **output, double *psd);

void calculate_fisher_elements_batch(double *frequency, 
	int length, 
	int base_dimension, 
	int full_dimension, 
	std::complex<double> **response_deriv, 
	double **output,
	double *psd);

void calculate_derivatives_old(double  **amplitude_deriv, 
       	double **phase_deriv,
       	double *amplitude,
       	double *frequencies,
       	int length, 
       	string detector, 
       	string  gen_method,
       	gen_params *parameters);

void calculate_derivatives(std::complex<double>  **response_deriv, 
       	double *frequencies,
       	int length, 
       	int dimension, 
       	string detector, 
       	string  gen_method,
       	gen_params_base<double> *parameters, 
	int order);

void fisher_autodiff(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	double **output,
	int dimension, 
	//double *parameters,
	gen_params *parameters,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL
	);
void fisher_autodiff_interp(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	double **output,
	int dimension, 
	//double *parameters,
	gen_params *parameters,
	int downsampling_factor	,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL
	);

void fisher_autodiff_batch_mod(double *frequency, 
	int length,
	string generation_method, 
	string detector, 
	double **output,
	int base_dimension, 
	int full_dimension, 
	gen_params *parameters,
	int *amp_tapes = NULL,
	int *phase_tapes = NULL,
	double *noise = NULL
	);


void PhenomP_fisher(double *frequency,
	int length,
	gen_params *parameters,
	std::complex<double> **waveform_derivative,
	int *amp_tapes,
	int *phase_tapes,
	std::string detector
	);

void construct_waveform_derivative(double *frequency, 
	int length,
	int dimension,
	std::complex<double> **waveform_deriv,
	source_parameters<double> *input_params,
	int *waveform_tapes
	);

void calculate_derivatives_autodiff(double *frequency,
	int length,
	int dimension,
	std::string generation_method,
	gen_params *parameters,
	std::complex<double> **waveform_deriv,
	int *waveform_tapes,/*<< Waveform tapes -- length=6*/
	std::string detector
	);
void time_phase_corrected_derivative_autodiff_full_hess(double **dt, 
	int length, 
	double *frequencies,
	gen_params_base<double> *params, 
	std::string generation_method, 
	int dimension, 
	bool correct_time);
void time_phase_corrected_derivative_autodiff(double **dt, 
	int length, 
	double *frequencies,
	gen_params_base<double> *params, 
	std::string generation_method, 
	int dimension, 
	bool correct_time);
template <class T>
void time_phase_corrected_derivative_numerical(T **dt, 
	int length, 
	T *frequencies,
	gen_params_base<T> *params, 
	std::string generation_method, 
	int dimension, 
	bool correct_time);
void time_phase_corrected_derivative_autodiff_numerical(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time);

std::string local_generation_method(std::string generation_method);


void prep_fisher_calculation(double *parameters, bool *, double *,double*, int,gen_params_base<double> *input_params, std::string generation_method, int dim);

void detect_adjust_parameters( double *freq_boundaries,
	double *grad_freqs, 
	int *boundary_num,
	gen_params_base<double> *input_params, 
	std::string generation_method, 
	std::string detector, 
	int dim);

void unpack_parameters(double *parameters, 
	gen_params_base<double> *input_params, 
	std::string generation_method, 
	int dimension, 
	bool *log_factors);

template<class T>
void repack_parameters(T *avec_parameters, gen_params_base<T> *a_params, std::string generation_method, int dim, gen_params_base<double> *original_params);

template<class T>
void repack_non_parameter_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method);

template<class T>
void deallocate_non_param_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method);
#endif
