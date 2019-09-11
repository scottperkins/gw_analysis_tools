#ifndef FISHER_H
#define FISHER_H
#include "util.h"

/*! \file 
 */
using namespace std;


void fisher(double *frequency, 
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

void calculate_derivatives(double  **amplitude_deriv, 
       	double **phase_deriv,
       	double *amplitude,
       	double *frequencies,
       	int length, 
       	string detector, 
       	string  gen_method,
       	gen_params *parameters);

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

void unpack_parameters(adouble *parameters, gen_params *input_params, std::string generation_method, int dim);
//void repack_parameters(adouble *parameters, gen_params_ad *input_params, std::string generation_method, int dim);

#endif
