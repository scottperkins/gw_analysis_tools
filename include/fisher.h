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



//##############################################################
//outdated
//void amplitude_derivative(double *frequencies, 	
//			int length, 
//			double **amp_derivative,
//			int dimension,
//			double *params,
//			adouble (*amp_fun)(adouble f, adouble *parameters)
//			);
//
//adouble IMRPhenomD_amp_parameter_tranform(adouble f, 
//					adouble *parameters);
//
//std::complex<adouble> IMRPhenomD_waveform_parameter_tranform(adouble f, 					
//					adouble *parameters);
//void waveform_derivative(double *frequencies, 	
//			int length, 
//			std::complex<double> **waveform_derivative,
//			int dimension,
//			double *params,
//			std::complex<adouble> (*waveform)(adouble f, adouble *parameters)
//			);
//void intialize_tape();

#endif
