#ifndef WAVEFORM_GENERATOR_H
#define WAVEFORM_GENERATOR_H
#include <math.h>
//#include "general_parameter_structures.h"
#include "util.h"
#include <complex>
#include <string>

/*! \file 
 */

template<class  T>
int fourier_waveform(T *frequencies, 
			int length,
			std::complex<T> *waveform_plus, 
			std::complex<T> *waveform_cross, 
			std::string generation_method,
			gen_params_base<T> *parameters
			);

int fourier_waveform(double *frequencies, 
			int length,
			double *waveform_plus_real, 
			double *waveform_plus_imag, 
			double *waveform_cross_real, 
			double *waveform_cross_imag,
			std::string generation_method,
			gen_params *parameters
			);


int fourier_waveform(double *frequencies, 
			int length,
			std::complex<double> *waveform,
			std::string generation_method,
			gen_params *parameters);

int fourier_waveform(double *frequencies, 
			int length,
			double *waveform_real,
			double *waveform_imag,
			std::string generation_method,
			gen_params *parameters);

template<class T>
int fourier_amplitude(T *frequencies, 
			int length,
			T *amplitude, 
			std::string generation_method,
			gen_params_base<T> *parameters
			);

template<class T>
int fourier_phase(T *frequencies, 
			int length,
			T *phase, 
			std::string generation_method,
			gen_params_base<T> *parameters
			);
template<class T>
int fourier_phase(T *frequencies, 
			int length,
			T *phase_plus, 
			T *phase_cross, 
			std::string generation_method,
			gen_params_base<T> *parameters
			);

template<class T>
std::string prep_source_parameters(source_parameters<T> *out, gen_params_base<T> *in,std::string generation_method);
template<class T>
void cleanup_source_parameters(source_parameters<T> *params,std::string generation_method);
//int fourier_waveform_polarizations(double *frequencies,
//				int length,
//				std::complex<double> hcross,
//				std::complex<double> hplus,
//				string generation_method,
//				gen_params *parameters
//				)
#endif
