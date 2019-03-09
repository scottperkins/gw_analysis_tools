#ifndef WAVEFORM_GENERATOR_H
#define WAVEFORM_GENERATOR_H
#include <math.h>
//#include "general_parameter_structures.h"
#include "util.h"
#include <complex>
#include <string>

/*! \file 
 */

int fourier_waveform(double *frequencies, 
			int length,
			std::complex<double> *waveform_plus, 
			std::complex<double> *waveform_cross, 
			std::string generation_method,
			gen_params *parameters
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

int fourier_amplitude(double *frequencies, 
			int length,
			double *amplitude, 
			std::string generation_method,
			gen_params *parameters
			);

int fourier_phase(double *frequencies, 
			int length,
			double *phase, 
			std::string generation_method,
			gen_params *parameters
			);

//int fourier_waveform_polarizations(double *frequencies,
//				int length,
//				std::complex<double> hcross,
//				std::complex<double> hplus,
//				string generation_method,
//				gen_params *parameters
//				)
#endif
