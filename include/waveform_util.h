#ifndef WAVEFORM_UTIL_H
#define WAVEFORM_UTIL_H
#include "waveform_generator.h"
#include "util.h"
#include <string>



double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				std::complex<double> *data,
				double *psd,
				std::string detector,
				std::string generation_method,
				gen_params *param
				);
double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				double *data_real,
				double *data_imag,
				double *psd,
				std::string detector,
				std::string generation_method,
				gen_params *param
				);
double calculate_snr(std::string detector,
                        std::complex<double> *waveform,
                        double *frequencies,
                        int length);

int fourier_detector_response(double *frequencies, 
			int length,
			std::complex<double> *hplus, 
			std::complex<double> *hcross, 
			std::complex<double> *detector_response, 
			double theta, 
			double phi, 
			std::string detector);

int fourier_detector_response(double *frequencies, 
			int length,
			std::complex<double> *response, 
			std::string detector,
			std::string generation_method,
			gen_params *parameters
			);

#endif

