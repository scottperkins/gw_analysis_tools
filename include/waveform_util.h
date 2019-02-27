#ifndef WAVEFORM_UTIL_H
#define WAVEFORM_UTIL_H
#include "waveform_generator.h"
#include "util.h"


std::complex<double> Q(double theta, double phi, double iota);

double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				std::complex<double> *data,
				std::string detector,
				std::string generation_method,
				gen_params param
				);
double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				double *data_real,
				double *data_imag,
				std::string detector,
				std::string generation_method,
				gen_params param
				);
double calculate_snr(std::string detector,
                        std::complex<double> *waveform,
                        double *frequencies,
                        int length);
#endif

