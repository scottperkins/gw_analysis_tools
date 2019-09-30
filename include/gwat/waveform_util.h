#ifndef WAVEFORM_UTIL_H
#define WAVEFORM_UTIL_H
#include "waveform_generator.h"
#include "util.h"
#include <string>

/*! \file 
 * Header file for waveform specific utilites
 */

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

template<class T>
int fourier_detector_response(T *frequencies, 
			int length,
			std::complex<T> *hplus, 
			std::complex<T> *hcross, 
			std::complex<T> *detector_response, 
			T theta, 
			T phi, 
			std::string detector);
template<class T>
int fourier_detector_response(T *frequencies, 
			int length,
			std::complex<T> *hplus, 
			std::complex<T> *hcross, 
			std::complex<T> *detector_response, 
			T theta, 
			T phi, 
			T psi, 
			std::string detector);
template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
			int length,
			std::complex<T> *hplus, 
			std::complex<T> *hcross, 
			std::complex<T> *detector_response, 
			T ra, 
			T dec, 
			T psi, 
			double gmst, 
			T *times,
			T LISA_alpha0,
			T LISA_phi0,
			T LISA_thetal,
			T LISA_phil,
			std::string detector);

template<class T>
int fourier_detector_response(T *frequencies, 
			int length,
			std::complex<T> *response, 
			std::string detector,
			std::string generation_method,
			gen_params_base<T> *parameters
			);
template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
			int length,
			std::complex<T> *response, 
			std::string detector,
			std::string generation_method,
			gen_params_base<T> *parameters
			);
template<class T>
int fourier_detector_response_equatorial(T *frequencies, 
			int length,
			std::complex<T> *response, 
			std::string detector,
			std::string generation_method,
			gen_params_base<T> *parameters,
			T *times
			);
template<class T>
void time_phase_corrected(T *times, 
	int length, 
	T *frequencies, 
	gen_params_base<T> *params,
	std::string generation_method,
	bool correct_time
	);
int fourier_detector_amplitude_phase(double *frequencies, 
			int length,
			double *amplitude, 
			double *phase, 
			std::string detector,
			std::string generation_method,
			gen_params *parameters
			);
template<class T>
void map_extrinsic_angles(gen_params_base<T> *params);

#endif

