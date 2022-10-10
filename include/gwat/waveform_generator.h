#ifndef WAVEFORM_GENERATOR_H
#define WAVEFORM_GENERATOR_H
#include <math.h>
//#include "general_parameter_structures.h"
#include "util.h"
#include <complex>
#include <string>

/*! \file 
 */

template<class T>
class waveform_polarizations
{
public:
	//Tensor modes
	std::complex<T> *hplus= NULL;
	std::complex<T> *hcross= NULL;
	//Vector modes
	std::complex<T> *hx= NULL;
	std::complex<T> *hy= NULL;
	//Scalar modes
	std::complex<T> *hb= NULL;
	std::complex<T> *hl= NULL;
	//Track what polarizations are active 
	//Order: {Tensor plus, Tensor cross, vector x, vector y, scalar breathing, scalar longitudinal}
	bool active_polarizations[6]={true,true,false,false,false,false};

	void allocate_memory(int length);
	void deallocate_memory();
};
template<class T>
void assign_polarizations(std::string generation_method, waveform_polarizations<T> *wp);
bool check_extra_polarizations(std::string generation_method);

template<class T>
int time_waveform(T *times, 
	int length,
	waveform_polarizations<T> *wp,
	std::string generation_method,
	gen_params_base<T> *parameters
	);
template<class  T>
int fourier_waveform(T *frequencies, 
			int length,
			waveform_polarizations<T> *wp,
			std::string generation_method,
			gen_params_base<T> *parameters
			);

int fourier_waveform(double *frequencies, 
			int length,
			double *waveform_plus_real, 
			double *waveform_plus_imag, 
			double *waveform_cross_real, 
			double *waveform_cross_imag,
			double *waveform_x_real, 
			double *waveform_x_imag,
			double *waveform_y_real, 
			double *waveform_y_imag,
			double *waveform_b_real, 
			double *waveform_b_imag,
			double *waveform_l_real, 
			double *waveform_l_imag,
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

//double tidal_error(double tidal_s, double tidal_a, double q);
//adouble tidal_error(adouble tidal_s, adouble tidal_a, adouble q);

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
