#include "waveform_generator_C.h"
#include "waveform_generator.h"
#include <complex>
#include <string>
#include <math.h>
#include "util.h"

int fourier_waveformC(double *frequencies,
		int length,
		double *waveform_plus_real,
		double *waveform_plus_imag,
		double *waveform_cross_real,
		double *waveform_cross_imag,
		char *generation_method,
		double mass1,
		double mass2,
		double DL,
		double spin1x,
		double spin1y,
		double spin1z,
		double spin2x,
		double spin2y,
		double spin2z,
		double phiRef,
		double tc,
		double f_ref,
		double *ppE_beta,
		int *ppE_b,
		int Nmod,
		double incl_angle,
		double theta,
		double phi
		)
{
	std::string method = generation_method;
	gen_params params;
	params.mass1 = mass1;
	params.mass2 = mass2;
	params.Luminosity_Distance = DL;
	params.spin1[0] = spin1x;
	params.spin1[1] = spin1y;
	params.spin1[2] = spin1z;
	params.spin2[0] = spin2x;
	params.spin2[1] = spin2y;
	params.spin2[2] = spin2z;
	params.tc = tc;
	params.bppe = ppE_b;
	params.Nmod = Nmod;
	params.betappe = ppE_beta;
	params.incl_angle = incl_angle;
	params.theta = theta;
	params.phi = phi;
	params.phiRef = phiRef;
	params.f_ref = f_ref;
	params.NSflag1=false;
	params.NSflag2=false;
	params.sky_average=false;
	std::complex<double> *waveform_plus = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveform_cross = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	fourier_waveform(frequencies, length, waveform_plus, waveform_cross, method,&params);	
	//fourier_waveform(frequencies, length, waveform_plus_real, waveform_plus_imag,waveform_cross_real,waveform_cross_imag, method,&params);	
	//std::cout<<mass1<<" "<<mass2<<" "<<DL<<" "<<spin1z<<" "<<spin2z<<" "<<tc<<" "<<phic<<" "<<ppE_b<<" "<<ppE_beta<<" "<<incl_angle<<" "<<theta<<" "<<phi<<" "<<phiRef<<" "<<f_ref<<" "<<generation_method<<" "<<std::endl;
	for (int i = 0; i<length; i ++){
		waveform_plus_real[i] = real(waveform_plus[i]);
		waveform_plus_imag[i] = imag(waveform_plus[i]);
		waveform_cross_real[i] = real(waveform_cross[i]);
		waveform_cross_imag[i] = imag(waveform_cross[i]);
		//waveform_plus_real[i] =1;
		//waveform_plus_imag[i] = 1;
		//waveform_cross_real[i] = 1;
		//waveform_cross_imag[i] =1;
		//std::cout<<i<<std::endl;
		//std::cout<<waveform_cross_real[i]<<std::endl;
		//std::cout<<waveform_cross_imag[i]<<std::endl;
		//std::cout<<waveform_plus_real[i]<<std::endl;
		//std::cout<<waveform_plus_imag[i]<<std::endl;
	}
	free(waveform_plus);
	free(waveform_cross);
	return 1;
}
int fourier_amplitudeC(double *frequencies,
		int length,
		double *amplitude,
		char *generation_method,
		double mass1,
		double mass2,
		double DL,
		double spin1x,
		double spin1y,
		double spin1z,
		double spin2x,
		double spin2y,
		double spin2z,
		double incl_angle,
		double theta,
		double phi
		)
{
	std::string method = generation_method;
	gen_params params;
	params.mass1 = mass1;
	params.mass2 = mass2;
	params.Luminosity_Distance = DL;
	params.spin1[0] = spin1x;
	params.spin1[1] = spin1y;
	params.spin1[2] = spin1z;
	params.spin2[0] = spin2x;
	params.spin2[1] = spin2y;
	params.spin2[2] = spin2z;
	params.incl_angle = incl_angle;
	params.theta = theta;
	params.phi = phi;
	params.NSflag1=false;
	params.NSflag2=false;
	params.sky_average=false;
	amplitude[0] = 10;
	fourier_amplitude(frequencies, length, amplitude, method,&params);	
	return 1;
}
int fourier_phaseC(double *frequencies,
		int length,
		double *phase,
		char *generation_method,
		double mass1,
		double mass2,
		double DL,
		double spin1x,
		double spin1y,
		double spin1z,
		double spin2x,
		double spin2y,
		double spin2z,
		double tc,
		double f_ref,
		double phiRef,
		double *ppE_beta,
		int *ppE_b,
		int Nmod,
		double incl_angle,
		double theta,
		double phi
		)
{
	std::string method = generation_method;
	gen_params params;
	params.mass1 = mass1;
	params.mass2 = mass2;
	params.spin1[0] = spin1x;
	params.spin1[1] = spin1y;
	params.spin1[2] = spin1z;
	params.spin2[0] = spin2x;
	params.spin2[1] = spin2y;
	params.spin2[2] = spin2z;
	params.tc = tc;
	params.Luminosity_Distance = DL;
	params.phiRef = phiRef;
	params.bppe = ppE_b;
	params.betappe = ppE_beta;
	params.Nmod = Nmod;
	params.incl_angle = incl_angle;
	params.theta = theta;
	params.phi = phi;
	params.f_ref = f_ref;
	params.NSflag1=false;
	params.NSflag2=false;
	params.sky_average=false;
	fourier_phase(frequencies, length, phase, method,&params);	
	return 1;
}

//void initiate_LumD_Z_interp_C()
//{
//	initiate_LumD_Z_interp();
//}
//void free_LumD_Z_interp_C()
//{
//	free_LumD_Z_interp();
//}
