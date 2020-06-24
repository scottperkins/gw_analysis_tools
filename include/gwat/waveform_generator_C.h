#ifndef WAVEFORM_GENERATOR_C_H
#define WAVEFORM_GENERATOR_C_H

/*! \file 
 * Header file for the C wrapping of the waveform_generation.cpp 
 *
 */
#ifdef __cplusplus
extern "C"
{
#endif
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
		double tc,
		double f_ref,
		double phiRef,
		double *ppE_beta,
		double *ppE_b,
		int Nmod,
		double incl_angle,
		double theta,
		double phi
		);
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
		);
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
		double *ppE_b,
		int Nmod,
		double incl_angle,
		double theta,
		double phi
		);
void initiate_LumD_Z_interp_C();
void free_LumD_Z_interp_C();
#ifdef __cplusplus
}
#endif
#endif 
