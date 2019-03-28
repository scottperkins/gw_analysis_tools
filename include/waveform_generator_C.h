#ifndef WAVEFORM_GENERATOR_C_H
#define WAVEFORM_GENERATOR_C_H

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
		double phic,
		double tc,
		double ppE_beta,
		int ppE_b,
		double incl_angle,
		double theta,
		double phi
		);
#ifdef __cplusplus
}
#endif
#endif 
