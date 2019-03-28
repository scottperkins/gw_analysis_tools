#include "waveform_generator_C.h"
#include "waveform_generator.h"

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
		)
{
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
	params.phic = phic;
	params.bppe = ppE_b;
	params.betappe = ppE_beta;
	params.incl_angle = incl_angle;
	params.theta = theta;
	params.phi = phi;
	
	return 1;
}
