
#ifndef GWATPY_WRAPPING_H
#define GWATPY_WRAPPING_H
#include "util.h"
#include "detector_util.h"


#ifdef __cplusplus
extern "C"
{
#endif

int fourier_waveform_py(double *frequencies,
	int length, 
	double  *wf_plus_real,
	double  *wf_plus_imaginary,
	double  *wf_cross_real,
	double  *wf_cross_imaginary,
	char * generation_method, 
	gen_params_base<double> *parameters);

gen_params_base<double>* gen_params_base_py(double mass1, double mass2);
void gen_params_base_py_destructor(gen_params_base<double> *p);

int DL_from_Z_py(double z, char * COSMOLOGY, double *out);
int calculate_chirpmass_py(double mass1, double mass2,double *out);
int calculate_eta_py(double mass1, double mass2,double *out);
int calculate_mass1_py(double chirpmass, double eta,double *out);
int calculate_mass2_py(double chirpmass, double eta,double *out);

int get_detector_parameters(char *detector, double *LAT,double *LON, double *location, double *response_tensor);

void populate_noise_py(double *frequencies, char * detector, double *noise_root, int length, double integration_time);

void ppE_theory_fisher_transformation_py(double m1,
	double m2,
	double *spin1,
	double *spin2,
	double chip,
	double phip,
	double Luminosity_Distance,
	double phiRef,
	double tc,
	double RA,
	double DEC,
	double phi_l,
	double theta_l,
	double psi,
	double incl_angle,
	double gmst,
	bool reduced_spin,
	bool sky_average ,
	bool NSflag1 ,
	bool NSflag2 ,
	bool horizon_coord ,
	bool equatorial_orientation ,
	int Nmod,
	double *betappe,
	double **original_fisher,
	double **new_fisher,
	char * original_method,
	char * new_method,
	int dimension
	);

#ifdef __cplusplus
}
#endif

#endif
