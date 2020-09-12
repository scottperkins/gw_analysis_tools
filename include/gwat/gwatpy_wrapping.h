
#ifndef GWATPY_WRAPPING_H
#define GWATPY_WRAPPING_H
#include "util.h"
#include "detector_util.h"


#ifdef __cplusplus
extern "C"
{
#endif

int fourier_detector_response_py(double *frequencies,
	int length, 
	double  *response_real,
	double  *response_imaginary,
	char * detector,
	char *generation_method, 
	gen_params_base<double> *parameters);

int fourier_waveform_py(double *frequencies,
	int length, 
	double  *wf_plus_real,
	double  *wf_plus_imaginary,
	double  *wf_cross_real,
	double  *wf_cross_imaginary,
	char * generation_method, 
	gen_params_base<double> *parameters);

gen_params_base<double>* gen_params_base_py(
	double mass1, 
	double mass2,
	double *spin1,
	double *spin2,
	double Luminosity_Distance,
	double incl_angle,
	double RA,
	double DEC,
	double psi,
	double gmst,
	double f_ref,
	double theta_l,
	double phi_l,
	double theta,
	double phi,
	char * cosmology,
	bool equatorial_orientation,
	bool horizon_coord,
	bool NSflag1,
	bool NSflag2,
	bool dep_postmerger,
	bool shift_time,
	bool shift_phase,
	bool sky_average,
	double LISA_alpha0,
	double LISA_phi0,
	int Nmod_phi,
	int Nmod_sigma,
	int Nmod_beta,
	int Nmod_alpha,
	int *phii,
	int *sigmai,
	int *betai,
	int *alphai,
	double *delta_phi,
	double *delta_sigma,
	double *delta_beta,
	double *delta_alpha,
	int Nmod,
	double *bppe,
	double *betappe
 );
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
