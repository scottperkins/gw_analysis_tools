
#ifndef GWATPY_WRAPPING_H
#define GWATPY_WRAPPING_H
#include "util.h"
#include "detector_util.h"
#include "mcmc_gw.h"
#include "mcmc_sampler.h"


#ifdef __cplusplus
extern "C"
{
#endif

double gps_to_GMST_radian_py(double gps);
double t_0PN_py(double f, double chirpmass);
double f_0PN_py(double t, double chirpmass);
double DTOA_DETECTOR_py(double RA, double DEC, double GMST_rad, char *det1, char *det2);

void mcmc_data_interface_destructor_py(mcmc_data_interface *interface);
mcmc_data_interface * mcmc_data_interface_py(
	int min_dim,
	int max_dim,
	int chain_id,
	int nested_model_number,
	int chain_number,
	double RJ_step_width,
	bool burn_phase);

MCMC_modification_struct * MCMC_modification_struct_py( 
	int ppE_Nmod, 
	double *bppe,
	int gIMR_Nmod_phi,
	int *gIMR_phii,
	int gIMR_Nmod_sigma,
	int *gIMR_sigmai,
	int gIMR_Nmod_beta,
	int *gIMR_betai,
	int gIMR_Nmod_alpha,
	int *gIMR_alphai,
	bool NSflag1,
	bool NSflag2
	);

void pack_local_mod_structure_py(
	mcmc_data_interface *interface, 
	double *param, 
	int *status, 
	char * waveform_extended, 
	//void * parameters, 
	MCMC_modification_struct *full_struct,
	MCMC_modification_struct *local_struct);
char * MCMC_prep_params_py(
	double *param, 
	double *temp_params, 
	gen_params_base<double> *gen_params, 
	int dimension, char * generation_method, 
	MCMC_modification_struct *mod_struct,
	bool save_gmst);

void MCMC_modification_struct_py_destructor(MCMC_modification_struct *mod_struct);
void repack_parameters_py(
	double *parameters, 
	gen_params_base<double> *gen_param, 
	char * generation_method, 
	int dim );

int time_detector_response_py(double *times,
	int length, 
	double  *response_real,
	double  *response_imaginary,
	char * detector,
	char *generation_method, 
	gen_params_base<double> *parameters);
int fourier_detector_response_py(double *frequencies,
	int length, 
	double  *response_real,
	double  *response_imaginary,
	char * detector,
	char *generation_method, 
	gen_params_base<double> *parameters);

int time_waveform_full_py(double *times,
	int length, 
	double  *wf_plus_real,
	double  *wf_plus_imaginary,
	double  *wf_cross_real,
	double  *wf_cross_imaginary,
	double  *wf_x_real,
	double  *wf_x_imaginary,
	double  *wf_y_real,
	double  *wf_y_imaginary,
	double  *wf_b_real,
	double  *wf_b_imaginary,
	double  *wf_l_real,
	double  *wf_l_imaginary,
	char *generation_method, 
	gen_params_base<double> *parameters);

int fourier_waveform_full_py(double *frequencies,
	int length, 
	double  *wf_plus_real,
	double  *wf_plus_imaginary,
	double  *wf_cross_real,
	double  *wf_cross_imaginary,
	double  *wf_x_real,
	double  *wf_x_imaginary,
	double  *wf_y_real,
	double  *wf_y_imaginary,
	double  *wf_b_real,
	double  *wf_b_imaginary,
	double  *wf_l_real,
	double  *wf_l_imaginary,
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
