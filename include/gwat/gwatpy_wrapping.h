
#ifndef GWATPY_WRAPPING_H
#define GWATPY_WRAPPING_H
#include "util.h"
#include "detector_util.h"


#ifdef __cplusplus
extern "C"
{
#endif

int calculate_chirpmass_py(double mass1, double mass2,double *out);

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
