#include "gwatpy_wrapping.h"
#include "util.h"
#include "fisher.h"

int DL_from_Z_py(double z, char * COSMOLOGY, double *out)
{
	*out = DL_from_Z(z,std::string(COSMOLOGY));
	return 0;
}
int calculate_chirpmass_py(double mass1, double mass2,double *out)
{
	*out = calculate_chirpmass(mass1,mass2);
	return 0;
}

int calculate_eta_py(double mass1, double mass2,double *out)
{
	*out = calculate_eta(mass1,mass2);
	return 0;
}
int calculate_mass1_py(double chirpmass, double eta,double *out)
{
	*out = calculate_mass1(chirpmass,eta);
	return 0;
}
int calculate_mass2_py(double chirpmass, double eta,double *out)
{
	*out = calculate_mass2(chirpmass,eta);
	return 0;
}
void populate_noise_py(double *frequencies, char * detector, double *noise_root, int length, double integration_time){
	populate_noise(frequencies, std::string(detector), noise_root, length, integration_time);
	return;
}

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
	)
{
	gen_params params;
	params.mass1 = m1;
	params.mass2 = m2;
	if(!reduced_spin){
		params.spin1[0] = spin1[0];
		params.spin1[1] = spin1[1];
		params.spin1[2] = spin1[2];
		params.spin2[0] = spin2[0];
		params.spin2[1] = spin2[1];
		params.spin2[2] = spin2[2];
	}
	else{
		params.chip=chip;
		params.chip=phip;
		params.spin1[2]=spin1[2];
		params.spin2[2]=spin2[2];
	}
	params.Luminosity_Distance = Luminosity_Distance;
	params.phiRef = phiRef;
	params.tc = tc;
	params.sky_average = sky_average;
	params.NSflag1 = NSflag1;
	params.NSflag2 = NSflag2;
	params.horizon_coord = horizon_coord;
	params.equatorial_orientation = equatorial_orientation;
	params.RA = RA;
	params.DEC = DEC;
	params.phi_l = phi_l;
	params.theta_l = theta_l;
	params.psi = psi;
	params.incl_angle = incl_angle;
	params.gmst = gmst;
	params.Nmod = Nmod;
	params.betappe = new double[Nmod];
		
	for(int i = 0 ; i<Nmod; i++){
		params.betappe[i]=betappe[i];
	}
	
	ppE_theory_covariance_transformation(std::string(original_method),std::string(new_method),dimension, &params, original_fisher, new_fisher);
	delete [] params.betappe;
}
