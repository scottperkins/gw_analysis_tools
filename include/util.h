#ifndef UTIL_H
#define UTIL_H
//#include "general_parameter_structures.h"
#include <string>
#include <complex>
#include <adolc/adouble.h>
//#define adouble double
/*! \file 
 *General utilities (functions and structures) independent of modelling method
 */

/*! Euler number*/
const double gamma_E = 0.5772156649015328606065120900824024310421;
/*!Speed of light m/s*/
const double c = 299792458.;
/*!Gravitational constant in m**3/(s**2 SolMass)*/
const double G =6.674e-11*(1.98855e30);
/*! G/c**3 seconds per solar mass*/
const double MSOL_SEC =492549095.e-14; 
/*!consts.kpc.to('m')*1000/c Mpc in sec*/
const double MPC_SEC = 3085677581.e13/c; 

//GSL versions of the constants, but G seems off..
//const double gamma_E = M_EULER;
//const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
//const double G =GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;
//const double MSOL_SEC = GSL_CONST_MKSA_SOLAR_MASS*(GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT/(c*c*c));
//const double MPC_SEC = GSL_CONST_MKSA_PARSEC*1e6/c; 

/*!\struct
 * \brief Structure for interfacing with the libraries 
 * 
 * Structure to interface with the libraries - Units are in solar masses and mpc 
 *  contains the generation parameters including source parameters and theory parameters
 */
struct gen_params
{	
	/*!mass of the larger body in Solar Masses*/
	double mass1;
	/*!mass of the smaller body in Solar Masses*/
	double mass2;
	/*!Luminosity distance to the source*/
	double Luminosity_Distance;
	/*!Spin vector of the larger mass [Sx,Sy,Sz]*/
	double spin1[3];
	/*!Spin vector of the smaller mass [Sx,Sy,Sz]*/
	double spin2[3];
	/*!coalescence phase of the binary*/
	double phic;
	/*!coalescence time of the binary*/
	double tc;
	/*!ppE b parameter (power of the frequency)*/
	int bppe;
	/*!ppE coefficient for the phase modification*/
	double betappe;
	/*!*angle between angular momentum and the total momentum */
	double incl_angle;
	/*! spherical angles for the source location relative to the detector*/
	double theta;
	double phi;
	/*! BOOL flag for early termination of NS binaries*/
	bool NSflag;
};

/*!\brief To speed up calculations within the for loops, we pre-calculate reoccuring powers of M*F and Pi, since the pow() function is prohibatively slow
 *
 * Powers of PI are initialized once, and powers of MF need to be calculated once per for loop (if in the inspiral portion).
 *
 * use the functions precalc_powers_ins_amp, precalc_powers_ins_phase, precalc_powers_pi to initialize
 */
template <class T>
struct useful_powers
{
	T MFthird;
	T MF2third;
	T MF4third;
	T MF5third;
	T MFsquare;
	T MF7third;
	T MF8third;
	T MFcube;
	T MFminus_5third;
	T MF3fourth;
	double PIsquare;
	double PIcube;
	double PIthird;
	double PI2third;
	double PI4third;
	double PI5third;
	double PI7third;
	double PIminus_5third;
};

/*!\struct 
 * \brief For internal data transfers
 *
 * Structure to facililate parameter tranfers - All dimensionful quantities are in seconds
 */
template <class T>
struct source_parameters
{
	/*! mass of the larger component*/
	T mass1;
	/*! mass of the smaller component*/
	T mass2;
	/*! Total mass*/
	T M;
	/*! z-Spin component of the larger body*/
	T spin1z;
	/*! z-Spin component of the smaller body*/
	T spin2z;
	/*! x-Spin component of the larger body*/
	T spin1x;
	/*! x-Spin component of the smaller body*/
	T spin2x;
	/*! y-Spin component of the larger body*/
	T spin1y;
	/*! y-Spin component of the smaller body*/
	T spin2y;
	/*!Chirp mass of the binary*/
	T chirpmass;
	/*!Symmetric mass ratio*/
	T eta;
	/*!Symmetric spin combination*/
	T chi_s;
	/*!Antisymmetric spin combination*/
	T chi_a;
	/*!Effective spin */
	T chi_eff;
	/*! PN spin */
	T chi_pn;
	/*!Luminoisity Distance*/
	T DL;
	/*! Delta mass comibination*/
	T delta_mass;
	/*!Ringdown frequency after merger*/
	T fRD;
	/*!Dampening frequency after merger*/
	T fdamp;
	/*! Transition Frequency 1 for the amplitude*/	
	T f1;
	/*! Transition Frequency 2 for the amplitude*/	
	T f3;
	/*! Transition frequency 1 for the phase*/
	T f1_phase;
	/*! Transition frequency 2 for the phase*/
	T f2_phase;
	/*! Coalescence phase*/
	T phic;
	/*! Coalescence time*/
	T tc;
	/*overall amplitude factor*/	
	T A0;

	//######### ppE parameters ##############
	/*Beta factor for ppE formalism*/
	T betappe;
	/*b power for ppE formalism*/
	int bppe;

static source_parameters<T> populate_source_parameters(
			T mass1, 
			T mass2, 
			T Luminosity_Distance, 
			T *spin1,
			T *spin2, 
			T phi_c,
			T t_c) ;
};
double calculate_eta(double mass1, double mass2);
adouble calculate_eta(adouble mass1, adouble mass2);

double calculate_chirpmass(double mass1, double mass2);
adouble calculate_chirpmass(adouble mass1, adouble mass2);

double calculate_mass1(double chirpmass, double eta);
adouble calculate_mass1(adouble chirpmass, adouble eta);
	
double calculate_mass2(double chirpmass, double eta);
adouble calculate_mass2(adouble chirpmass, adouble eta);

/*!\brief Trapezoidal sum rule to approximate discrete integral - Uniform spacing
 *
 * This version is faster than the general version, as it has half the function calls
 *
 * Something may be wrong with this function - had an overall offset for real data that was
 * fixed by using the simpsons rule - not sure if this was because of a boost in accuracy or 
 * because something is off with the trapezoidal sum
 */
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += 2*integrand[i];
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/2;
}

/*!\brief Trapezoidal sum rule to approximate discrete integral - Non-Uniform spacing
 *
 * This version is slower than the uniform version, but will handle non-uniform spacing
 */
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += (integrand[i]+integrand[i-1])*(delta_x[i]-delta_x[i-1])/2;
	return sum;
}

/*!\brief Simpsons sum rule to approximate discrete integral - Uniform spacing
 *
 * More accurate than the trapezoidal rule, but must be uniform
 */
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand);
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
	{
		if( i % 2 == 0)
			sum += 2*integrand[i];
		else
			sum+=4*integrand[i];
	}
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/3;
}
long factorial(long num);
#endif
