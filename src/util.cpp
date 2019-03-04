#include "util.h"
//#include "general_parameter_structures.h"
#include "noise_util.h"
#include <math.h>
#include <string>
#include <complex>
#include <iostream>
#include <adolc/adouble.h>

/*! \file
 *
 * General utilities that are method independent
 */

/*! \brief Builds the structure that shuttles source parameters between functions
 *
 * Populates the structure that is passed to all generation methods - contains all relavent source parameters 
 */
template <class T>
source_parameters<T> source_parameters<T>::populate_source_parameters(
			T mass1, /**< mass of the larger body - in Solar Masses*/ 
			T mass2, /**< mass of the smaller body - in Solar Masses*/
			T Luminosity_Distance,/**< Luminosity Distance in Mpc*/ 
			T *spin1,/** spin vector of the larger body  {sx,sy,sz}*/
			T *spin2, /** spin vector of the smaller body  {sx,sy,sz}*/
			T phi_c,/** coalescence phase*/
			T t_c /** coalescence time*/
			) 
{

	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters params;
	params.mass1 = mass1*MSOL_SEC;
	params.mass2 = mass2*MSOL_SEC;
	params.spin1x = spin1[0];
	params.spin2x = spin2[0];
	params.spin1y = spin1[1];
	params.spin2y = spin2[1];
	params.spin1z = spin1[2];
	params.spin2z = spin2[2];
	params.chi_s = (1./2)*(params.spin1z+params.spin2z);
	params.chi_a = (1./2)*(params.spin1z-params.spin2z);
	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
	params.chirpmass = calculate_chirpmass(params.mass1,params.mass2);
	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
	params.eta = calculate_eta(params.mass1,params.mass2);	
	params.M = params.mass1 + params.mass2;
	params.chi_eff = (params.mass1*(params.spin1z)+ params.mass2*(params.spin2z))/(params.M);
	params.chi_pn = params.chi_eff - (38*params.eta/113)*(2*params.chi_s);
	params.DL = Luminosity_Distance*MPC_SEC;
	params.delta_mass = sqrt(1.-4*params.eta);
	params.phic = phi_c;
	params.tc = t_c;
	params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	return params;
}


/*! \brief Calculates the chirp mass from the two component masses
 *
 * The output units are whatever units the input masses are
 */
double calculate_chirpmass(double mass1, double mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}
adouble calculate_chirpmass(adouble mass1, adouble mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}

/*!\brief Calculates the symmetric mass ration from the two component masses
 */
double calculate_eta(double mass1, double mass2)
{
	return (mass1 * mass2) / pow(mass1 + mass2 ,2);
}
adouble calculate_eta(adouble mass1, adouble mass2)
{
	return (mass1 * mass2) / pow(mass1 + mass2 ,2);
}

/*! \brief Calculates the larger mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass1(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}

adouble calculate_mass1(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}
/*! \brief Calculates the smaller mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass2(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}
adouble calculate_mass2(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}


long factorial(long num)
{
	int prod = 1;
	int step = num;
	while (step>0)
	{
		prod *= step;
		step -=1;
	}
	return prod;
}

double pow_int(double base, int power)
{
	double prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++)
		prod = prod * base;
	if (pow>0)
		return prod;
	else
		return 1./prod;
}
adouble pow_int(adouble base, int power)
{
	adouble prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++)
		prod = prod * base;
	if (pow>0)
		return prod;
	else
		return 1./prod;
}

template class source_parameters<double>;
template class source_parameters<adouble>;
