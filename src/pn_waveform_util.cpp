#include "pn_waveform_util.h"
#include <math.h>
#include <adolc/adouble.h>
#include "util.h"
/*! \file
 *
 * PN waveform utilities 
 */
/*! \brief Time t as a function of f up to 2nd PN order
 *
 * Taken from https://arxiv.org/pdf/gr-qc/0411129.pdf
 *
 * Non-precessing for now
 */
template<class T>
T t_2PN(T f, T eta, T chirpmass, T chi1, T chi2, T tc)
{
	T Ldots1 = chi1;
	T Ldots2 = chi2;
	T S1hatdotS2hat = chi1*chi2;
	T m1 = calculate_mass1(chirpmass, eta);
	T m2 = calculate_mass2(chirpmass, eta);
	T M = m1+m2;
	T beta = (1./12.) * ( chi1 * (113.* m1*m1/(M*M) + 75.*eta)*Ldots1 + chi2 * (113. * m2*m2/(M*M) + 75.*eta)*Ldots2); //Spin orbit
	T sigma = (eta/48.) * chi1 * chi2 * (-247*S1hatdotS2hat + 721. * (Ldots1)*(Ldots2));//Spin spin

	T t = tc - 0.5e1 / 0.256e3 * pow(0.3141592654e1 * f, -0.8e1 / 0.3e1) * pow(chirpmass, -0.5e1 / 0.3e1) * (0.1e1 + 0.4e1 / 0.3e1 * (0.743e3 / 0.336e3 + 0.11e2 / 0.4e1 * eta) * pow(0.3141592654e1 * chirpmass * pow(eta, -0.3e1 / 0.5e1) * f, 0.2e1 / 0.3e1) - 0.8e1 / 0.5e1 * (0.4e1 * 0.3141592654e1 - beta) * 0.3141592654e1 * chirpmass * pow(eta, -0.3e1 / 0.5e1) * f + 0.2e1 * (0.3058673e7 / 0.1016064e7 + 0.5429e4 / 0.1008e4 * eta + 0.617e3 / 0.144e3 * pow(eta, 0.2e1) - sigma) * pow(0.3141592654e1 * chirpmass * pow(eta, -0.3e1 / 0.5e1) * f, 0.4e1 / 0.3e1));

	return t;
}

/*! \brief utility to get frequency from time before merger
 *
 * See https://arxiv.org/pdf/1902.00021.pdf
 */
template<class T>
T f_0PN(T t, T chirpmass)
{
	//5^(3/8)/(8 pi) 
	T factor = 0.07275685064;
	return factor * pow(chirpmass, -5./8.) * pow(t,-3./8);
}

/*! \brief Inverse of f_0PN
 */
template<class T>
T t_0PN(T f, T chirpmass)
{
	T factor = 0.07275685064;
	return pow(f/factor * pow(chirpmass, 5./8.) , -8./3.);
}
/*! \brief Utility function for the frequency at the innermost stable circular orbit (ISCO)
 *
 */
template<class T>
T FISCO(T mass)
{
	T inverse_fisco = pow_int(std::sqrt(6),3)*M_PI*mass;
	return 1./inverse_fisco;
}

//#####################################################################################
template double t_2PN<double>(double, double, double, double, double, double);
template adouble t_2PN<adouble>(adouble, adouble, adouble, adouble, adouble, adouble);
template double f_0PN<double>(double, double);
template adouble f_0PN<adouble>(adouble, adouble);
template double t_0PN<double>(double, double);
template adouble t_0PN<adouble>(adouble, adouble);
template adouble FISCO(adouble mass);
template double FISCO(double mass);
