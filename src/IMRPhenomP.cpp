#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>
#include <math.h>
#include <algorithm>



template<class T>
T IMRPhenomPv2<T>::alpha(T omega, T q,T chi2l, T chi2){

	T alpha;
	alpha = (-5*(338688*pow_int(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow_int(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow_int(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q + 168*pow_int(chi2,2)*pow_int(1 + q,2)*(3 + 4*q)))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow_int(1 + q,4)*(3 + 4*q) - 
             84672*pow_int(chi2l,2)*q*pow_int(1 + q,2)*(3 + 4*q)*
              (75 + q*(113 + 6*pow_int(chi2,2)*q*pow_int(1 + q,2))) + 
             q*(-271003598 + 127008*pow_int(chi2,4)*pow_int(q,2)*pow_int(1 + q,4)*(3 + 4*q) - 
                q*(327403764 + q*(181442579 + 39432548*q)) - 
                1512*pow_int(chi2,2)*pow_int(1 + q,2)*(213 + q*(1802 + q*(2909 + 956*q)))))) + 
       1008*omega*pow_int(1 + q,2)*(1344*M_PI*pow_int(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(-5253 + q*(-18854 + q*(-18197 - 2972*q + 168*pow_int(chi2,2)*pow_int(1 + q,2)*(3 + 4*q)))))*
        log(omega)))/(6.5028096e7*omega*q*pow_int(1 + q,4));
	return alpha;
}


template<class T>
T IMRPhenomPv2<T>::epsilon(T omega, T q, T chi2l, T chi2)
{
	T epsilon;
	epsilon = (-5*(338688*pow_int(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow_int(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow_int(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow_int(1 + q,4)*(3 + 4*q) - 
             84672*pow_int(chi2l,2)*q*pow_int(1 + q,2)*(3 + 4*q)*(75 + 113*q) - 
             q*(271003598 + q*(327403764 + q*(181442579 + 39432548*q))))) - 
       1008*omega*pow_int(1 + q,2)*(-1344*M_PI*pow_int(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(5253 + q*(18854 + q*(18197 + 2972*q))))*log(omega)))/(6.5028096e7*omega*q*pow_int(1 + q,4));
	return epsilon;
}

template<class T>
T IMRPhenomPv2<T>::d(int l, int m_prime, int m,T s)
{
	T sqrt2 = sqrt(2);
	T s_sqrt = sqrt(1+s*s);
	/*cos(beta/2)*/
	T cb = (1./sqrt2)*sqrt(1 + 1./s_sqrt);
	/*sin(beta/2)*/
	T sb = (1./sqrt2) * sqrt( 1. - 1./s_sqrt);
	
	T overall_factor = sqrt(factorial(l+m) * factorial(l-m) * factorial(l +m_prime) *factorial(l-m_prime));
	
	/*Limits of the sum (determined by keeping all factorials positive)*/
	int kmax;	
	int kmax_options[3];
	if (l == std::abs(m))
		kmax = 0;
	else
	{
		if(m < m_prime) return 0;
		else
		{
			kmax_options[0] = m-m_prime;	
			kmax_options[1] = l+m;
			kmax_options[2] = l-m;
			kmax = *std::min_element(kmax_options, kmax_options+3);
		}
	}

	/*Compute the rotation matrix*/
	int k = 0;
	T sum = 0;
	while(k<=kmax)
	{
		sum=sum+ pow_int(-1.,k+m_prime-m)/ ( factorial(l+m-k)*factorial(l-m-k)*
					factorial(k)*factorial(m_prime-m+k) ) *
					pow_int(cb,2*l-2*k-m_prime+m)*pow_int(sb,2*k+m_prime-m);
		k++;
	}
	return sum;
	
	
}

/*! \brief Constructs the waveform for IMRPhenomPv2 - uses IMRPhenomD, then twists up
 *
 * arguments:
 * 	array of frequencies, length of that array, a complex array for the output waveform, and a source_parameters structure
 */
template <class T>
int IMRPhenomPv2<T>::construct_waveform(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length, /**< integer length of the array of frequencies and the waveform*/	
				std::complex<T> *waveform,/**< complex T array for the waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{
	T M = params-> M;
	T chirpmass = params->chirpmass;
	T DL = params->DL;
	lambda_parameters<T> lambda, *lambda_ptr;
	this->assign_lambda_param(params, &lambda);

	/*Initialize the post merger quantities*/
	this->post_merger_variables(params);

	params->f1_phase = 0.018/(params->M);
	params->f2_phase = params->fRD/2.;

	params->f1 = 0.014/(params->M);
	params->f3 = this->fpeak(params, &lambda);

	useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

	T deltas[6];
	T pn_amp_coeffs[7];
	T pn_phase_coeffs[8];

	this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	

	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);


	//T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
	T A0 = params->A0;
	T s = .1;
	T chi2l = .1;
	T chi2 = .1;
	int m = 2;
	T q = params->mass2/params->mass1;

	T f;
	std::complex<T> amp, phase;
	//T amp, phase;
	std::complex<T> i;
	i = std::complex<T> (0,1.);
	for (int j =0; j< length; j++)
	{
		f = frequencies[j];
		if (f<params->f1_phase)
		{
			this->precalc_powers_ins(f, M, &pows);
		}
		amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
		phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
		amp = amp * this->d(2,2,2,s);
		phase = phase + (std::complex<T>)(2 * this->epsilon(M_PI*f, q, chi2l,chi2) 
				+ m * this->alpha(M_PI*f,q,chi2l,chi2));
		waveform[j] = amp * std::exp(-i * phase);

	}
	//}
	return 1;
}

template class IMRPhenomPv2<double>;
template class IMRPhenomPv2<adouble>;
