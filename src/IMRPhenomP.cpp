#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>
#include <math.h>
#include <algorithm>
/*! \file 
 *
 * Source code for IMRPhenomP
 */

//Shamelessly stolen from lalsuite
/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
tmp1 = vx*cos(angle) - vy*sin(angle);\
tmp2 = vx*sin(angle) + vy*cos(angle);\
vx = tmp1;\
vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
tmp1 = vx*cos(angle) + vz*sin(angle);\
tmp2 = - vx*sin(angle) + vz*cos(angle);\
vx = tmp1;\
vz = tmp2



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
		//Probably mulitply frequency by M here..
		phase = phase + (std::complex<T>)(2 * this->epsilon(M_PI*f, q, chi2l,chi2) 
				+ m * (this->alpha(M_PI*f,q,chi2l,chi2)  +params->alpha0));
		waveform[j] = amp * std::exp(-i * phase);

	}
	//}
	return 1;
}


/*! /Brief Parameter transformtion to precalculate needed parameters for PhenomP from source parameters
 *
 * Pretty much stolen verbatim from lalsuite
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_Param_Transform(source_parameters<T> *params /*< Source Parameters*/
						)
{
	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1z;
	T chi2_l = params->spin2z;
	
	T m1_2 = params->mass1 * params->mass1;
	T m2_2 = params->mass2 * params->mass2;
	
	T S1_perp = m1_2 * sqrt( params->spin1y* params->spin1y +
				params->spin1x * params->spin1x);
	T S2_perp = m2_2 * sqrt( params->spin2y* params->spin2y +
				params->spin2x * params->spin2x);

	T A1 = 2 + (3*params->mass2)/ ( 2 * params->mass1);
	T A2 = 2 + (3*params->mass1)/ ( 2 * params->mass2);
	T ASp1 = A1*S1_perp;
	T ASp2 = A2*S2_perp;
	T num = (ASp2>ASp1) ? ASp2 : ASp1;
	T denom = (params->mass2 > params->mass1)? A2*m2_2 : A1*m1_2;
	params->chip = num/denom;

	//Compute the rotation operations for L, J0 at fref	
	//T v_ref = pow(M_PI * params->f_ref * (params->M),1./3);
	T L0 = 0.0;
	useful_powers<T> pows;
	IMRPhenomD<T> temp;
	temp.precalc_powers_PI(&pows);
	temp.precalc_powers_ins(params->f_ref, params->M, &pows);
	
	L0 = params->M * params-> M * this->L2PN(params->eta, &pows);
	
	//_sf denotes source frame - ie L_hat propto z_hat
	T J0x_sf = m1_2 * params->spin1x + m2_2 * params->spin2x;	
	T J0y_sf = m1_2 * params->spin1y + m2_2 * params->spin2y;	
	T J0z_sf = L0 + m1_2* params->spin1z + m2_2 * params->spin2z;
		
	T J0 = sqrt(J0x_sf * J0x_sf + J0y_sf * J0y_sf + J0z_sf * J0z_sf ) ;
	
	//thetaJ_sf is the angle between J0 and L (zhat)
	T thetaJ_sf;
	thetaJ_sf = acos(J0z_sf/J0);

	//azimuthal angle of J0 in the source frame
	T phiJ_sf;
	phiJ_sf = atan(J0y_sf/J0x_sf); //*NOTE* lalsuite uses "atan2" - not standard
	params->phi_aligned = - phiJ_sf;

	//Rotation of the system s.t. the total J is pointed in zhat
	T tmp1,tmp2;
	T incl = params->incl_angle;
	T phiRef = params->phiRef;
	T Nx_sf = sin(incl) * cos(M_PI/2. - phiRef);
	T Ny_sf = sin(incl) * sin(M_PI/2. - phiRef);
	T Nz_sf = cos(incl);
	T tmp_x = Nx_sf;
	T tmp_y = Ny_sf;
	T tmp_z = Nz_sf;
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	T kappa;
	kappa = -atan(tmp_y/tmp_x);

	//alpha0
	tmp_x = 0.;
	tmp_y = 0.;
	tmp_z = 1.;
	
	ROTATEZ(-phiJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEY(-thetaJ_sf, tmp_x,tmp_y, tmp_z);
	ROTATEZ(kappa, tmp_x,tmp_y, tmp_z);
	params->alpha0 = atan(tmp_y/tmp_x);

	tmp_x = Nx_sf;
  	tmp_y = Ny_sf;
  	tmp_z = Nz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	T Nx_Jf = tmp_x; // let's store those two since we will reuse them later (we don't need the y component)
  	T Nz_Jf = tmp_z;
  	params->thetaJN = acos(Nz_Jf);

	/* Finally, we need to redefine the polarizations :
	   PhenomP's polarizations are defined following Arun et al (arXiv:0810.5336)
	   i.e. projecting the metric onto the P,Q,N triad defined with P=NxJ/|NxJ| (see (2.6) in there).
	   By contrast, the triad X,Y,N used in LAL
	   ("waveframe" in the nomenclature of T1500606-v6)
	   is defined in e.g. eq (35) of this document
	   (via its components in the source frame; note we use the defautl Omega=Pi/2).
	   Both triads differ from each other by a rotation around N by an angle \zeta
	   and we need to rotate the polarizations accordingly by 2\zeta
	  */

	T Xx_sf = -cos(incl)*sin(phiRef);
  	T Xy_sf = -cos(incl)*cos(phiRef);
  	T Xz_sf = sin(incl);
  	tmp_x = Xx_sf;
  	tmp_y = Xy_sf;
  	tmp_z = Xz_sf;
  	ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
  	ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
  	//now the tmp_a are the components of X in the J frame
  	//we need the polar angle of that vector in the P,Q basis of Arun et al
  	// P=NxJ/|NxJ| and since we put N in the (pos x)z half plane of the J frame
  	T PArunx_Jf = 0.;
  	T PAruny_Jf = -1.;
  	T PArunz_Jf = 0.;
  	// Q=NxP
  	T QArunx_Jf = Nz_Jf;
  	T QAruny_Jf = 0.;
  	T QArunz_Jf = -Nx_Jf;
  	T XdotPArun = tmp_x*PArunx_Jf+tmp_y*PAruny_Jf+tmp_z*PArunz_Jf;
  	T XdotQArun = tmp_x*QArunx_Jf+tmp_y*QAruny_Jf+tmp_z*QArunz_Jf;
  	params->zeta_polariz = atan2(XdotQArun , XdotPArun);
	
}

template<class T>
T IMRPhenomPv2<T>::L2PN( T eta, useful_powers<T> *pow)
{
	T x = pow->MF2third * pow->PI2third;
	return eta/(sqrt(x))*(1 + (3./2 + eta/6)*x + (3.373 - 19.*eta/8 - eta*eta/24 )*x*x);
}

template class IMRPhenomPv2<double>;
template class IMRPhenomPv2<adouble>;
