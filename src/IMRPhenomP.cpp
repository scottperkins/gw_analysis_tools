#include "IMRPhenomP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include "IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>
#include <math.h>
#include <algorithm>
#include <type_traits>
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

const double sqrt_6 = 2.44948974278317788;

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
				std::complex<T> *waveform_plus,/**< complex T array for the plus polariaztion waveform to be output*/ 
				std::complex<T> *waveform_cross,/**< complex T array for the cross polarization waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				)
{
	//Initialize Spherical harmonics for polarization construction
	sph_harm<T> harmonics;
	T phiHarm = 0.;
	harmonics.Y22 = std::complex<T>(0.0,0.0);
	harmonics.Y21 = std::complex<T>(0.0,0.0);
	harmonics.Y20 = std::complex<T>(0.0,0.0);
	harmonics.Y2m1 =std::complex<T>(0.0,0.0);
	harmonics.Y2m2 =std::complex<T>(0.0,0.0);
	harmonics.Y22 = XLALSpinWeightedSphericalHarmonic(params->incl_angle,phiHarm, -2,2,2);
	harmonics.Y21 = XLALSpinWeightedSphericalHarmonic(params->incl_angle,phiHarm, -2,2,1);
	harmonics.Y20 = XLALSpinWeightedSphericalHarmonic(params->incl_angle,phiHarm, -2,2,0);
	harmonics.Y2m1 = XLALSpinWeightedSphericalHarmonic(params->incl_angle,phiHarm, -2,2,-1);
	harmonics.Y2m2 = XLALSpinWeightedSphericalHarmonic(params->incl_angle,phiHarm, -2,2,-2);
	
	//if(std::is_same<T,double>::value){
	//	std::cout<<"sph harm: "<<harmonics.Y22<<" "<<harmonics.Y21<<" "<<harmonics.Y20<<" "<<harmonics.Y2m2<<" "<<harmonics.Y2m1<<" "<<std::endl;
	//}
	
	//if(std::is_same<T,double>::value){
	//	std::cout<<params->alpha0<<std::endl;
	//}
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


	//Rescale amplitude because lalsuite does
	T A0 = params->A0 / (2. * sqrt(5. / (64.*M_PI)) );
	T q = params->mass1/params->mass2;
	T d2[5] ;
	T dm2[5];
	std::complex<T> hp_factor = std::complex<T>(0.0,0.0);
	std::complex<T> hc_factor = std::complex<T>(0.0,0.0);
	T epsilon;
	T alpha;

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
		else{
			pows.MF2third = pow(M * f, 2./3.);//PhenomP requires this for WignerD matrices
		}
		amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
		phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
		//Calculate WignerD matrices -- See mathematica nb for the forms: stolen from lalsuite
		this->WignerD(d2,dm2, &pows, params);
		//Calculate Euler angles alpha and epsilon
		//epsilon = this->epsilon(M_PI*f, q, params->chil,params->chip);
		this->calculate_euler_angles(&alpha, &epsilon, params->M *M_PI*f, q, params->chil, params->chip);
		//Add offset to alpha
		alpha = alpha - params->alpha0;
		//if(std::is_same<T,double>::value){
		//	std::cout<<d2[0]<<" "<<d2[1]<<" "<<d2[2]<<" "<<d2[3]<<" "<<d2[4]<<std::endl;
		//}

		//Twist it up
		calculate_twistup(alpha, &hp_factor, &hc_factor, d2, dm2, &harmonics);
		//if(std::is_same<T,double>::value){
		//	std::cout<<hp_factor<<hc_factor<<std::endl;
		//}

		//if(std::is_same<T,double>::value){
		//	std::cout<<"hp: "<<hp_factor<<" hc: "<<hc_factor<<std::endl;
		//}
		//Probably mulitply frequency by M here..
		//hp_factor=1;
		//hc_factor=1;
		phase = phase + (std::complex<T>)(2. * epsilon);
		
		waveform_plus[j] = amp *hp_factor *  std::exp(-i * phase)/std::complex<T>(2.,0.0);
		waveform_cross[j] = amp *hc_factor *  std::exp(-i * phase)/std::complex<T>(2.,0.0);
		hp_factor = 0.;
		hc_factor = 0.;

	}
	//}
	return 1;
}

template<class T>
void IMRPhenomPv2<T>::WignerD(T d2[5],T dm2[5], useful_powers<T> *pows, source_parameters<T> *params)
{
	T L = this->L2PN(params->eta, pows) * params->M * params->M;
	T s = params->SP / ( L + params-> SL);
	T s_2 = s*s;
	//if(std::is_same<T,double>::value){
	//	std::cout<<"L: "<<L<<" s: "<<s<<" SP: "<<params->SP<<std::endl;
	//}
	T cos_beta = 1./sqrt(1.0 + s_2);
	T cos_beta_half = sqrt( (1.0 + cos_beta) / 2.0);
	T sin_beta_half = sqrt( (1.0 - cos_beta) / 2.0);
	T c2 = cos_beta_half * cos_beta_half;
	T s2 = sin_beta_half * sin_beta_half;
	T c3 = c2 * cos_beta_half;
	T s3 = s2 * sin_beta_half;
	T c4 = c3 * cos_beta_half;
	T s4 = s3 * sin_beta_half;
	
	d2[0] = s4;
	d2[1] = 2*cos_beta_half * s3;
	d2[2] = sqrt_6 * s2*c2 ;
	d2[3] = 2 * c3* sin_beta_half;
	d2[4] = c4;
	//Exploit Symmetry
	dm2[0] = d2[4];
	dm2[1] = -d2[3];
	dm2[2] = d2[2];
	dm2[3] = -d2[1];
	dm2[4] = d2[0];
	
}

template<class T>
void IMRPhenomPv2<T>::calculate_twistup( T alpha, std::complex<T> *hp_factor, std::complex<T> *hc_factor, T d2[5], T dm2[5], sph_harm<T> *sph_harm)
{
	std::complex<T> T2m;
	std::complex<T> Tm2m;
	std::complex<T> exp_a = std::exp(std::complex<T>(0,1) *alpha);
	std::complex<T> exp_ma = std::complex<T>(1.,0.0)/exp_a;
	std::complex<T> exp_2a = exp_a*exp_a;
	std::complex<T> exp_m2a = exp_ma*exp_ma;
	std::complex<T> exp_a_vec[5] = {exp_m2a, exp_ma,std::complex<T>(1.0,0.0) , exp_a, exp_2a};
	std::complex<T> harmonics[5] = 
			{sph_harm->Y2m2, sph_harm->Y2m1, sph_harm->Y20, sph_harm->Y21, sph_harm->Y22};
	for (int m = -2; m<=2; m++)
	{
		T2m = exp_a_vec[-m+2] * dm2[m+2] * harmonics[m+2];
		Tm2m = exp_a_vec[m+2] * d2[m+2] * conj(harmonics[m+2]); 
		*hp_factor += T2m + Tm2m;
		*hc_factor += std::complex<T>(0,1.) * (T2m - Tm2m);
	}
	
}
template<class T>
void IMRPhenomPv2<T>::calculate_euler_angles(T *alpha, T *epsilon, T omega, T q, T chil, T chip)
{
	*alpha = this->alpha(omega,q, chil, chip) ;	
	*epsilon = this->epsilon(omega,q,chil,chip);
}
/*! /Brief Parameter transformtion to precalculate needed parameters for PhenomP from source parameters
 *
 * Pretty much stolen verbatim from lalsuite
 */
template<class T>
void IMRPhenomPv2<T>::PhenomPv2_Param_Transform(source_parameters<T> *params /*< Source Parameters*/
						)
{
	//TESTING
	params->phiRef=0.0;
	params->f_ref=10.;


	//Calculate spin parameters chil and chip
	T chi1_l = params->spin1z;
	T chi2_l = params->spin2z;

	T q = params->mass1/params->mass2;
	T chi_eff = (params->mass1*chi1_l + params->mass2*chi2_l) / params->M; /* Effective aligned spin */
  	T chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
	params->chil = chil;


	
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
	params->SP = params->chip * params->mass1 * params->mass1;
	params->SL = chi1_l * params->mass1 * params->mass1 + chi2_l *params->mass2 * params->mass2;
	//if(std::is_same<T,double>::value){
	//	std::cout<<params->chil<<std::endl;
	//}
	
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
	
	//if(std::is_same<T,double>::value){
	//	std::cout<<tmp_y<<tmp_x<<std::endl;
	//}
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
  	params->zeta_polariz = atan(XdotQArun / XdotPArun);
	
}

template<class T>
T IMRPhenomPv2<T>::L2PN( T eta, useful_powers<T> *pow)
{
	T x = pow->MF2third * pow->PI2third;
	T x2 = x*x;
	T eta2 = eta*eta;
	//return eta/(sqrt(x))*(1 + (3./2 + eta/6)*x + (3.373 - 19.*eta/8 - eta*eta/24 )*x*x);
	return (eta*(1.0 + (1.5 + eta/6.0)*x + (3.375 - (19.0*eta)/8. - eta2/24.0)*x2)) / sqrt(x);
}

template class IMRPhenomPv2<double>;
template class IMRPhenomPv2<adouble>;
