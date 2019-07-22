#include "ppE_IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>

/*! \file
 *
 * Source code file for parameterized post Einsteinian Modifications to the precessing waveform model IMRPhenomP
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

//#############################################################################
template<class T>
T ppE_IMRPhenomPv2_Inspiral<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow) 
{
	ppE_IMRPhenomD_Inspiral<T> co_prec_model;
	T out = co_prec_model.phase_ins(f, param, pn_coeff, lambda, pow);
	return out;
}

template<class T>
T ppE_IMRPhenomPv2_Inspiral<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> co_prec_model;
	T out =  co_prec_model.Dphase_ins(f, param, pn_coeff, lambda);
	return out;

}

//Direct copy from PhenomPv2, except the model for the reference phase is ppE
template<class T>
void ppE_IMRPhenomPv2_Inspiral<T>::PhenomPv2_Param_Transform(source_parameters<T> *params)
{
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
	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	params->SP = params->chip * m1*m1;
	params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	T L0 = 0.0;
	useful_powers<T> pows;
	ppE_IMRPhenomD_Inspiral<T> temp;
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
	phiJ_sf = atan(J0y_sf/J0x_sf); //*NOTE* lalsuite uses "atan2" - doesn't work with adolc
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
  	params->zeta_polariz = atan(XdotQArun / XdotPArun);
	//if(std::is_same< double, T>::value){
	//	std::cout<<params->spin1z<<std::endl;
	//
	//}
	
}
//#############################################################################
template<class T>
T ppE_IMRPhenomPv2_IMR<T>::phase_mr(T f, source_parameters<T> *param,
	lambda_parameters<T> *lambda) 
{
	ppE_IMRPhenomD_IMR<T> co_prec_model;
	return co_prec_model.phase_mr(f, param, lambda);
}

template<class T>
T ppE_IMRPhenomPv2_IMR<T>::Dphase_mr(T f, source_parameters<T> *param,
	lambda_parameters<T> *lambda) 
{
	ppE_IMRPhenomD_IMR<T> co_prec_model;
	return co_prec_model.Dphase_mr(f, param, lambda);
}
template<class T>
T ppE_IMRPhenomPv2_IMR<T>::phase_int(T f, source_parameters<T> *param,
	lambda_parameters<T> *lambda) 
{
	ppE_IMRPhenomD_IMR<T> co_prec_model;
	return co_prec_model.phase_int(f, param, lambda);
}

template<class T>
T ppE_IMRPhenomPv2_IMR<T>::Dphase_int(T f, source_parameters<T> *param,
	lambda_parameters<T> *lambda) 
{
	ppE_IMRPhenomD_IMR<T> co_prec_model;
	return co_prec_model.Dphase_int(f, param, lambda);
}
//Direct copy from PhenomPv2, except the model for the reference phase is ppE
template<class T>
void ppE_IMRPhenomPv2_IMR<T>::PhenomPv2_Param_Transform(source_parameters<T> *params)
{
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
	T m1 = q/(1+q);
	T m2 = 1./(1+q);
	params->SP = params->chip * m1*m1;
	params->SL = chi1_l * m1 * m1 + chi2_l *m2 * m2;
	
	//Compute the rotation operations for L, J0 at fref	
	T L0 = 0.0;
	useful_powers<T> pows;
	ppE_IMRPhenomD_IMR<T> temp;
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
	phiJ_sf = atan(J0y_sf/J0x_sf); //*NOTE* lalsuite uses "atan2" - doesn't work with adolc
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
  	params->zeta_polariz = atan(XdotQArun / XdotPArun);
	//if(std::is_same< double, T>::value){
	//	std::cout<<params->spin1z<<std::endl;
	//
	//}
	
}
template class ppE_IMRPhenomPv2_Inspiral<double>;
template class ppE_IMRPhenomPv2_Inspiral<adouble>;
template class ppE_IMRPhenomPv2_IMR<double>;
template class ppE_IMRPhenomPv2_IMR<adouble>;
