#include "IMRPhenomD_NRT.h"
#include "IMRPhenomD.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"

/*! \file 
 * File for the addition of tidal effects to the waveform. 
 *
 *Extends the IMRPhenomD template to include tidal effects, following the 
 * NRTidal model of arXiv:1905.06011v2. Specifically, equations 17, 18, 19, 20, 21, 
 * and 24 were used. 
 */

//##############################################################################
//##############################################################################
template<class T>
T IMRPhenomD_NRT<T>::Pade(T f, source_parameters<T> *param, char deriv)
{ 
  T x = pow((M_PI * param->M *f), 2./3.);
  
  T xpowers[5];
  for(int i = 0; i<5; i++)
    {
      xpowers[i] = pow(x, 1 + i/2.);
    }
  
  T P_NRT = 1, P_NRTdenom = 1;
  for(int i = 0; i<5; i++)
    {
      P_NRT += (n_NRT[i] * xpowers[i]);
      if(i<3) P_NRTdenom += (d_NRT[i] * xpowers[i]); 
    }

  if(!deriv)
    {
      P_NRT = xpowers[3] * P_NRT/P_NRTdenom; //Note that this is the Pade function times x^(5/2).
      //std::cout<<P_NRT<<","<<std::endl; 
      
      return P_NRT;
    }
  else
    {
      T DP_NRT;
      T term1, term2;

      term1 = n_NRT[0] + (3./2.)*n_NRT[1]*pow(x,1./2.) + 2*n_NRT[2]*xpowers[0] + (5./2.)*n_NRT[3]*xpowers[1] + 3*n_NRT[4]*xpowers[2]; 
      term2 = d_NRT[0] + (3./2.)*d_NRT[1]*pow(x,1./2.) + 2*d_NRT[2]*xpowers[0];
  
      DP_NRT = xpowers[3]*term1/P_NRTdenom + xpowers[3]* term2*P_NRT/(P_NRTdenom*P_NRTdenom) + (5./2.)*xpowers[1]* (P_NRT/P_NRTdenom);
      //This is a derivative with respect to x. To get a derivative with respect to f, need to multiply by dx/df
      DP_NRT = DP_NRT* (2/3.)*pow((param->M*param->M* M_PI*M_PI/f), 1/3.);

      return DP_NRT;
    }
}

/*Old version*/
//template<class T> 
//T IMRPhenomD_NRT<T>::phase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *powers)
//{
//  IMRPhenomD<T> model;
//  T gr_ins = model.phase_ins(f, param, pn_coeff, lambda, powers);
//  T phaseout = gr_ins;
//  
//  T phasefix; 
//  bool deriv = false; //tells the Pade function not to take a derivative
//
//  phasefix = - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
//  phaseout = phaseout + phasefix;
//  
//  return phaseout;
//}

/* I rewrote this so it's NOT an overloaded version of phase_ins -- completely new. Only calculates the NRT phase, not the total phase_ins
 */
template<class T> 
T IMRPhenomD_NRT<T>::phase_ins_NRT(T f, source_parameters<T> *param)
{
  
  /*Note that this is just the NRT part now*/
  T phaseout= 0;
  T phasefix; 
  bool deriv = false; //tells the Pade function not to take a derivative

  phasefix = - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
  phaseout = phaseout + phasefix;
  
  return phaseout;
}


/* This shouldn't be needed anymore, but I'll leave it in. At the moment, it's not called at all (changed the name so its not overloaded)
 */
template<class T>
T IMRPhenomD_NRT<T>::Dphase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Dphase_ins(f, param, pn_coeff, lambda);
  T phaseout = gr_ins;
  bool deriv = true;
  
  phaseout = phaseout - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
  return gr_ins; 
}

/*For the amplitude, you may want to call build_amp instead of each amplitude component (ins,int,mr) separately, since LALsuite seems to append the correction
 * to the entire waveform. Up to you though. You would overload that function exactly like the others
 */

/*Still be used as overloaded function*/
template<class T>
T IMRPhenomD_NRT<T>::amp_ins(T f, source_parameters<T> *param, T *pn_coeff,
                        lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
  IMRPhenomD<T> model;
  T gr_ins = model.amp_ins(f, param, pn_coeff, lambda, powers);
  T ampout = gr_ins;

  T x = pow((M_PI * param->M *f), 2./3.);

  T amp_NRT = -sqrt(5*M_PI*param->eta / 24.) * (9 * param->M * param->M / param->DL) * (3./16.) * param->tidal_weighted * pow(x, 13./4.) * (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) ) / (1 + 13477.8* pow(x, 4));
  //has to be scaled by f^(7/6)/A0 to be consistent with the rest of GW analysis tools
  amp_NRT = amp_NRT*powers->MF7sixth/ (param->A0*pow(param->M, 7/6.));

  ampout += amp_NRT;
  //std::cout<<"Called IMRPhenomD_NRT::amp_ins"<<std::endl; 
  return ampout; 
}

/*Still be used as overloaded function*/
template<class T>
T IMRPhenomD_NRT<T>::Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Damp_ins(f, param, pn_coeff, lambda);
  T ampout = gr_ins; 

  T x = pow((M_PI * param->M *f), 2./3.);
  
  T aNRT = (1 + (449./108.)*x + (22672./9.) * pow(x, 2.89) );
  T aNRTdenom = (1 + 13477.8* pow(x, 4));
  T coefficient =  -sqrt(5*M_PI*param->eta / 24.) * (9 * pow(param->M, 2.0) / param->DL) * (3./16.) * param->tidal_weighted;

  T amp_NRT = coefficient* pow(x,13./4.)* aNRT/aNRTdenom;

  T DaNRT = coefficient * ( (-53911.2 * pow(x, 25./4.) * aNRT /(aNRTdenom * aNRTdenom)) + (pow(x, 13./4.) * ((449./108.) + 7280.23* pow(x, 1.89))/aNRTdenom) + (13*pow(x, 9./4.) * aNRT/(4.*aNRTdenom)) );
  //This is actually the derivative of the amplitude with respect to x
  //to make it the derivative of the amplitude with respect to f, need to multiply by dx/df, done in the next line  
  DaNRT = DaNRT* (2/3.)*pow((param->M*param->M* M_PI*M_PI/f), 1/3.);

  DaNRT = DaNRT *pow(param->M * f, 7/6.)/ (param->A0*pow(param->M, 7/6.)) + amp_NRT * (7/6.)*pow(f,1/6.)/ (param->A0);
  //with appropriate scaling factors
  std::cout<<"DaNRT ="<<DaNRT<<std::endl;
  debugger_print(__FILE__,__LINE__,DaNRT);
  ampout = ampout + DaNRT; 
  return ampout; 
}

template<class T>
T IMRPhenomD_NRT<T>::taper(T f, int length, source_parameters<T> *params)
{
  //std::cout<<"Called IMRPhenomD_NRT::taper"<<std::endl; 
  IMRPhenomD<T> model;
  T sigma;
  T fmerger, kappa, kappa_eff;
  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
  
  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.
  
  kappa = kappa_eff;
  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa); //Equation 11 of arXiv:1804.02235
  T fmerger12 = fmerger*1.2;
  T z = (fmerger - fmerger12)/(f - fmerger) + (fmerger - fmerger12)/(f - fmerger12); 
  
  if(f < fmerger)
    {
      return 1.0; 
    }
  else if(fmerger < f && f < fmerger12)
    {
      sigma = 1.0/(exp(z) + 1.0);
      return sigma; 
    }
  else if(f > fmerger12)
    {
      return 0.0; 
    }

  return -1.0; //If it gets to this point, something is wrong and the change in sign should indicate that
}

/*Just call the phenomD construct_waveform (still calls your version of amp though, because that's still overloaded) and tacks on the phase at the end with the taper
 * Up to you on how to do the amplitude portion
 */
template<class T>
int IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  //int status = model.construct_waveform(frequencies, length, waveform, params);
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

  /* T fmerger, kappa, kappa_eff;
  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
  
  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.
  
  kappa = kappa_eff;
  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa);

  params->f3 = fmerger;
  */
  //Wanted to try using the fmerger defined in the NR Tidal paper instead of fpeak. This made a huge difference...
	
  useful_powers<T> pows;
  this->precalc_powers_PI(&pows);

  T deltas[6];
  T pn_amp_coeffs[7];
  T pn_phase_coeffs[12];
  
  this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
  this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	
  
  this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
  this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
  
  //################################################################
  //Calculate phase and coalescence time variables
  T phic, f_ref, tc, phi_shift, tc_shift;
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, not phenomD_NRT
 * Arbitrary constant, not really important
 */
  if(params->shift_phase ){
    f_ref = params->f_ref;
    this->precalc_powers_ins(f_ref, M, &pows);
    phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
    phic = 2*params->phiRef + phi_shift;
  }
  //If phic is specified, ignore f_ref phiRef and use phic
  else{
    f_ref = 0;
    phic = params->phiRef;
  }
  
  //Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
  //This aligns more with the physical meaning of tc, but the phase is NO LONGER just
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, not phenomD_NRT
 * Arbitrary constant, not really important
 */
  if(params->shift_time){
    T alpha1_offset = this->assign_lambda_param_element(params,14);
    tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
  }
  else{
    tc_shift=0;
  }
  tc = 2*M_PI*params->tc + tc_shift;
  
  T A0 = params->A0* pow(M,7./6.);
  
  T f;
  std::complex<T> amp, phase;
  std::complex<T> i;
  i = std::complex<T> (0,1.);
  T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
  for (size_t j =0; j< length; j++)
    {
      f = frequencies[j];
      if(f>fcut){
	amp = 0.0;
	waveform[j] = 0.0;
      }
      else{	
	if (f<params->f1_phase)
	  {
	    this->precalc_powers_ins(f, M, &pows);
	  }
	else
	  {
	    pows.MFsixth= pow(M*f,1./6.);	
	    pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
	  }
	amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
	phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
	/*Append phase_ins_NRT to the entire waveform*/
	{
		T phaseNRT = this->phase_ins_NRT(f, params);
		phase += phaseNRT;
	}
	//phase +=   (T)(tc*(f-f_ref) - phic);
	phase -=   (T)(tc*(f-f_ref) + phic);
	waveform[j] = amp * std::exp(-i * phase);
      }
      
    }
	
  //###################################### The next part applies the taper.  
  for(int i = 0; i<length; i++)
    {
      waveform[i] = waveform[i] * taper(frequencies[i], length, params);
    }
  
  return 1;
}

/*Old version*/

//template<class T>
//int IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
//{
////TESTING 
//  T kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
//  T kappa = kappa_eff;
//  T fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa); //Equation 11 of arXiv:1804.02235
////#############################
//  IMRPhenomD<T> model;
//  //int status = model.construct_waveform(frequencies, length, waveform, params);
//  T M = params-> M;
//  T chirpmass = params->chirpmass;
//  T DL = params->DL;
//  lambda_parameters<T> lambda, *lambda_ptr;
//  this->assign_lambda_param(params, &lambda);
//  
//  /*Initialize the post merger quantities*/
//  this->post_merger_variables(params);
//  params->f1_phase = 0.018/(params->M);
//  params->f2_phase = params->fRD/2.;
//  
//  params->f1 = 0.014/(params->M);
//  params->f3 = this->fpeak(params, &lambda);
//
//  /* T fmerger, kappa, kappa_eff;
//  kappa_eff = (3./16.) * params->tidal_weighted; //kappa effective as defined in arXiv:1804.02235 and arXiv:1905.06011v2.
//  
//  //kappa = 3.*(params->mass2 * pow(params->mass1, 4.)* params->tidal1/ pow(params->M,5.) + params->mass1 * pow(params->mass2, 4.)* params->tidal2/ pow(params->M,5.));
//  //This is the definition for kappa_2^T given in arXiv:1804.02235. As described in that paper,
//  //it is accurate enough to use kappa effective. Can test by turning this back on. No change in output.
//  
//  kappa = kappa_eff;
//  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa);
//
//  params->f3 = fmerger;
//  */
//  //Wanted to try using the fmerger defined in the NR Tidal paper instead of fpeak. This made a huge difference...
//	
//  useful_powers<T> pows;
//  this->precalc_powers_PI(&pows);
//
//  T deltas[6];
//  T pn_amp_coeffs[7];
//  T pn_phase_coeffs[12];
//  
//  this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
//  this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	
//  
//  this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
//  this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);
//  
//  //################################################################
//  //Calculate phase and coalescence time variables
//  T phic, f_ref, tc, phi_shift, tc_shift;
//  //If phic is unspecified - use f_ref and phiRef
//  if(params->shift_phase ){
//    f_ref = params->f_ref;
//    this->precalc_powers_ins(f_ref, M, &pows);
//    //NOTE : this is weird, but seems to agree with LIGO
//    //phi_shift = (this->build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
//    phi_shift = (model.build_phase(f_ref,&lambda,params,&pows,pn_phase_coeffs));
//    phic = 2*params->phiRef + phi_shift;
//  }
//  //If phic is specified, ignore f_ref phiRef and use phic
//  else{
//    f_ref = 0;
//    phic = params->phiRef;
//  }
//  
//  //Assign shift: first shift so coalescence happens at t=0, then shift from there according to tc
//  //This aligns more with the physical meaning of tc, but the phase is NO LONGER just
//  if(params->shift_time){
//    T alpha1_offset = this->assign_lambda_param_element(params,14);
//    //NOTE : this is weird, but seems to agree with LIGO
//    //tc_shift = this->Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
//    tc_shift = model.Dphase_mr(params->f3, params, &lambda)+(-lambda.alpha[1]+alpha1_offset)*params->M/params->eta;
//  }
//  else{
//    tc_shift=0;
//  }
//  //tc = 2*M_PI*params->tc - tc_shift;
//  tc = 2*M_PI*params->tc + tc_shift;
//  
//  //T A0 = sqrt(M_PI/30)*chirpmass*chirpmass/DL * pow(M_PI*chirpmass,-7./6);
//  T A0 = params->A0* pow(M,7./6.);
//  
//  T f;
//  std::complex<T> amp, phase;
//  std::complex<T> i;
//  i = std::complex<T> (0,1.);
//  T fcut = .2/M; //Cutoff frequency for IMRPhenomD - all higher frequencies return 0
//  for (size_t j =0; j< length; j++)
//    {
//      f = frequencies[j];
//      if(f>fcut){
//	amp = 0.0;
//	waveform[j] = 0.0;
//      }
//      else{	
//	if (f<params->f1_phase)
//	  {
//	    this->precalc_powers_ins(f, M, &pows);
//	  }
//	else
//	  {
//	    pows.MFsixth= pow(M*f,1./6.);	
//	    pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
//	  }
//	amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
//	phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
//	if(true){
//		T phaseNRT = this->phase_ins_NRT(f, params,pn_phase_coeffs,&lambda,&pows);
//		phase += phaseNRT;
//	}
//	//phase +=   (T)(tc*(f-f_ref) - phic);
//	phase -=   (T)(tc*(f-f_ref) + phic);
//	waveform[j] = amp * std::exp(-i * phase);
//      }
//      
//    }
//	
//  //###################################### The next part applies the taper.  
//  for(int i = 0; i<length; i++)
//    {
//      waveform[i] = waveform[i] * taper(frequencies[i], length, params);
//    }
//  
//  return 1;
//}

template class IMRPhenomD_NRT<double>;
template class IMRPhenomD_NRT<adouble>;
