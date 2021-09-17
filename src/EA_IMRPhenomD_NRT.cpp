#include "EA_IMRPhenomD_NRT.h"
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
 * File for the construction of waveforms in Einstein AEther theory, 
 * specifically for neutron stars. Tidal effects are included. 
 *
 * Add relevant references.
 */

template<class T>
T EA_IMRPhenomD_NRT<T>::EA_phase_ins(T f, useful_powers<T> *powers, source_parameters<T> *param)
{
  T phaseout;
  phaseout = 0; //temporary assignment while developing
  return phaseout; 
}

template<class T>
T EA_IMRPhenomD_NRT<T>::EA_amp_ins(T f, useful_powers<T> *powers, source_parameters<T> *param)
{
  T ampout;
  ampout = 0; //temporary assignment while developing
  return ampout; 
}

template<class T>
int EA_IMRPhenomD_NRT<T>::EA_construct_waveform(T *frequencies, int length, waveform_polarizations<T> *waveform, source_parameters<T> *params)
{
  std::cout<<"Used EA_construct_waveform. Print statement in line 39 of EA_IMRPhenomD_NRT.cpp"<<std::endl; 
  //IMRPhenomD_NRT<T> model;
  params->NRT_phase_coeff = - (3./16.) * params->tidal_weighted * (39./(16. * params->eta));
  if(params->tidal1<=0){params->oct1=1;params->quad1 = 1;}
  else{
    params->quad1 = this->calculate_quad_moment(params->tidal1);
    params->oct1 = this->calculate_oct_moment(params->quad1);
  }
  if(params->tidal2<=0){params->oct2=1;params->quad2 = 1;}
  else{
    params->quad2 = this->calculate_quad_moment(params->tidal2);
    params->oct2 = this->calculate_oct_moment(params->quad2);
  }

  this->calculate_spin_coefficients_3p5(params);
  params->NRT_amp_coefficient = this->calculate_NRT_amp_coefficient(params);
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
  T pn_phase_coeffs[12];
  
  this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
  this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	
  
  this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
  this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);

  //################################################################
  //Calculate phase and coalescence time variables
  T phic, f_ref, tc, phi_shift, tc_shift;
  /*Note -- to match LALSuite, this sets phiRef equal to phi_phenomD at f_ref, 
   * not phenomD_NRT. Arbitrary constant, not really important
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
	waveform->hplus[j] = 0.0;
	waveform->hcross[j] = 0.0;
	waveform->hx[j] = 0.0;
	waveform->hy[j] = 0.0;
	waveform->hb[j] = 0.0;
	waveform->hl[j] = 0.0;
      }
      else{	
	//if (f<params->f1_phase)
	//This is always needed for NRT
	if (true)
	  {
	    this->precalc_powers_ins(f, M, &pows);
	  }
	else
	  {
	    //pows.MFsixth= pow(M*f,1./6.);	
	    //pows.MF7sixth= pow_int(pows.MFsixth,7);//*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth*pows.MFsixth;
	    //Needed for Pade function
	    //pows.MFthird = pows.MFsixth*pows.MFsixth;
	    //pows.MF2third = pows.MFthird*pows.MFthird;
	  }
	amp = (A0 * this->build_amp(f,&lambda,params,&pows,pn_amp_coeffs,deltas));
	phase = (this->build_phase(f,&lambda,params,&pows,pn_phase_coeffs));
	/*Append phase_ins_NRT and amp_ins_NRT to the entire waveform*/
	{
	  T phaseNRT = this->phase_ins_NRT(f,&pows,params);
	  phase += phaseNRT;
	  T phaseSpinNRT = this->phase_spin_NRT(f, &pows,params);
	  phase += phaseSpinNRT;
	  T ampNRT = (A0*this->amp_ins_NRT(f,&pows, params));
	  amp +=ampNRT;
	}
	/*Append EA correction to the entire waveform*/
	{
	  T EAphase = this->EA_phase_ins(f, &pows, params);
	  phase += EAphase;
	  T EAamp = (A0*this->EA_amp_ins(f,&pows, params));
	  amp += EAamp; 
	}
	phase -= (T)(tc*(f-f_ref) + phic);
	waveform->hplus[j] = amp * std::exp(-i * phase);
	waveform->hcross[j] = amp * std::exp(-i * phase);
      }
      
    }
	
  //###################################### The next part applies the taper.  
  for(int i = 0; i<length; i++)
    {
      waveform->hplus[i] = waveform->hplus[i] * std::complex<T>((T)(1.0) - this->taper(frequencies[i], length, params),0);
      waveform->hcross[i] = waveform->hcross[i] * std::complex<T>((T)(1.0) - this->taper(frequencies[i], length, params),0);
    }
  
  return 1; 
}

/* Below, the construct_phase, construct_amplitude, and construct_waveform 
 * functions are overloaded with the earlier versions from parent classes.
 * Because they will not have the correct arguments and we had to define new 
 * functions, these simply call the original and warn the user that they have 
 * not been updated for Einstein AEther theory. 
 */

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  model.construct_phase(frequencies, length, phase, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_phase function. Not an updated EA version."<<std::endl;
  return 1;
}

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  model.construct_amplitude(frequencies, length, amplitude, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_amplitude function. Not an updated EA version."<<std::endl;
  return 1;
}

template<class T>
int EA_IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
  IMRPhenomD_NRT<T> model;
  model.construct_waveform(frequencies, length, waveform, params);
  std::cout<<"WARNING: code is using IMRPhenomD_NRT version of construct_waveform function. Not EA version. Search code for call to construct_waveform function and change to EA_construct_waveform."<<std::endl; 
  return 1; 
}

template class EA_IMRPhenomD_NRT<double>;
template class EA_IMRPhenomD_NRT<adouble>;
