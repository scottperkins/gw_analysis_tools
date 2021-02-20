#include "IMRPhenomD_NRT.h"
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
  //T x = pow((M_PI * param->M *f*1000), 2./3.);
  //x = 0.0;
  //std::cout<<x<<","<<std::endl; 
  
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


template<class T> 
T IMRPhenomD_NRT<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
  IMRPhenomD<T> model;
  T gr_ins = model.phase_ins(f, param, pn_coeff, lambda, powers);
  T phaseout = gr_ins;
  bool deriv = false; //tells the Pade function not to take a derivative

  //std::cout<<f<<","<<std::endl; 
  //std::cout<<Pade(f, param, deriv)<<","<<std::endl; 

  //std::cout<<"{"<<f<<","<<Pade(f, param, deriv)<<"},"<<std::endl; 
  phaseout = phaseout - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv)); 
  //phaseout = phaseout + (1e6)*((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
  std::cout<<"Called IMRPhenomD_NRT::phase_ins"<<std::endl; 
  
  return phaseout;
}


template<class T>
T IMRPhenomD_NRT<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Dphase_ins(f, param, pn_coeff, lambda);
  T phaseout = gr_ins;
  bool deriv = true;

  /*
  for(int i = 0; i<5; i++)
    {
      std::cout<<"n_NRT["<<i<<"]:"<<n_NRT[i]<<std::endl;
    }
  for(int i = 0; i<3; i++)
    {
      std::cout<<"d_NRT["<<i<<"]:"<<d_NRT[i]<<std::endl;
    }
  */
  
  phaseout = phaseout - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv));
  
  return phaseout; 
}

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
  std::cout<<"Called IMRPhenomD_NRT::amp_ins"<<std::endl; 
  return ampout; 
}

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
  fmerger = (1./(2*params->M * M_PI))* 0.3586* sqrt(params->mass2 / params->mass1) *(1 + 3.354e-2 * kappa + 4.315e-5 * kappa * kappa)/(1 + 7.542e-2 * kappa + 2.236e-4* kappa * kappa);
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

template<class T>
int IMRPhenomD_NRT<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
  IMRPhenomD<T> model;
  int status = model.construct_waveform(frequencies, length, waveform, params);
  T f; 
  for(int i = 0; i<length; i++)
    {
      waveform[i] = waveform[i] * taper(frequencies[i], length, params);
    }
  return status;
}

template class IMRPhenomD_NRT<double>;
template class IMRPhenomD_NRT<adouble>;
