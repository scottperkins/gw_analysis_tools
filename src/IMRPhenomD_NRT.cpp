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
  T x = pow((M_PI * param->M *f/2.), 2./3.); 
  T xpowers[5];
  for(int i = 0; i<5; i++)
    {
      xpowers[i] = pow(x, 1 + i/2.);
      //std::cout<< "x^"<< 1 + i/2.<<": "<< xpowers[i]<<"  "; 
    }
  /* std::cout<<std::endl;
   for(int j = 0; j<5; j++)
    std::cout<< "n_NRT["<< j<<"]: "<< n_NRT[j]<<std::endl;
  for(int i=0; i<3; i++) std::cout<< "d_NRT["<<i<<"]: "<< d_NRT[i]<<std::endl; 
  */
  
  T P_NRT = 1, P_NRTdenom = 1;
  for(int i = 0; i<5; i++)
    {
      P_NRT += (n_NRT[i] * xpowers[i]);
      if(i<3) P_NRTdenom += (d_NRT[i] * xpowers[i]); 
    }

  if(!deriv)
    {
      P_NRT = xpowers[3] * P_NRT/P_NRTdenom;
      std::cout<<"Pade function: " << P_NRT<<std::endl; 
      return P_NRT;
    }
  else
    {
      T DP_NRT;
      T term1, term2;

      term1 = n_NRT[0] + (3./2.)*n_NRT[1]*pow(x,1./2.) + 2*n_NRT[2]*xpowers[0] + (5./2.)*n_NRT[3]*xpowers[1] + 3*n_NRT[4]*xpowers[2]; 
      term2 = d_NRT[0] + (3./2.)*d_NRT[1]*pow(x,1./2.) + 2*d_NRT[2]*xpowers[0];
  
      DP_NRT = xpowers[3]*term1/P_NRTdenom + xpowers[3]* term2*P_NRT/(P_NRTdenom*P_NRTdenom) + (5./2.)*xpowers[1]* (P_NRT/P_NRTdenom); 

      return DP_NRT;
    }
}


template<class T> 
T IMRPhenomD_NRT<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
  IMRPhenomD<T> model;
  //T PIMFcube = pow(M_PI * param->chirpmass *f, 1./3.);
  T gr_ins = model.phase_ins(f, param, pn_coeff, lambda, powers);
  T phaseout = gr_ins;
  bool deriv = false; //tells the Pade function not to take a derivative

  phaseout = phaseout - ((3./16.) * param->tidal_weighted * (39./(16. * param->eta)) * Pade(f, param, deriv)); 
  
  return phaseout;
}


template<class T>
T IMRPhenomD_NRT<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Dphase_ins(f, param, pn_coeff, lambda);
  T phaseout = gr_ins;
  bool deriv = true;

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

  T x = pow((M_PI * param->M *f/2.), 2./3.);
  T amp_NRT = -sqrt(5*M_PI*param->eta / 24) * (9 * pow(param->M, 2.0) / param->DL) * (3./16.) * param->tidal_weighted * pow(x, 13./4.) * (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) ) / (1 + 13477.8* pow(x, 4));

  ampout += amp_NRT;
  return ampout; 
}

template<class T>
T IMRPhenomD_NRT<T>::Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
  IMRPhenomD<T> model;
  T gr_ins = model.Damp_ins(f, param, pn_coeff, lambda);
  T ampout = gr_ins; 

  T x = pow((M_PI * param->M *f/2.), 2./3.);
  T aNRT = (1 + (449./108)*x + (22672./9.) * pow(x, 2.89) );
  T aNRTdenom = (1 + 13477.8* pow(x, 4));
  T coefficient =  -sqrt(5*M_PI*param->eta / 24) * (9 * pow(param->M, 2.0) / param->DL) * (3./16.) * param->tidal_weighted;

  T DaNRT = coefficient * ( (-53911.2 * pow(x, 25./4.) * aNRT /(aNRTdenom * aNRTdenom)) + (pow(x, 13./4.) * ((449/108) + 7280.23* pow(x, 1.89))/aNRTdenom) + (13*pow(x, 9./4.) * aNRT/(4*aNRTdenom)) );

  ampout = ampout + DaNRT; 
  return ampout; 
}
template class IMRPhenomD_NRT<double>;
template class IMRPhenomD_NRT<adouble>;
