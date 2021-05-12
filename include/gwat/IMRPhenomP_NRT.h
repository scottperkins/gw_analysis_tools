#ifndef IMRPHENOMP_NRT_H
#define IMRPHENOMP_NRT_H
#include "IMRPhenomP.h"
#include "util.h"

/*! \file
 */

/*! Class that extends the IMRPhenomP waveform to include tidal effects in the 
 *inspiral portion of the phase.
 */

template<class T>
class IMRPhenomPv2_NRT: public IMRPhenomPv2<T>
{
public:
  // virtual T Pade(T f, source_parameters<T> *param, char deriv); 
  
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *pow); 
virtual T Dphase_ins(T r, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda); 

virtual T amp_ins(T f, source_parameters<T> *param, T *pn_coeff,
                        lambda_parameters<T> *lambda,useful_powers<T> *pow);
virtual T Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

};

//##############################################################################

/*!Numerically calibrated coefficients of Pade approximant (PNRTidal_v2) 
 *from arXiv:1905.06011 equations 19, 20, and 21.
 */
//const double c_NRT[4] = {3115./1248, -M_PI, 28024205./3302208., -4283.* M_PI/1092.};
//const double d_NRT[3]={-15.111208, -(c_NRT[3] + (c_NRT[1] * d_NRT[0]) - 90.550822)/c_NRT[0], 8.0641096};
//const double n_NRT[5]={c_NRT[0] + d_NRT[0], ((c_NRT[0] * c_NRT[1]) - c_NRT[3] - (c_NRT[1] * d_NRT[0]) + 90.550822)/c_NRT[0], c_NRT[2] + (c_NRT[0] * d_NRT[0]) + d_NRT[2], 90.550822, -60.253578};


#endif
