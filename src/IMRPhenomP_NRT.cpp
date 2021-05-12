#include "IMRPhenomP_NRT.h"
#include "IMRPhenomD_NRT.h"
//#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
//#include <iostream>
//#include <cmath>
//#include <complex>
#include "util.h"

/*! \file 
 * File for the addition of tidal effects to the precessing waveform model IMRPhenomP. 
 *
 *Extends the IMRPhenomP template to include tidal effects, through inheritance 
 * from IMRPhenomD_NRT, following the NRTidal model of arXiv:1905.06011v2.
 * Specifically, equations 17, 18, 19, 20, 21, and 24 were used. 
 */

//#############################################################################
template<class T>
T IMRPhenomPv2_NRT<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow) 
{
	IMRPhenomD_NRT<T> co_prec_model;
	T out = co_prec_model.phase_ins(f, param, pn_coeff, lambda, pow);
	return out;
}

template<class T>
T IMRPhenomPv2_NRT<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda)
{
	IMRPhenomD_NRT<T> co_prec_model;
	T out =  co_prec_model.Dphase_ins(f, param, pn_coeff, lambda);
	return out;

}

template<class T>
T IMRPhenomPv2_NRT<T>::amp_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow) 
{
	IMRPhenomD_NRT<T> co_prec_model;
	T out = co_prec_model.amp_ins(f, param, pn_coeff, lambda, pow);
	return out;
}

template<class T>
T IMRPhenomPv2_NRT<T>::Damp_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda) 
{
	IMRPhenomD_NRT<T> co_prec_model;
	T out = co_prec_model.Damp_ins(f, param, pn_coeff, lambda);
	return out;
}
template class IMRPhenomPv2_NRT<double>;
template class IMRPhenomPv2_NRT<adouble>;
