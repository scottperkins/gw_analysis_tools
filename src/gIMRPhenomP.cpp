#include "gIMRPhenomP.h"
#include "gIMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>

/*! \file
 *
 * Source code file for parameterized post Einsteinian Modifications to the precessing waveform model IMRPhenomP
 */

//#############################################################################
template<class T>
T gIMRPhenomPv2<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow)
{
	gIMRPhenomD<T> model;
	return model.phase_ins(f,param, pn_coeff, lambda, pow);
}
template<class T>
T gIMRPhenomPv2<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda)
{
	gIMRPhenomD<T> model;
	return model.Dphase_ins(f,param, pn_coeff, lambda);
}
template<class T>
void gIMRPhenomPv2<T>::assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda)
{
	gIMRPhenomD<T> model;
	model.assign_lambda_param(source_param, lambda );
	return;
}


template<class T>
void gIMRPhenomPv2<T>::assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff)
{
	gIMRPhenomD<T> model;
	model.assign_static_pn_phase_coeff(source_param, coeff );
	return;

}
//#############################################################################

template class gIMRPhenomPv2<double>;
template class gIMRPhenomPv2<adouble>;
