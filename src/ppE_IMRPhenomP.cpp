#include "ppE_IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>

/*! \file
 *
 * Source code file for parameterized post Einsteinian Modifications to the precessing waveform model IMRPhenomP
 */

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
template class ppE_IMRPhenomPv2_Inspiral<double>;
template class ppE_IMRPhenomPv2_Inspiral<adouble>;
template class ppE_IMRPhenomPv2_IMR<double>;
template class ppE_IMRPhenomPv2_IMR<adouble>;
