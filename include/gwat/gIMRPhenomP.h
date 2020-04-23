#ifndef GIMRPHENOMP_H
#define GIMRPHENOMP_H
#include "util.h"
#include "IMRPhenomP.h"
/*! \file
 *
 */

template<class T>
class gIMRPhenomPv2: public IMRPhenomPv2<T>
{
public:
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow);
virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda);
virtual void assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda);

virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

//virtual void PhenomPv2_Param_Transform(source_parameters<T> *params);

};
#endif
