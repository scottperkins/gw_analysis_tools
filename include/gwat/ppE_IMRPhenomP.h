#ifndef PPE_IMRPHENOMP_H
#define PPE_IMRPHENOMP_H
#include "util.h"
#include "IMRPhenomP.h"
/*! \file
 *
 */

template<class T>
class ppE_IMRPhenomPv2_Inspiral: public IMRPhenomPv2<T>
{
public:
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda, useful_powers<T> *pow);
virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff,
	lambda_parameters<T> *lambda);

//virtual void PhenomPv2_Param_Transform(source_parameters<T> *params);

};
template<class T>
class ppE_IMRPhenomPv2_IMR: public ppE_IMRPhenomPv2_Inspiral<T>
{
public:
virtual T phase_mr(T f, source_parameters<T> *param,lambda_parameters<T> *lambda);
virtual T Dphase_mr(T f, source_parameters<T> *param,lambda_parameters<T> *lambda);
virtual T phase_int(T f, source_parameters<T> *param,lambda_parameters<T> *lambda);
virtual T Dphase_int(T f, source_parameters<T> *param,lambda_parameters<T> *lambda);
//virtual void PhenomPv2_Param_Transform(source_parameters<T> *params);

};

#endif
