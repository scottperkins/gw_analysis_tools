#ifndef GIMRPHENOMD_H
#define GIMRPHENOMD_H
#include "IMRPhenomD.h"
#include "util.h"

template<class T>
class gIMRPhenomD: public IMRPhenomD<T>
{
public:
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
	lambda_parameters<T> *lambda, useful_powers<T> *pow);
virtual T Dphase_ins(T f, source_parameters<T> *params, 
	T *pn_coeff, lambda_parameters<T> *lambda);

virtual void assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda);

virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

};

#endif
