#ifndef GIMRPHENOMD_H
#define GIMRPHENOMD_H
#include "IMRPhenomD.h"
#include "util.h"

/*! \file 
 *
 * Header file for gIMRPhenomD
 */

template<class T>
class gIMRPhenomD: public IMRPhenomD<T>
{
public:
virtual void ppE_gIMR_mapping(gen_params_base<T> *parameters, int PN_order,T *beta, int *b);
virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
	lambda_parameters<T> *lambda, useful_powers<T> *pow);
virtual T Dphase_ins(T f, source_parameters<T> *params, 
	T *pn_coeff, lambda_parameters<T> *lambda);

virtual void assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda);

virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

};

#endif
