#include "ppE_IMRPhenomD_NRT.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"

/*!\brief Overloaded method for the inspiral portion of the phase
 */
template<class T>
T ppE_IMRPhenomD_NRT<T>::phase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
	IMRPhenomD_NRT<T> model;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins_NRT(f, param, pn_coeff, lambda,powers);
	T phaseout= gr_ins;

	//return the model itself, no ppE corrections test
	//return phaseout;
	
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  pow((PIMFcube),param->bppe[i]) * param->betappe[i];

	return phaseout;

}
/*
template<class T>
T ppE_IMRPhenomD_NRT_Inspiral<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda)
{
	IMRPhenomD_NRT<T> model;
	T gr_ins = model.Dphase_ins(f, param, pn_coeff,lambda);
	T phaseout= gr_ins;

	//return the model itself, no ppE corrections test
	//return phaseout;

	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  ((T)param->bppe[i]) / 3. * pow( f ,(((T)param->bppe[i])/3.-1.))* pow((param->chirpmass *M_PI ),((T)param->bppe[i])/3.) * param->betappe[i];
	return phaseout;

	}*/
//################################################################################
//################################################################################
// Merger/Ringdown not supported for IMRPhenomD_NRT. Use Planck Taper to take amplitude to zero.


template class ppE_IMRPhenomD_NRT<double>;
template class ppE_IMRPhenomD_NRT<adouble>;
//template class ppE_IMRPhenomD_NRT_Inspiral<double>;
//template class ppE_IMRPhenomD_NRT_Inspiral<adouble>;
//template class ppE_IMRPhenomD_NRT_IMR<double>;
//template class ppE_IMRPhenomD_NRT_IMR<adouble>;
