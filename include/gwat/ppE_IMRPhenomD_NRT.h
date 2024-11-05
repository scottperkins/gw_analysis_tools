#ifndef PPE_IMRPHENOMD_NRT_H
#define PPE_IMRPHENOMD_NRT_H
#include "IMRPhenomD_NRT.h"
#include "util.h"

/*! \file 
 */

/*! Class that extends the IMRPhenomD_NRT waveform to include non-GR terms in the inspiral portion of
 * the phase. We do not have merger/ringdown included in the IMRPhenomD_NRT waveform template. Instead, a taper is applied to the amplitude. 
 */
 
template<class T> 
class ppE_IMRPhenomD_NRT: public IMRPhenomD_NRT<T>
{
public:
virtual T phase_ins_NRT(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow);

  //virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

};


#endif
