#ifndef IMRPHENOMP_H
#define IMRPHENOMP_H
#include "IMRPhenomD.h"
#include "util.h"


template<class T> 
class IMRPhenomPv2: public IMRPhenomD<T>
{
public:

virtual T alpha(T omega, T q,T chi2l, T chi2);

virtual T epsilon(T omega, T q, T chi2l, T chi2);

virtual T d(int l, int mp, int m,T s);

virtual int construct_waveform(T *frequencies, 
				int length, 
				std::complex<T> *waveform_plus,
				std::complex<T> *waveform_cross,
				source_parameters<T> *params 
				);

virtual T calc_s(T f, source_parameters<T> *params);

virtual void PhenomPv2_Param_Transform(source_parameters<T> *params);

virtual T L2PN( T eta, useful_powers<T> *pow);

};

#endif
