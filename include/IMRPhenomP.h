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

virtual int construct_waveform(T *frequencies, /**< T array of frequencies the waveform is to be evaluated at*/
				int length, /**< integer length of the array of frequencies and the waveform*/	
				std::complex<T> *waveform,/**< complex T array for the waveform to be output*/ 
				source_parameters<T> *params /*Structure of source parameters to be initialized before computation*/
				);
};

#endif
