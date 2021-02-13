#ifndef TAYLORT2_H
#define TAYLORT2_H
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <complex>
#include "util.h"

/*! \file 
 * Header file for utilities
 */



//##########################################################################
//function declarations 
template <class T>
class TaylorT2
{
public:


virtual int construct_waveform(T *times, int length, std::complex<T> *hplus, std::complex<T> *hcross, source_parameters<T> *params);

virtual T THETA(T time, source_parameters<T> *params);

virtual T X(T theta , source_parameters<T> *params);

virtual T PHI(T theta , source_parameters<T> *params);

virtual T GAMMA(T x, source_parameters<T> *params);

virtual T M_ADM(T gamma, source_parameters<T> *params);

virtual void contract_H(T x, T psi, source_parameters<T> *params, std::complex<T> *hplus, std::complex<T> *hcross);


};
//###########################################################################
#endif
