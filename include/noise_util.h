#ifndef NOISE_UTIL_H
#define NOISE_UTIL_H
#include <string>
/*! \file
 */
void populate_noise(double *frequencies,std::string detector,double *noise_root,  int length=0);

double aLIGO_analytic(double f);

double Hanford_O1_fitted(double f);
#endif 
