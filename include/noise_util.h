#ifndef NOISE_UTIL_H
#define NOISE_UTIL_H
#include <string>
#include <complex>
/*! \file
 */
void populate_noise(double *frequencies,std::string detector,double *noise_root,  int length=0);

double aLIGO_analytic(double f);

std::complex<double> Q(double theta, double phi, double iota);

double right_interferometer_cross(double theta, double phi);

double right_interferometer_plus(double theta, double phi);

double Hanford_O1_fitted(double f);
#endif 
