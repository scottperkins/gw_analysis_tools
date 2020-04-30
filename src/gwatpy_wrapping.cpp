#ifndef GWATPY_WRAPPING_H
#define GWATPY_WRAPPING_H
#include "util.h"
#include "detector_util.h"


#ifdef __cplusplus
extern "C"
{
#endif

int calculate_chirpmass_py(double mass1, double mass2,double *out)
{
	*out = calculate_chirpmass(mass1,mass2);
	return 0;
}

void populate_noise_py(double *frequencies, char * detector, double *noise_root, int length, double integration_time){
	populate_noise(frequencies, std::string(detector), noise_root, length, integration_time);
	return;
}

#ifdef __cplusplus
}
#endif

#endif
