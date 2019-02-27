#include "noise_util.h"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

/*! \file
 * Routines to construct noise curves for various detectors*/

/*! \brief Function to populate the squareroot of the noise curve for various detectors
 *
 * If frequencies are left as NULL, standard frequency spacing is applied and the frequencies are returned, in which case the frequencies argument becomes an output array
 *
 * Detector names must be spelled exactly
 *
 * Detectors include: aLIGO_analytic, Hanford_O1_fitted
 */
void populate_noise(double *frequencies, /**< double array of frquencies (NULL)*/
		std::string detector, /**< String to designate the detector noise curve to be used */
		double *noise_root, /**< ouptput double array for the square root of the PSD of the noise of the specified detector*/
		int length/**< integer length of the output and input arrays*/
		)
{
	if(detector == "aLIGO_analytic")
	{
		if(frequencies == NULL)
		{
			int len = 1000;
			for (int i =0; i<len;i++)
			{
				frequencies[i] = 10.+ i;
				noise_root[i] = aLIGO_analytic(frequencies[i]);
			}
		}
		
		else
		{
			for (int i =0; i<length;i++)
			{
				noise_root[i] = aLIGO_analytic(frequencies[i]);
			}
		}
	}
	if(detector == "Hanford_O1_fitted")
	{
		if(frequencies == NULL)
		{
			int len = 1000;
			for (int i =0; i<len;i++)
			{
				frequencies[i] = 10.+ i;
				noise_root[i] = Hanford_O1_fitted(frequencies[i]);
			}
		}
		
		else
		{
			for (int i =0; i<length;i++)
			{
				noise_root[i] = Hanford_O1_fitted(frequencies[i]);
			}
		}
	}
			
		
}


double aLIGO_analytic(double f)
{
	double S= 3e-48;
	double fknee = 70.;
	double x = fknee/f;
	double x4 = x*x*x*x;
	return sqrt( S * (x4 + 2 + 2*x*x)/5 );
}

double Hanford_O1_fitted(double f)
{
	double avec[7]={47.8466,-92.1896,35.9273,-7.61447,0.916742,-0.0588089,0.00156345};
	double S0 = .8464;
	double x = log(f);
	return sqrt(S0) * exp( avec[0] + avec[1]*x + avec[2]*x*x +
            avec[3] * x*x*x + avec[4]*x*x*x*x + avec[5]*x*x*x*x*x + avec[6]*x*x*x*x*x*x);
}
