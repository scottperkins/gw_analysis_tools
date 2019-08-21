#include "detector_util.h"
#include "util.h"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

/*! \file
 * Routines to construct noise curves for various detectors and for detector specific utilities for response functions and coordinate transformations
 */

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

/*! \brief Analytic function approximating the PSD for aLIGO 
 *
 * CITE (Will?)
 */
double aLIGO_analytic(double f)
{
	double S= 3e-48;
	double fknee = 70.;
	double x = fknee/f;
	double x4 = x*x*x*x;
	return sqrt( S * (x4 + 2 + 2*x*x)/5 );
}

/*! \brief Numerically fit PSD to the Hanford Detector's O1
 *
 * CITE (Yunes?)
 */
double Hanford_O1_fitted(double f)
{
	double avec[7]={47.8466,-92.1896,35.9273,-7.61447,0.916742,-0.0588089,0.00156345};
	double S0 = .8464;
	double x = log(f);
	return sqrt(S0) * exp( avec[0] + avec[1]*x + avec[2]*x*x +
            avec[3] * x*x*x + avec[4]*x*x*x*x + avec[5]*x*x*x*x*x + avec[6]*x*x*x*x*x*x);
}

/*! \brief Utility for the overall amplitude and phase shift for spin-aligned systems
 *
 * For spin aligned, all the extrinsic parameters have the effect of an overall amplitude modulation and phase shift
 */
std::complex<double> Q(double theta, double phi, double iota, double psi)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);
	double sp2 = sin(2.*phi);
	double ci = cos(iota);
	double cpsi2 = cos(2*psi);
	double spsi2 = sin(2*psi);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	double Fcross = ct * sp2;
	//std::complex<double> Q = (1+ci*ci)/2. *Fplus + std::complex<double>(0,Fcross*ci);
	std::complex<double> Q= (1+ci*ci)/2. *(Fplus*cpsi2 - Fcross*spsi2) + std::complex<double>(0,(Fcross *cpsi2 + Fplus*spsi2)*ci);
	return Q;
}

/*! \brief Utility for the overall amplitude and phase shift for spin-aligned systems
 *
 * For spin aligned, all the extrinsic parameters have the effect of an overall amplitude modulation and phase shift
 */
std::complex<double> Q(double theta, double phi, double iota)
{
	std::complex<double> Qout= Q(theta, phi, iota, 0.);
	return Qout;
}

/*! \brief Response function of a 90 deg interferometer for plus polarization
 *
 * Theta and phi are local, horizontal coordinates relative to the detector
 */
double right_interferometer_plus(double theta, double phi)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	return Fplus;
}

/*! \brief Response function of a 90 deg interferometer for cross polarization
 *
 * Theta and phi are local, horizontal coordinates relative to the detector
 */
double right_interferometer_cross(double theta, double phi)
{
	double ct = cos(theta);
	double sp2 = sin(2.*phi);
	double Fcross = ct * sp2;
	return Fcross;
}

/*! \brief Transform from celestial coordinates to local horizontal coords
 *
 * (RA,DEC) -> (altitude, azimuth)
 *
 * Need gps_time of transformation, as the horizontal coords change in time
 *
 * detector is used to specify the lat and long of the local frame
 */
void celestial_horizon_transform(double RA, /**< in RAD*/
		double DEC, /**< in RAD*/
		double gps_time, 
		std::string detector, 
		double *phi, /**< in RAD*/
		double *theta /**< in RAD*/
		)
{
	double LAT, LONG, azimuth_offset;
	if (detector =="Hanford" || detector == "hanford")
	{
		LAT =  H_LAT;		
		LONG =  H_LONG;		
		azimuth_offset = H_azimuth_offset;
		
	}
	else if (detector =="Livingston" || detector == "livingston")
	{
		LAT =  L_LAT;		
		LONG =  L_LONG;		
		azimuth_offset = L_azimuth_offset;

	}
	else if (detector =="Virgo" || detector == "virgo")
	{
		LAT =  V_LAT;		
		LONG =  V_LONG;		
		azimuth_offset = V_azimuth_offset;

	}
	else {
		std::cout<<"Invalid detector"<<std::endl;
		exit(1);
	}
	celestial_horizon_transform(RA,DEC, gps_time, LONG, LAT, phi, theta);
	*phi += azimuth_offset;
	if(*phi>2*M_PI) *phi -=2*M_PI;
}

/*! \brief Numerical derivative of the transformation
 *
 * Planned for use in Fisher calculations, but not currently implemented anywhere
 */
void derivative_celestial_horizon_transform(double RA, /**< in RAD*/
		double DEC, /**< in RAD*/
		double gps_time, 
		std::string detector, 
		double *dphi_dRA, 
		double *dtheta_dRA, 
		double *dphi_dDEC, 
		double *dtheta_dDEC 
		)
{
	double phip, phim, thetap, thetam;
	double epsilon = 1e-5;

	celestial_horizon_transform(RA+epsilon, DEC,gps_time,detector,&phip,&thetap);
	celestial_horizon_transform(RA-epsilon, DEC,gps_time,detector,&phim,&thetam);
	*dtheta_dRA = (thetap-thetam)/(2*epsilon);
	*dphi_dRA = (phip-phim)/(2*epsilon);

	celestial_horizon_transform(RA, DEC+epsilon,gps_time,detector,&phip,&thetap);
	celestial_horizon_transform(RA, DEC-epsilon,gps_time,detector,&phim,&thetam);
	*dtheta_dDEC = (thetap-thetam)/(2*epsilon);
	*dphi_dDEC = (phip-phim)/(2*epsilon);
}

/*! \brief calculate difference in time of arrival (DTOA) for a given source location and 2 different detectors
 */
double DTOA(double theta1, /**< spherical polar angle for detector 1 in RAD*/
	double theta2, /**<spherical polar angle for detector 2 in RAD*/ 
	std::string detector1, /**< name of detector one*/
	std::string detector2 /**<name of detector two*/
	)
{
	double R1, R2;
	//detector one
	if(detector1 == "Hanford" || detector1 == "hanford")
	{
		R1 = H_radius;	
	}
	else if(detector1 == "Livingston" || detector1 == "livingston")
	{
		R1 = L_radius;	
	}
	else if(detector1 == "Virgo" || detector1 == "virgo")
	{
		R1 = V_radius;	
	}
	//detector 2
	if(detector2 == "Hanford" || detector2 == "hanford")
	{
		R2 = H_radius;	
	}
	else if(detector2 == "Livingston" || detector2 == "livingston")
	{
		R2 = L_radius;	
	}
	else if(detector2 == "Virgo" || detector2 == "virgo")
	{
		R2 = V_radius;	
	}
	return (cos(theta1)*R1 - cos(theta2)*R2)/c;
}

/*! /brief Analytic approximation of the radius from the center of earth to a given location
 *
 * Just the raidus as a function of angles, modelling an oblate spheroid
 */
double radius_at_lat(double latitude, /**< latitude in degrees*/
			double elevation /**<elevation in meters*/
			)
{
	double numerator = pow(RE_equatorial*RE_equatorial * cos(latitude*M_PI/180.) ,2 ) 
			+ pow(RE_polar*RE_polar * sin(latitude*M_PI/180. ), 2);
	double denominator = pow(RE_equatorial * cos(latitude*M_PI/180.) ,2 ) 
			+ pow(RE_polar * sin(latitude*M_PI/180. ), 2);
	return sqrt(numerator/denominator) + elevation;
}

/*! \brief Calculates the response coefficients for a detector with response tensor D for a source at RA, Dec, and psi
 *
 * Taken from LALSuite
 *
 * The response tensor for each of the operational detectors is precomputed in detector_util.h, but to create a new tensor, follow the outline in Anderson et al 36 PRD 63 042003 (2001) Appendix B
 */
void detector_response_functions_equatorial(double D[3][3],/**< Detector Response tensor (3x3)*/
	double ra,/**<Right ascension in rad*/
	double dec,/**<Declination in rad*/
	double psi,/**< polarization angle in rad*/
	double gmst,/**<Greenwich mean sidereal time (rad)*/
	double *Fplus,/**< [out] Fplus response coefficient*/
	double *Fcross	/**<[out] Fcross response coefficient*/
	)
{
	int i;
	double X[3];
	double Y[3];
	
	/* Greenwich hour angle of source (radians). */
	const double gha = gmst - ra;
	
	/* pre-compute trig functions */
	const double cosgha = cos(gha);
	const double singha = sin(gha);
	const double cosdec = cos(dec);
	const double sindec = sin(dec);
	const double cospsi = cos(psi);
	const double sinpsi = sin(psi);
	
	/* Eq. (B4) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
	X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
	X[2] =  sinpsi * cosdec;
	
	/* Eq. (B5) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	Y[0] =  sinpsi * singha - cospsi * cosgha * sindec;
	Y[1] =  sinpsi * cosgha + cospsi * singha * sindec;
	Y[2] =  cospsi * cosdec;
	
	/* Now compute Eq. (B7) of [ABCF] for each polarization state, i.e.,
	 * with s+=1 and sx=0 to get F+, with s+=0 and sx=1 to get Fx */
	*Fplus = *Fcross = 0.0;
	for(i = 0; i < 3; i++) {
	        const double DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
	        const double DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
	        *Fplus  += X[i] * DX - Y[i] * DY;
	        *Fcross += X[i] * DY + Y[i] * DX;
	}

}
/*! \brief Same as the other function, but for active detectors
 */
void detector_response_functions_equatorial(std::string detector,/**< Detector */
	double ra,/**<Right ascension in rad*/
	double dec,/**<Declination in rad*/
	double psi,/**< polarization angle in rad*/
	double gmst,/**<Greenwich mean sidereal time (rad)*/
	double *Fplus,/**< [out] Fplus response coefficient*/
	double *Fcross	/**<[out] Fcross response coefficient*/
	)
{
	double responseM[3][3];
	if(detector =="Hanford" || detector=="hanford"){
		for(int i =0; i<3; i++){
			for(int j =0 ;j<3; j++){
				responseM[i][j] = Hanford_D[i][j];
			}
		}
	}
	else if(detector =="Livingston" || detector=="livingston"){
		for(int i =0; i<3; i++){
			for(int j =0 ;j<3; j++){
				responseM[i][j] = Livingston_D[i][j];
			}
		}
	}
	else if(detector =="Virgo" || detector=="virgo"){
		for(int i =0; i<3; i++){
			for(int j =0 ;j<3; j++){
				responseM[i][j] = Virgo_D[i][j];
			}
		}
	}
	else{
		std::cout<<"ERROR -- unsupported detector"<<std::endl;
		exit(1);
	}
	detector_response_functions_equatorial(responseM, ra, dec, psi, gmst, Fplus, Fcross);
}
