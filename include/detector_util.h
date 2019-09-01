#ifndef DETECTOR_UTIL_H
#define DETECTOR_UTIL_H
#include <string>
#include <complex>
/*! \file
 * Header file for all detector-specific utilities
 */

/*Detector Orientation  variables*/
/*https://www.ligo.org/scientists/GW100916/GW100916-geometry.html*/
/*Also see https://www.ligo.org/scientists/GW100916/detectors.txt*/
/* All in degrees*/
/*East is (-)*/
//const double H_LAT =46+ 27*60 + 19*3600 ;
//const double H_LONG =119+24*60 + 28*3600 ;
//const double H_azimuth_offset = 36 + 90;
//const double H_radius = 6367299.93401105;
//const double H_elevation = 123.139;//in meters
//const double L_LAT =30+33*60 + 46*3600 ;
//const double L_LONG =90+46*60 + 27*3600 ;
//const double L_azimuth_offset = 180 + 18;
//const double L_radius = 6372795.50144497;
//const double L_elevation = 13.1064;//in meters
//const double V_LAT =43+37*60 + 53*2600 ;
//const double V_LONG =-(10+30*60+16*3600) ;
//const double V_azimuth_offset = 71;
//const double V_radius = 6374824.24470673;//in meters
//const double V_elevation = 8;//in meters

/*Found LIGOs version see https://www.ligo.org/scientists/GW100916/detectors.txt*/
const double H_LAT =0.81079526383 ;//in rad
const double H_LONG =-2.08405676917 ;//in rad
//const double H_azimuth_offset = 5.65487724844;//in rad
//const double H_azimuth_offset = 4.08408092164;//in rad
const double H_azimuth_offset = 2.199;//in rad
const double H_radius = 6367299.93401105;
//const double H_elevation = 123.139;//in meters
const double H_elevation = 142.554;//in meters

const double L_LAT =0.53342313506 ;//in rad
const double L_LONG =-1.58430937078 ;//in rad
//const double L_azimuth_offset = 4.40317772346;//in rad
//const double L_azimuth_offset = 2.83238139666;//in rad
const double L_azimuth_offset = 3.4557;//in rad
const double L_radius = 6372795.50144497;
//const double L_elevation = 13.1064;//in meters
const double L_elevation = -6.574;//in meters

const double V_LAT =0.76151183984 ;//in rad
const double V_LONG =0.18333805213 ;//in rad
//const double V_azimuth_offset = 0.33916285222;//in rad
//const double V_azimuth_offset = 5.05155183261;//in rad
const double V_azimuth_offset = 1.239;//in rad
//const double V_radius = 6374824.24470673;//in meters
const double V_radius = 6368051.92301;//in meters
//const double V_elevation = 8;//in meters
const double V_elevation = 51.884;//in meters

//Earth radii at the poles and the equator (Earth is an ellispoid)
const double RE_polar =6357e3 ;//in meters
const double RE_equatorial = 6378e3 ;//in meters

const double Hanford_D[3][3] = {{-0.392632, -0.0776099, -0.247384}, {-0.0776099, 0.319499, 
  0.227988}, {-0.247384, 0.227988, 0.0730968}};

const double Livingston_D[3][3] = {{0.411318, 0.14021, 
  0.247279}, {0.14021, -0.108998, -0.181597}, {0.247279, -0.181597, 
-0.302236}};

const double Virgo_D[3][3] = {{0.243903, -0.0990959, -0.232603}, {-0.0990959, -0.447841, 
  0.187841}, {-0.232603, 0.187841, 0.203979}};




void populate_noise(double *frequencies,std::string detector,double *noise_root,  int length=0);

double aLIGO_analytic(double f);

std::complex<double> Q(double theta, double phi, double iota);
std::complex<double> Q(double theta, double phi, double iota, double psi);

double right_interferometer_cross(double theta, double phi);

double right_interferometer_plus(double theta, double phi);

double Hanford_O1_fitted(double f);

void celestial_horizon_transform(double RA, double DEC, double gps_time, std::string detector, double *phi, double *theta);

void derivative_celestial_horizon_transform(double RA, 
		double DEC, 
		double gps_time, 
		std::string detector, 
		double *dphi_dRA, 
		double *dtheta_dRA, 
		double *dphi_dDEC, 
		double *dtheta_dDEC 
		);

double DTOA(double theta1, double theta2, std::string detector1, std::string detector2);

double radius_at_lat(double latitude, double elevation);

void detector_response_functions_equatorial(double D[3][3],
	double ra,
	double dec,
	double psi,
	double gmst,
	double *Fplus,
	double *Fcross
	);
void detector_response_functions_equatorial(std::string detector,
	double ra,
	double dec,
	double psi,
	double gmst,
	double *Fplus,
	double *Fcross
	);

template<class T>
T LISA_response_plus( T theta_s, T phi_s, T theta_l, T phi_l, T alpha_0, T phi_0);
template<class T>
T LISA_response_cross( T theta_s, T phi_s, T theta_l, T phi_l, T alpha_0, T phi_0);

#endif 
