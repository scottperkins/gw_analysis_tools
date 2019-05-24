#ifndef NOISE_UTIL_H
#define NOISE_UTIL_H
#include <string>
#include <complex>
/*! \file
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
const double V_radius = 6374824.24470673;//in meters
//const double V_elevation = 8;//in meters
const double V_elevation = 51.884;//in meters

//Earth radii at the poles and the equator (Earth is an ellispoid)
const double RE_polar =6357e3 ;//in meters
const double RE_equatorial = 6378e3 ;//in meters




void populate_noise(double *frequencies,std::string detector,double *noise_root,  int length=0);

double aLIGO_analytic(double f);

std::complex<double> Q(double theta, double phi, double iota);

double right_interferometer_cross(double theta, double phi);

double right_interferometer_plus(double theta, double phi);

double Hanford_O1_fitted(double f);

void celestial_horizon_transform(double RA, double DEC, double gps_time, std::string detector, double *phi, double *theta);

double DTOA(double theta1, double theta2, std::string detector1, std::string detector2);

double radius_at_lat(double latitude, double elevation);
#endif 
