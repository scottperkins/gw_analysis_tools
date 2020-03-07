#ifndef DETECTOR_UTIL_H
#define DETECTOR_UTIL_H
#include <string>
#include <complex>
#include "util.h"
/*! \file
 * Header file for all detector-specific utilities
 */


/*Found LIGOs version see https://www.ligo.org/scientists/GW100916/detectors.txt*/
const double H_LAT =0.81079526383 ;//in rad
const double H_LONG =-2.08405676917 ;//in rad
//const double H_azimuth_offset = 5.65487724844;//in rad
//const double H_azimuth_offset = 4.08408092164;//in rad
const double H_azimuth_offset = 2.199;//in rad
const double H_radius = 6367299.93401105;
//const double H_elevation = 123.139;//in meters
const double H_elevation = 142.554;//in meters
const double H_location[3] = {-2.16141492636e+06,-3.83469517889e+06,4.60035022664e+06};

const double L_LAT =0.53342313506 ;//in rad
const double L_LONG =-1.58430937078 ;//in rad
//const double L_azimuth_offset = 4.40317772346;//in rad
//const double L_azimuth_offset = 2.83238139666;//in rad
const double L_azimuth_offset = 3.4557;//in rad
const double L_radius = 6372795.50144497;
//const double L_elevation = 13.1064;//in meters
const double L_elevation = -6.574;//in meters
const double L_location[3] = {-7.42760447238e+04,-5.49628371971e+06,3.22425701744e+06};

const double V_LAT =0.76151183984 ;//in rad
const double V_LONG =0.18333805213 ;//in rad
//const double V_azimuth_offset = 0.33916285222;//in rad
//const double V_azimuth_offset = 5.05155183261;//in rad
const double V_azimuth_offset = 1.239;//in rad
//const double V_radius = 6374824.24470673;//in meters
const double V_radius = 6368051.92301;//in meters
//const double V_elevation = 8;//in meters
const double V_elevation = 51.884;//in meters
const double V_location[3] = {4.54637409900e+06,8.42989697626e+05,4.37857696241e+06 };


const double K_LAT =0.6355068497 ;//in rad
const double K_LONG =2.396441015 ;//in rad
const double K_radius =6371060 ;//in meters
//const double K_azimuth_offset = ;//in rad
//######################################
//WRONG
const double K_azimuth_offset = 1.239;//in rad
//######################################
const double K_elevation = 414.181;//in meters
const double K_location[3] = {-3777336.024,3484898.411, 3765313.697};

const double I_LAT =0.2484185302005262;//in rad
const double I_LONG =1.3340133249409993 ;//in rad
const double I_radius =6376850 ;//in meters
//######################################
//WRONG
const double I_azimuth_offset = 1.239;//in rad
//######################################
//const double I_azimuth_offset = ;//in rad
const double I_elevation = 3.9624;//in meters
const double I_location[3] = {1450526.82294155,6011058.39047265,1558018.27884102 };

const double CE_LAT = 0.706553;//in rad
const double CE_LONG = -1.99869;//in rad
const double CE_radius =6370201 ;//in meters
const double CE_azimuth_offset = 0;//in rad
//const doublCE I_azimuth_offset = ;//in rad
const double CE_elevation = 1190.854;//in meters
const double CE_location[3] = {-2.01262e6,-4.41298e6,4.13995e6};

//Earth radii at the poles and the equator (Earth is an ellispoid)
const double RE_polar =6357e3 ;//in meters
const double RE_equatorial = 6378e3 ;//in meters

/*! Response Tensor for Hanford
 *
 * Calculated using arXiv:gr-qc/0008066, equation B6 with the table at the end -- see mathematica script
 * */
const double Hanford_D[3][3] = {{-0.392614, -0.0776134, -0.247389}, {-0.0776134, 0.319524, 0.227998}, {-0.247389, 0.227998, 0.07309}};

/*! Response Tensor for Livingston
 *
 * Calculated using arXiv:gr-qc/0008066, equation B6 with the table at the end -- see mathematica script
 */
const double Livingston_D[3][3] = {{0.411281, 0.14021, 0.247295}, {0.14021, -0.109006, -0.181616}, {0.247295, -0.181616, -0.302275}};

/*! Response Tensor for Virgo
 *
 * Calculated using arXiv:gr-qc/0008066, equation B6 with the table at the end -- see mathematica script
 */
const double Virgo_D[3][3] = {{0.243874, -0.0990838, -0.232576}, {-0.0990838, -0.447826, 
  0.187833}, {-0.232576, 0.187833, 0.203952}};

const double Kagra_D[3][3]={{-0.18599, 0.153167, -0.324951}, {0.153167, 
  0.349518, -0.170874}, {-0.324951, -0.170874, -0.163529}};

const double Indigo_D[3][3] ={{0.470823664959826499, -0.120908243530907122, 
  0.0279526040438164719}, {-0.120908243530907122, -0.001050025852002534, 0.115836952558600455}, {0.0279526040438164719, 
  0.115836952558600455, -0.469773644459096742}};

/*! Response Tensor for CE placed at Great Basin desert with the yarm pointing north and the x arm pointing east
 */
const double CE_D[3][3] = {{0.377622, -0.268333, -0.102451}, {-0.268333, -0.0883621, -0.224639}, {-0.102451, -0.224639, -0.289259}};


const int analytic_PSD_models_N = 6;
const std::string analytic_PSD_models[6] = {"aLIGO_analytic","Hanford_O1_fitted","LISA_SADC","LISA_SADC_CONF","LISA_CONF","LISA"};
const int interp_PSD_models_N = 20;
const std::string interp_PSD_models[20] = {"AdLIGOMidHigh","_AdLIGODesign","AdLIGODesign","AdLIGODesign_smoothed","AdLIGOAPlus","AdLIGOAPlus_smoothed","CE1","CE1_smoothed","CE2","CE2_smoothed","AdVIRGOPlus2_opt","AdVIRGOPlus2_opt_smoothed","AdVIRGOPlus2_pess","AdVIRGOPlus2_pess_smoothed","AdVIRGOPlus1","AdVIRGOPlus1_smoothed","KAGRA_opt","KAGRA_opt_smoothed","KAGRA_pess","KAGRA_pess_smoothed"};
void populate_noise(double *frequencies,std::string detector,double *noise_root,  int length, double integration_time=48);

double aLIGO_analytic(double f);

double LISA_analytic_SADC(double f);
double LISA_analytic(double f);
double LISA_POMS(double f);
double LISA_PACC(double f);
double LISA_SC(double f, double alpha, double beta, double kappa, double gamma, double fk);
void sort_LISA_SC_coeffs( double *alpha, double *beta, double  *kappa, double *gamma, double *fk,  double integration_time);

std::complex<double> Q(double theta, double phi, double iota);
std::complex<double> Q(double theta, double phi, double iota, double psi);

template<class T>
T right_interferometer_cross(T theta, T phi);

template<class T>
T right_interferometer_plus(T theta, T phi);

template<class T>
void right_interferometer(T *fplus, T *fcross, T theta, T phi, T psi);

double Hanford_O1_fitted(double f);

template<class T>
void celestial_horizon_transform(T RA, T DEC, double gps_time, std::string detector, T *phi, T *theta);

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

template<class T>
void detector_response_functions_equatorial(double D[3][3],
	T ra,
	T dec,
	T psi,
	double gmst,
	T *Fplus,
	T *Fcross
	);
template<class T>
void detector_response_functions_equatorial(std::string detector,
	T ra,
	T dec,
	T psi,
	double gmst,
	T *Fplus,
	T *Fcross
	);

template<class T>
void detector_response_functions_equatorial(std::string detector,
	T ra,
	T dec,
	T psi,
	double gmst,
	T *times,
	int length,
	T LISA_alpha0,
	T LISA_phi0,
	T theta_j_ecl,
	T phi_j_ecl,
	T *Fplus,
	T *Fcross
	);

double DTOA_DETECTOR(double RA, 
	double DEC, 
	double GMST_rad,
	std::string detector1, 
	std::string detector2 );
double DTOA_earth_centered_coord(double RA, 
	double DEC, 
	double GMST_rad,
	const double *loc1, 
	const double *loc2 );
template<class T>
T LISA_response_plus_time( T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T t);
template<class T>
T LISA_response_cross_time( T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T t);
template<class T>
T LISA_response_plus( source_parameters<T> *params,T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T f);
template<class T>
T LISA_response_cross( source_parameters<T> *params,T theta_s, T phi_s, T theta_j, T phi_j, T alpha_0, T phi_0, T f);

double p_single_detector(double omega, int samples);
double p_N_detector(double omega, int samples,int N_detectors,std::string *detectors,int rand_seed);
double p_single_detector_fit(double omega);
double p_triple_detector_fit(double omega);
double pdet_triple_detector_fit(double rho_thresh,double rho_opt);
double p_triple_detector_interp(double omega);
double p_single_detector_interp(double omega);
#endif 
