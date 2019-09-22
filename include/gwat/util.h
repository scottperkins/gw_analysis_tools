#ifndef UTIL_H
#define UTIL_H
//#include "general_parameter_structures.h"
#include <string>
#include <complex>
#include "adolc/adouble.h"
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
//#define adouble double
/*! \file 
 *General utilities (functions and structures) independent of modelling method
 */

/*! Euler number*/
const double gamma_E = 0.5772156649015328606065120900824024310421;
/*!Speed of light m/s*/
const double c = 299792458.;
/*!Gravitational constant in m**3/(s**2 SolMass)*/
const double G =6.674e-11*(1.98855e30);
/*! G/c**3 seconds per solar mass*/
//const double MSOL_SEC =492549095.e-14; 
//const double MSOL_SEC =492549095.e-14; 
const double MSOL_SEC =4.925491025543575903411922162094833998e-6 ;
/*!consts.kpc.to('m')*1000/c Mpc in sec*/
//const double MPC_SEC = 3085677581.e13/c; 
const double MPC_SEC = 3.085677581491367278913937957796471611e22/c; 
/*!1 year in seconds*/
const double T_year = 31557600.;


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
//GSL versions of the constants, but G seems off..
//const double gamma_E = M_EULER;
//const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
//const double G =GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;
//const double MSOL_SEC = GSL_CONST_MKSA_SOLAR_MASS*(GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT/(c*c*c));
//const double MPC_SEC = GSL_CONST_MKSA_PARSEC*1e6/c; 

struct fftw_outline
{
	fftw_complex *in, *out;
	fftw_plan p;
};

/*!\struct
 *
 * Container for set of spin weighted spherical harmonics (all spins are -2)
 */
template <class T>
struct sph_harm
{
	std::complex<T> Y22;
	std::complex<T> Y21;
	std::complex<T> Y20;
	std::complex<T> Y2m1;
	std::complex<T> Y2m2;
};

/*!\struct
 * \brief Structure for interfacing with the libraries 
 * 
 * Structure to interface with the libraries - Units are in solar masses and mpc 
 *  contains the generation parameters including source parameters and theory parameters
 *
 *  *NOTE* not all the members of this structure need to be assigned for usage. In fact, some are reduntant. It's up to the user to determine what fields require an assignment. (Sorry)
 */

template<class T>
class gen_params_base
{
public:	
	/*!mass of the larger body in Solar Masses*/
	T mass1;
	/*!mass of the smaller body in Solar Masses*/
	T mass2;
	/*!Luminosity distance to the source*/
	T Luminosity_Distance;
	/*!Spin vector of the larger mass [Sx,Sy,Sz]*/
	T spin1[3];
	/*!Spin vector of the smaller mass [Sx,Sy,Sz]*/
	T spin2[3];
	/*!coalescence phase of the binary*/
	T phic=0;
	/*!coalescence time of the binary*/
	T tc=0;
	/*!ppE b parameter (power of the frequency) - vector for multiple modifications*/
	int *bppe;
	/*!ppE coefficient for the phase modification - vector for multiple modifications*/
	T *betappe;
	/*!Number of phase modificatinos*/
	int Nmod;
	/*!*angle between angular momentum and the total momentum */
	T incl_angle;
	/*! spherical angles for the source location relative to the detector*/
	T theta;
	T phi;
	/*! Equatorial coordinates of source*/
	T RA;
	T DEC;
	/*! Greenwich Mean Sidereal time (for detector orientation - start of data*/
	double gmst;
	//Polarization angle
	T psi =0 ;
	/*! BOOL flag for early termination of NS binaries*/
	bool NSflag;

	/*! Reference frequency for PhenomPv2*/
	T f_ref=0;
	
	T phiRef=0;

	//###################################################
	//Either define all these parameters for Pv2, or define
	//all the source frame parameters above
	/*!thetaJ -- optional domain is [0,M_PI]  */
	T thetaJN = -10;

	T alpha0 = 0;


	T chip = 0;

	T chi1_l = 0;
	T chi2_l = 0;
	//###################################################
	T zeta_polariz =0;
	T phi_aligned = 0;

	T chil = 0;


	bool sky_average;
	
	gsl_spline *Z_DL_spline_ptr=NULL;

	gsl_interp_accel *Z_DL_accel_ptr=NULL;
		
	std::string cosmology="PLANCK15";
	
	bool shift_time = true;
	
};

/*! \brief convience wrapper for the gen_params_base class
 *
 * If using the code in the intended way, this is all the user should ever have to use. Just allows the user to drop the template parameter
 *
 * Also implemented for backwards compatibility with previous versions of the code
 */
class gen_params:public gen_params_base<double>
{

};
/*!\brief To speed up calculations within the for loops, we pre-calculate reoccuring powers of M*F and Pi, since the pow() function is prohibatively slow
 *
 * Powers of PI are initialized once, and powers of MF need to be calculated once per for loop (if in the inspiral portion).
 *
 * use the functions precalc_powers_ins_amp, precalc_powers_ins_phase, precalc_powers_pi to initialize
 */
template <class T>
struct useful_powers
{
	T MFthird;
	T MFsixth;
	T MF7sixth;
	T MF2third;
	T MF4third;
	T MF5third;
	T MFsquare;
	T MF7third;
	T MF8third;
	T MFcube;
	T MFminus_5third;
	T MF3fourth;
	double PIsquare;
	double PIcube;
	double PIthird;
	double PI2third;
	double PI4third;
	double PI5third;
	double PI7third;
	double PIminus_5third;
};

/*!\struct 
 * \brief For internal data transfers
 *
 * Structure to facililate parameter tranfers - All dimensionful quantities are in seconds
 */
template <class T>
struct source_parameters
{
	/*! mass of the larger component*/
	T mass1;
	/*! mass of the smaller component*/
	T mass2;
	/*! Total mass*/
	T M;
	/*Mass ratio*/	
	T q;
	/*! z-Spin component of the larger body*/
	T spin1z;
	/*! z-Spin component of the smaller body*/
	T spin2z;
	/*! x-Spin component of the larger body*/
	T spin1x;
	/*! x-Spin component of the smaller body*/
	T spin2x;
	/*! y-Spin component of the larger body*/
	T spin1y;
	/*! y-Spin component of the smaller body*/
	T spin2y;
	/*!Chirp mass of the binary*/
	T chirpmass;
	/*!Symmetric mass ratio*/
	T eta;
	/*!Symmetric spin combination*/
	T chi_s;
	/*!Antisymmetric spin combination*/
	T chi_a;
	/*!Effective spin */
	T chi_eff;
	/*! PN spin */
	T chi_pn;
	/*!Luminoisity Distance*/
	T DL;
	/*! Delta mass comibination*/
	T delta_mass;
	/*!Ringdown frequency after merger*/
	T fRD;
	/*!Dampening frequency after merger*/
	T fdamp;
	/*! Transition Frequency 1 for the amplitude*/	
	T f1;
	/*! Transition Frequency 2 for the amplitude*/	
	T f3;
	/*! Transition frequency 1 for the phase*/
	T f1_phase;
	/*! Transition frequency 2 for the phase*/
	T f2_phase;
	/*! Coalescence phase*/
	T phic;
	/*! Coalescence time*/
	T tc;
	/*overall amplitude factor*/	
	T A0;
	
	//############################################
	//PhenomP
	T s ;

	T chil;
	
	T chip;

	T f_ref=0;
	
	T phi_aligned;

	T incl_angle;
	
	T phiRef;

	T alpha0;

	T thetaJN;

	T zeta_polariz;

	//######### ppE parameters ##############
	/*Beta factor for ppE formalism*/
	T *betappe;
	/*b power for ppE formalism*/
	int *bppe;

	/*! Number of modifications to phase*/	
	int Nmod;
	
	//Spherical polar angles for the sky location relative to the detector in question
	T phi;

	T theta;
	
	T SP;

	T SL;

	bool sky_average;

	/*! Boolean -- shift time to 0 before shifting to tc or not*/
	bool shift_time = true;

	gsl_spline *Z_DL_spline_ptr=NULL;

	gsl_interp_accel *Z_DL_accel_ptr=NULL;
	
	std::string cosmology;

static source_parameters<T> populate_source_parameters(gen_params_base<T> *param_in);
static source_parameters<T> populate_source_parameters_old(
			T mass1, 
			T mass2, 
			T Luminosity_Distance, 
			T *spin1,
			T *spin2, 
			T phi_c,
			T t_c, 
			bool sky_average) ;
};
template<class T>
T A0_from_DL(T chirpmass, T DL, bool sky_average);

template<class T>
T DL_from_A0(T chirpmass, T A0, bool sky_average);

void initiate_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr);
void free_LumD_Z_interp(gsl_interp_accel **Z_DL_accel_ptr, gsl_spline **Z_DL_spline_ptr);
adouble Z_from_DL_interp(adouble DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr);
double Z_from_DL_interp(double DL,gsl_interp_accel *Z_DL_accel_ptr, gsl_spline *Z_DL_spline_ptr);

double Z_from_DL(double DL, std::string cosmology);
double DL_from_Z(double Z, std::string cosmology);
double cosmology_interpolation_function(double x, double *coeffs, int interp_degree);
double cosmology_lookup(std::string cosmology);

template<class T>
bool check_list(T j, T *list, int length);
void rm_fisher_dim(double **input,int full_dim, double **output,  int reduced_dim, int *removed_dims);

template<class T>
void gsl_LU_matrix_invert(T **input, T **inverse, int dim);

int gsl_cholesky_matrix_invert(double **input, double **inverse, int dim);

int normalized_gsl_cholesky_matrix_invert(double **input, double **inverse, int dim);

adouble Z_from_DL(adouble DL, std::string cosmology);
adouble DL_from_Z(adouble Z, std::string cosmology);
adouble cosmology_interpolation_function(adouble x, double *coeffs, int interp_degree);

void printProgress (double percentage);

void allocate_FFTW_mem_forward(fftw_outline *plan,int length);
void allocate_FFTW_mem_reverse(fftw_outline *plan,int length);
void deallocate_FFTW_mem(fftw_outline *plan);

double** allocate_2D_array( int dim1, int dim2);
int** allocate_2D_array_int( int dim1, int dim2);
void deallocate_2D_array(double **array, int dim1, int dim2);
void deallocate_2D_array(int **array, int dim1, int dim2);
double*** allocate_3D_array( int dim1, int dim2, int dim3);
int*** allocate_3D_array_int( int dim1, int dim2, int dim3);

void deallocate_3D_array(double ***array, int dim1, int dim2, int dim3);
void deallocate_3D_array(int ***array, int dim1, int dim2, int dim3);

void read_file(std::string filename,double **output, int rows, int cols );
void read_file(std::string filename, double *output );

void read_LOSC_data_file(std::string filename, 
			double *output,
			double *data_start_time,
			double *duration,
			double *fs);

void read_LOSC_PSD_file(std::string filename, 
			double **output,
			int rows,
			int cols);

void allocate_LOSC_data(std::string *data_files, 
			std::string psd_file, 
			int num_detectors,
			int psd_length,
			int data_file_length,
			double trigger_time,
			std::complex<double> **data,
			double **psds,
			double **freqs
			);

void free_LOSC_data(std::complex<double> **data,
		double **psds,
		double **freqs,
		int num_detectors,
		int length
		);

void tukey_window(double *window,
		int length,
		double alpha);

void write_file(std::string filename, double **input, int rows, int cols);
void write_file(std::string filename, double *input, int length);

double calculate_eta(double mass1, double mass2);
adouble calculate_eta(adouble mass1, adouble mass2);

double calculate_chirpmass(double mass1, double mass2);
adouble calculate_chirpmass(adouble mass1, adouble mass2);

double calculate_mass1(double chirpmass, double eta);
adouble calculate_mass1(adouble chirpmass, adouble eta);
	
double calculate_mass2(double chirpmass, double eta);
adouble calculate_mass2(adouble chirpmass, adouble eta);

template<class T>
void celestial_horizon_transform(T RA, T DEC, double gps_time, T LONG, T LAT,
				T *phi, T *theta);

template<class T>
T gps_to_GMST(T gps_time);
template<class T>
T gps_to_JD(T gps_time);


template<class T>
void transform_cart_sph(T *cartvec, T *sphvec);

template<class T>
void transform_sph_cart(T *sphvec, T *cartvec);

/*!\brief Trapezoidal sum rule to approximate discrete integral - Uniform spacing
 *
 * This version is faster than the general version, as it has half the function calls
 *
 * Something may be wrong with this function - had an overall offset for real data that was
 * fixed by using the simpsons rule - not sure if this was because of a boost in accuracy or 
 * because something is off with the trapezoidal sum
 */
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum_uniform(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += 2*integrand[i];
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/2;
}

/*!\brief Trapezoidal sum rule to approximate discrete integral - Non-Uniform spacing
 *
 * This version is slower than the uniform version, but will handle non-uniform spacing
 */
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand);
template <class T>
T trapezoidal_sum(double *delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
		sum += (integrand[i]+integrand[i-1])*(delta_x[i]-delta_x[i-1])/2;
	return sum;
}

/*!\brief Simpsons sum rule to approximate discrete integral - Uniform spacing
 *
 * More accurate than the trapezoidal rule, but must be uniform
 */
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand);
template <class T>
T simpsons_sum(double delta_x, int length, T *integrand)
{
	T sum=0;
	for (int i =1;i<length-1;i++)
	{
		if( i % 2 == 0)
			sum += 2*integrand[i];
		else
			sum+=4*integrand[i];
	}
	sum += integrand[0]+integrand[length-1];
	return sum*delta_x/3;
}
long factorial(long num);
double pow_int(double base, int power);
adouble pow_int(adouble base, int power);

template <class T>
std::complex<T> cpolar(T mag, T phase);

template <class T>
std::complex<T> XLALSpinWeightedSphericalHarmonic(
                                   T theta,  /**< polar angle (rad) */
                                   T phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    );
double cbrt_internal(double base);
adouble cbrt_internal(adouble base);
#endif
