#ifndef IMRPHENOMD_H
#define IMRPHENOMD_H
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <complex>
#include "util.h"

/*! \file 
 * Header file for utilities
 */

/*! \struct
 * Structure to facilitate IMRPhenomD parameter transfers
 */
template <class T>
struct lambda_parameters
{
	T rho[4];
	T v2;
	T gamma[4];
	T sigma[5];
	T beta[5];
	T alpha[7];
};


//##########################################################################
//function declarations 
template <class T>
class IMRPhenomD
{
public:

virtual void fisher_calculation_sky_averaged(double *frequency, 
			int length, 
			//double *parameters,
			gen_params *parameters,
			double **amplitude_deriv, 
			double **phase_deriv, 
			double *amplitude, 
			int *amp_tapes, 
			int *phase_tapes
			);
virtual void change_parameter_basis(T *old_param,
					T *new_param,
					bool sky_average
					);
virtual void construct_amplitude_derivative(double *frequencies, 
				int length,
				int dimension, 
				double **amplitude_derivative,
				source_parameters<double> *input_params,
				int *tapes=NULL
				);
virtual void construct_phase_derivative(double *frequencies, 
				int length,
				int dimension, 
				double **phase_derivative,
				source_parameters<double> *input_params,
				int *tapes=NULL
				);
virtual void amplitude_tape(source_parameters<double> *input_params, int *tape);

virtual void phase_tape(source_parameters<double> *input_params, int *tape);


virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);

virtual std::complex<T> construct_waveform(T frequency, 				
			source_parameters<T> *params );

virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);

//virtual T construct_amplitude(T frequency,  source_parameters<T> *params);

virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);

virtual T build_amp(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *amp_coeff, 
		T *deltas);

virtual T build_phase(T f, 
		lambda_parameters<T> *lambda, 
		source_parameters<T> *params, 
		useful_powers<T> *pows,
		T *phase_coeff);

virtual T assign_lambda_param_element(source_parameters<T> *source_param,int i);

virtual void assign_lambda_param(source_parameters<T> *source_param, lambda_parameters<T> *lambda);

virtual void precalc_powers_ins(T f, T M, useful_powers<T> *Mf_pows);

virtual void precalc_powers_PI( useful_powers<T> *PI_pows);

virtual void precalc_powers_ins_phase(T f, T M, useful_powers<T> *Mf_pows);

virtual void precalc_powers_ins_amp(T f, T M, useful_powers<T> *Mf_pows);

virtual void assign_pn_amplitude_coeff(source_parameters<T> *source_param, T *coeff);

virtual void assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff);

virtual void assign_nonstatic_pn_phase_coeff(source_parameters<T> *source_param, 
					T *coeff, 
					T f);

virtual void assign_nonstatic_pn_phase_coeff_deriv(source_parameters<T> *source_param, 
					T *Dcoeff, 
					T f);
virtual void post_merger_variables(source_parameters<T> *source_param);

virtual T fpeak(source_parameters<T> *params, lambda_parameters<T> *lambda);

virtual T amp_ins(T f, source_parameters<T> *param, T *pn_coeff, 
			lambda_parameters<T> *lambda,useful_powers<T> *pow);

virtual T Damp_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

virtual T phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *pow);

virtual T Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda);

virtual T amp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T Damp_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T amp_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda, T *deltas);

virtual T phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda );

virtual T Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual void phase_connection_coefficients(source_parameters<T> *param, 
				lambda_parameters<T> *lambda, 
				T *pn_coeffs);

virtual T calculate_beta1(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs);

virtual T calculate_beta0(source_parameters<T> *param, lambda_parameters<T> *lambda, T *pn_coeffs);

virtual T calculate_alpha1(source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual T calculate_alpha0(source_parameters<T> *param, lambda_parameters<T> *lambda);

virtual void amp_connection_coeffs(source_parameters<T> *param, 
			lambda_parameters<T> *lambda, 
			T *pn_coeffs, 
			T *coeffs);

virtual T calculate_delta_parameter_0(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_1(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_2(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_3(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
virtual T calculate_delta_parameter_4(T f1,T f2,T f3,T v1,
				T v2,T v3,T dd1,T dd3,T M);
};
//###########################################################################

/*!Numerically calibrated parameters from arXiv:1508.07253
 * see the table in the data directory for labeled version
 */
const double lambda_num_params[19][11]=
{{ 3931.9 ,  -17395.8 ,  3132.38 ,  343966. ,  -1.21626e6 ,  -70698. ,
    1.38391e6 ,  -3.96628e6 ,  -60017.5 ,  803515. , 
   -2.09171e6 }, 
{ -40105.5 ,  112253. ,  23561.7 ,  -3.47618e6 , 
   1.1375900000000002e7 ,  754313. ,  -1.30848e7 ,  3.64446e7 , 
   596227. ,  -7.42779e6 ,  1.8929e7 }, 
{ 83208.4 ,  -191238. , 
   -210916. ,  8.71798e6 ,  -2.69149e7 ,  -1.98898e6 ,  3.0888e7 , 
   -8.39087e7 ,  -1.4535e6 ,  1.70635e7 ,  -4.27487e7 }, 
{ 0.814984 , 
   2.57476 ,  1.16102 ,  -2.36278 ,  6.77104 ,  0.757078 ,  -2.72569 ,
    7.11404 ,  0.176693 ,  -0.797869 ,  2.11624 }, 
{ 0.0069274 , 
   0.0302047 ,  0.00630802 ,  -0.120741 ,  0.262716 ,  0.00341518 , 
   -0.107793 ,  0.27099 ,  0.000737419 ,  -0.0274962 , 
   0.0733151 }, 
{ 1.01034 ,  0.000899312 ,  0.283949 ,  -4.04975 , 
   13.2078 ,  0.103963 ,  -7.02506 ,  24.7849 ,  0.030932 ,  -2.6924 ,
    9.60937 }, 
{ 1.30816 ,  -0.00553773 ,  -0.0678292 ,  -0.668983 , 
   3.40315 ,  -0.0529658 ,  -0.992379 ,  4.82068 ,  -0.00613414 , 
   -0.384293 ,  1.75618 }, 
{ 2096.55 ,  1463.75 ,  1312.55 , 
   18307.3 ,  -43534.1 ,  -833.289 ,  32047.3 ,  -108609. ,  452.251 ,
    8353.44 ,  -44531.3 }, 
{ -10114.1 ,  -44631. ,  -6541.31 , 
   -266959. ,  686328. ,  3405.64 ,  -437508. ,  1.63182e6 , 
   -7462.65 ,  -114585. ,  674402. }, 
{ 22933.7 ,  230960. , 
   14961.1 ,  1.19402e6 ,  -3.10422e6 ,  -3038.17 ,  1.87203e6 , 
   -7.30915e6 ,  42738.2 ,  467502. ,  -3.06485e6 }, 
{ -14621.7 , 
   -377813. ,  -9608.68 ,  -1.71089e6 ,  4.33292e6 ,  -22366.7 , 
   -2.50197e6 ,  1.02745e7 ,  -85360.3 ,  -570025. , 
   4.39684e6 }, 
{ 97.8975 ,  -42.6597 ,  153.484 ,  -1417.06 , 
   2752.86 ,  138.741 ,  -1433.66 ,  2857.74 ,  41.0251 ,  -423.681 , 
   850.359 }, 
{ -3.2827 ,  -9.05138 ,  -12.4154 ,  55.4716 , 
   -106.051 ,  -11.953 ,  76.807 ,  -155.332 ,  -3.41293 ,  25.5724 , 
   -54.408 }, 
{ -0.000025156400000000002 ,  0.000019750300000000003 , 
   -0.0000183707 ,  0.0000218863 ,  0.0000825024 ,  7.15737e-6 , 
   -0.000055780000000000005 ,  0.000191421 ,  5.44717e-6 , 
   -0.0000322061 ,  0.0000797402 }, 
{ 43.3151 ,  638.633 ,  -32.8577 ,
    2415.89 ,  -5766.88 ,  -61.8546 ,  2953.97 ,  -8986.29 , 
   -21.5714 ,  981.216 ,  -3239.57 }, 
{ -0.0702021 ,  -0.162698 , 
   -0.187251 ,  1.13831 ,  -2.83342 ,  -0.17138 ,  1.71975 , 
   -4.53972 ,  -0.0499834 ,  0.606207 ,  -1.68277 }, 
{ 9.59881 , 
   -397.054 ,  16.2021 ,  -1574.83 ,  3600.34 ,  27.0924 ,  -1786.48 ,
    5152.92 ,  11.1757 ,  -577.8 ,  1808.73 }, 
{ -0.0298949 , 
   1.40221 ,  -0.0735605 ,  0.833701 ,  0.224001 ,  -0.0552029 , 
   0.566719 ,  0.718693 ,  -0.0155074 ,  0.157503 , 
   0.210768 }, 
{ 0.997441 ,  -0.00788445 ,  -0.0590469 ,  1.39587 , 
   -4.51663 ,  -0.0558534 ,  1.75166 ,  -5.99021 ,  -0.0179453 , 
   0.59651 ,  -2.06089 }};
#endif
