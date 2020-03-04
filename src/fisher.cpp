#include <fisher.h>
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <adolc/adolc_sparse.h>
#include <math.h>
#include <string>
#include "util.h"
#include "detector_util.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "waveform_generator.h"
#include "waveform_util.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <fstream>


using namespace std;



/*! \file 
 *
 */
/*
 * All subroutines associated with waveform differentiation and Fisher analysis
 *
 * Fisher options (both autodiff and numerical):
 *
 * IMRPhenomD sky-averaged (7) -- ln A0, phic, tc, ln chirpmass, ln eta, chi_symm, chi_antisymm
 * 
 * ppE_IMRPhenomD_Inspiral sky-averaged (7 + mods)-- ln A0, phic, tc, ln chirpmass, ln eta, chi_symm, chi_antisymm, betas
 * 
 * ppE_IMRPhenomD_IMR sky-averaged (7 + mods)-- ln A0, phic, tc, ln chirpmass, ln eta, chi_symm, chi_antisymm, betas
 *
 * IMRPhenomD !sky_averaged (11) --  RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, eta, chi1, chi2
 * 
 * ppE_IMRPhenomD_Inspiral !sky_averaged (11+mods) -- RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, eta, chi1, chi2, betas
 * 
 * ppE_IMRPhenomD_IMR !sky_averaged (11+mods) -- RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, eta, chi1, chi2, betas
 * 
 * dCS_IMRPhenomD !sky_averaged (12) -- RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, eta, chi1, chi2, \alpha^2 (sec^4)
 * 
 * EdGB_IMRPhenomD !sky_averaged (12) --RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, eta, chi1, chi2, \alpha^2 (sec^4)
 * 
 * IMRPhenomPv2 !sky_averaged (13) --  RA, DEC, psi,phiRef, tc, \iota_L,  ln DL, ln chirpmass, eta, chi_1 \dot \hat{L}, chi_2 \dot \hat{L} ,chi_p,  \phi_p (all at f_ref)
 * 
 * ppE_IMRPhenomPv2_Inspiral !sky_averaged (13 + mods) --  RA, DEC, psi,phiRef, tc, \iota_L, ln DL, ln chirpmass, ln eta, chi_1 \dot \hat{L}, chi_2 \dot \hat{L} ,chi_p,  \phi_p, \vec{beta} (all at f_ref)
 * 
 * ppE_IMRPhenomPv2_IMR !sky_averaged (13 + mods) --  RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, ln eta, chi_1 \dot \hat{L}, chi_2 \dot \hat{L} ,chi_p,  \phi_p, \vec{beta} (all at f_ref)
 * 
 * dCS_IMRPhenomPv2 !sky_averaged (14) --  RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, ln eta, chi_1 \dot \hat{L}, chi_2 \dot \hat{L} ,chi_p,  \phi_p,  \alpha (all at f_ref)
 * 
 * EdGB_IMRPhenomPv2 !sky_averaged (14) --  RA, DEC, psi, phiRef, tc, \iota_L, ln DL, ln chirpmass, ln eta, chi_1 \dot \hat{L}, chi_2 \dot \hat{L} ,chi_p,  \phi_p,  \alpha (all at f_ref)
 *
 *
 * All MCMC options correspond to the base, minus the coalescence time (which is maximized over) -- Reduced MCMC option correspond to the options above minus t_c (and phic for sky_averaged) -- Non reduced correspond to replacing \chi_1 \dot \hat{L}, \chi_2 \dot \hat{L}, \chi_p and \phi_p with |\chi_1|, |\chi_2|, \theta_1, \theta_2, \phi_1, and \phi_2
 *
 * Both fisher functions are compatible with OpenMP, but not with pthreads (if ADOL-C was installed with the openmp option)
 *
 * NOTE about non-sky-averaged: Matrices will most likely be singular if only using one detector, and still largly singular with 2. 3 detectors will mostly guarantee that the matrices are invertible. 
 */

/*!\brief Calculates the fisher matrix for the given arguments
 *
 * Utilizes numerical derivatives -- non-skyaveraged supports up to 4th order finite difference (sky averaged supports second order only)
 */
void fisher_numerical(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	string generation_method, 
	string detector, 
	double **output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params_base<double> *parameters,
	int order,/**< Order of the numerical derivative (2 or 4)**/
	//double *parameters,
	int *amp_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method -- if using numerical derivatives or speed isn't that important, just set to NULL*/
	int *phase_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	double *noise
	)
{
	//populate noise and frequency
	double *internal_noise = new double[length];
	if (noise)
	{
		for(int i = 0 ; i < length;i++)
		{
			internal_noise[i] = noise[i];
		}
	}
	else{
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
	}
		
	//populate derivatives - Derivatives of DETECTOR RESPONSE
	std::complex<double> **response_deriv = new std::complex<double>*[dimension];
	for (int i = 0 ; i<dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}
	
	calculate_derivatives(response_deriv, 
			frequency,
			length, 
			dimension, 
			detector, 
			generation_method,
			parameters,
			order);


	//calulate fisher elements
	bool log10_f=false;
	std::string integration_method="SIMPSONS";
	double *weights=NULL;
	calculate_fisher_elements(frequency, length,dimension, response_deriv, output,  internal_noise,integration_method,weights,log10_f);
	//Factor of 2 for LISA's second arm
	if(detector == "LISA"){
		for(int i = 0 ; i<dimension;i++){
			for(int j = 0  ;j<dimension; j++){
				output[i][j]*=2;
			}	
		}
	}
	for (int i = 0 ; i<dimension; i++){
		delete [] response_deriv[i];	
	}
	delete [] response_deriv;
	delete [] internal_noise;
}


void calculate_derivatives(std::complex<double>  **response_deriv, 
       	double *frequencies,
       	int length, 
       	int dimension, 
       	string detector, 
       	string  gen_method,
       	gen_params_base<double> *parameters,
	int order)
{
	double epsilon = 1e-6;
	//Order of numerical derivative
	double parameters_vec[dimension];
	bool log_factors[dimension];
	double param_p[dimension];
	double param_m[dimension];
	double *param_pp;
	double *param_mm;
	if(order >= 4){
		param_pp = new double[dimension];
		param_mm = new double[dimension];
	}
	//##########################################################
	std::string local_gen_method = local_generation_method(gen_method);
	unpack_parameters(parameters_vec, parameters, gen_method, dimension, log_factors);
	//##########################################################
	gen_params waveform_params;
	repack_non_parameter_options(&waveform_params,parameters, gen_method);
	//##########################################################
	if(parameters->sky_average)
	{
		double *amplitude_plus = new double[length];
		double *phase_plus = new double[length];
		double *phase_plusc = new double[length];
		double *amplitude_minus = new double[length];
		double *phase_minus = new double[length];
		double *phase_minusc= new double[length];
		double *amplitude = new double[length];
		double *amplitude_plus_plus; 
		double *amplitude_minus_minus;
		double *phase_plus_plus;
		double *phase_minus_minus;
		double *phase_plus_plusc;
		double *phase_minus_minusc;
		if(order>=4){
			amplitude_plus_plus = new double[length];
			amplitude_minus_minus = new double[length];
			phase_plus_plus = new double[length];
			phase_minus_minus = new double[length];
			phase_plus_plusc = new double[length];
			phase_minus_minusc = new double[length];
		}
		fourier_amplitude(frequencies, 
			length,
			amplitude,
			local_gen_method,
			parameters);	
		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = parameters_vec[j] ;
				param_m[j] = parameters_vec[j] ;
			}
			param_p[i] = parameters_vec[i] + epsilon;
			param_m[i] = parameters_vec[i] - epsilon;
			if(order>=4){
				for( int j =0;j<dimension;j++){
					param_pp[j] = parameters_vec[j] ;
					param_mm[j] = parameters_vec[j] ;
				}
				param_pp[i] = parameters_vec[i] + 2*epsilon;
				param_mm[i] = parameters_vec[i] - 2*epsilon;

			}
			repack_parameters(param_p, &waveform_params, gen_method, dimension, parameters);
			fourier_amplitude(frequencies, 
				length,
				amplitude_plus,
				local_gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus,
				phase_plusc,
				local_gen_method,
				&waveform_params);	

			repack_parameters(param_m, &waveform_params, gen_method, dimension, parameters);
			fourier_amplitude(frequencies, 
				length,
				amplitude_minus,
				local_gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_minus,
				phase_minusc,
				local_gen_method,
				&waveform_params);	
			if(order>=4){
				repack_parameters(param_pp, &waveform_params, gen_method, dimension, parameters);
				fourier_amplitude(frequencies, 
					length,
					amplitude_plus_plus,
					local_gen_method,
					&waveform_params);	
				fourier_phase(frequencies, 
					length,
					phase_plus_plus,
					phase_plus_plusc,
					local_gen_method,
					&waveform_params);	

				repack_parameters(param_mm, &waveform_params, gen_method, dimension, parameters);
				fourier_amplitude(frequencies, 
					length,
					amplitude_minus_minus,
					local_gen_method,
					&waveform_params);	
				fourier_phase(frequencies, 
					length,
					phase_minus_minus,
					phase_minus_minusc,
					local_gen_method,
					&waveform_params);	
			}
			double amplitude_deriv, phase_deriv;
			if(order==2){
				for (int l =0;l<length;l++)
				{
					amplitude_deriv = (amplitude_plus[l] -amplitude_minus[l])/(2*epsilon);
					phase_deriv = (phase_plus[l] -phase_minus[l])/(2*epsilon);
					response_deriv[i][l] = amplitude_deriv + 
						std::complex<double>(0,1)*phase_deriv*amplitude[l];
				}
			}
			else if(order==4){
				for (int l =0;l<length;l++)
				{
					amplitude_deriv = (-amplitude_plus_plus[l]+8.*amplitude_plus[l] -8.*amplitude_minus[l]+amplitude_minus_minus[l])/(12.*epsilon);
					phase_deriv = (-phase_plus_plus[l]+8.*phase_plus[l] -8.*phase_minus[l]+phase_minus_minus[l])/(12.*epsilon);
					response_deriv[i][l] = amplitude_deriv - 
						std::complex<double>(0,1)*phase_deriv*amplitude[l];
				}
			}

				
		}
		delete [] amplitude_plus;
		delete [] amplitude_minus; 
		delete [] phase_plus; 
		delete [] phase_minus;
		delete [] phase_plusc; 
		delete [] phase_minusc;
		delete [] amplitude;
		if(order>=4){
			delete [] amplitude_plus_plus;
			delete [] amplitude_minus_minus;
			delete [] phase_plus_plus;
			delete [] phase_minus_minus;
			delete [] phase_plus_plusc;
			delete [] phase_minus_minusc;
		}

	
	}
	//##########################################################
	else {
		std::complex<double> *response_plus= new std::complex<double>[length];
		std::complex<double> *response_minus= new std::complex<double>[length];
		std::complex<double> *response_plus_plus;
		std::complex<double> *response_minus_minus;
		double *times=NULL;
		double **dt=NULL;
		bool corr_time;
		int local_dimension=dimension;
		if(detector=="LISA"){
			times = new double[length];
			time_phase_corrected(times, length, frequencies,parameters, gen_method, false);
			dt = allocate_2D_array(dimension, length);
			corr_time = false;
			time_phase_corrected_derivative_numerical(dt, length, frequencies,parameters, gen_method, dimension, corr_time);
			//local_dimension++;
		}
		if(order >=4){
			response_plus_plus= new std::complex<double>[length];
			response_minus_minus= new std::complex<double>[length];
		}
		for (int i =0; i<local_dimension; i++){
			for( int j =0;j<local_dimension;j++){
				param_p[j] = parameters_vec[j] ;
				param_m[j] = parameters_vec[j] ;
			}
			param_p[i] = parameters_vec[i] + epsilon;
			param_m[i] = parameters_vec[i] - epsilon;
			if(order>=4){
				for( int j =0;j<dimension;j++){
					param_pp[j] = parameters_vec[j] ;
					param_mm[j] = parameters_vec[j] ;
				}
				param_pp[i] = parameters_vec[i] +2 *epsilon;
				param_mm[i] = parameters_vec[i] -2 *epsilon;
			}
			repack_parameters(param_p, &waveform_params, gen_method, dimension, parameters);
			if(detector=="LISA"){
				//correct time needs to stay false for now
				//time_phase_corrected(times, length,frequencies,  &waveform_params, local_gen_method, corr_time);
				//map_extrinsic_angles(&waveform_params);
			}
			fourier_detector_response(frequencies, 
				length,
				response_plus,
				detector,
				local_gen_method,
				&waveform_params,
				times);	

			repack_parameters(param_m, &waveform_params, gen_method, dimension, parameters);
			if(detector=="LISA"){
				//map_extrinsic_angles(&waveform_params);
			}
			fourier_detector_response(frequencies, 
				length,
				response_minus,
				detector,
				local_gen_method,
				&waveform_params,
				times);	
			if(order>=4){
				repack_parameters(param_pp, &waveform_params, gen_method, dimension, parameters);
				if(detector=="LISA"){
					//map_extrinsic_angles(&waveform_params);
				}
				fourier_detector_response(frequencies, 
					length,
					response_plus_plus,
					detector,
					local_gen_method,
					&waveform_params,
					times);	

				repack_parameters(param_mm, &waveform_params, gen_method, dimension, parameters);
				if(detector=="LISA"){
					//map_extrinsic_angles(&waveform_params);
				}
				fourier_detector_response(frequencies, 
					length,
					response_minus_minus,
					detector,
					local_gen_method,
					&waveform_params,
					times);	
			}
			if(order == 2){
				for (int l =0;l<length;l++)
				{
					response_deriv[i][l] = 
						(response_plus[l]-response_minus[l])/(2.*epsilon);
				}
			}
			else if(order == 4){
				for (int l =0;l<length;l++)
				{
					response_deriv[i][l] = 
						(-response_plus_plus[l]+8.*response_plus[l]-8.*response_minus[l]+response_minus_minus[l])/(12.*epsilon);
				}
			}
				
				
		}
		if (detector=="LISA"){
		//if (false){
			//Calculate derivative wrt time of waveform
			repack_parameters(param_p, 
				&waveform_params, 
				gen_method, 
				dimension, 
				parameters);
			fourier_detector_response(frequencies, 
				length,
				response_plus,
				detector,
				local_gen_method,
				&waveform_params,
				times);	
			std::complex<double> *deriv_t = new std::complex<double>[length];
			//deriv_t[0] = (response_plus[1] - response_plus[0])/(times[1]-times[0]);
			//One sided difference for now, because the difference 
			//in time is not constant 
			//(central difference depends on symmetric spacing?
			for(int i = 0 ; i<length-1; i++){
				deriv_t[i] = 
					(response_plus[i+1]-response_plus[i])/(times[i+1]-times[i]);
			}
			deriv_t[length-1] = 
				(response_plus[length-1]-response_plus[length-2])/(times[length-1]-times[length-2]);
			for(int i = 0 ; i<dimension; i++){
				for(int j = 0 ; j<length; j++){
					response_deriv[i][j]+=deriv_t[j] * dt[i][j];
				}
			}
			//do chain rule with dt above
			
			delete [] deriv_t;
		}

		delete [] response_plus;
		delete [] response_minus;
		if(detector=="LISA"){
			delete [] times;
			deallocate_2D_array(dt, dimension, length);
		}
		if(order>=4){
			delete [] response_plus_plus, response_minus_minus;
		}
	}
	if(order>= 4){
		delete [] param_pp;
		delete [] param_mm;
	}
	for(int l =0 ; l<dimension; l++){
		if(log_factors[l]){
			for(int j = 0 ; j<length; j++){
				response_deriv[l][j] *=parameters_vec[l] ;
			}
		}
	}
	deallocate_non_param_options(&waveform_params, parameters, gen_method);

}

/*!\brief Calculates the fisher matrix for the given arguments to within numerical error using automatic differention in ``batch'' mode for modifications to GR
 *
 * Built  around  ADOL-C -- A. Walther und A. Griewank: Getting started with ADOL-C. In U. Naumann und O. Schenk, Combinatorial Scientific Computing, Chapman-Hall CRC Computational Science, pp. 181-202 (2012).
 *
 * This constructs the fisher for a list of modifications in the usual way, but skips elements of the fisher that correspond to (mod, mod) elements. Specifically, it calculates all the fisher elements except (i, j ) i>GR_dimension && j > GR_dimension && i!=j.
 *
 * To find the fisher for one of the modifications, simply remove all the other  dimensions associated with the extra modifications using rm_fisher_dim in util.h
 *
 * !NOTE!:This routine only works as intended when GR is the injected value, that is all the  betas are evaluated at 0. And since the covariances between modifications are not computed, this should only be used to look at one modification at a time.
 */
void fisher_autodiff_batch_mod(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string detector, 
	double **output,/**< double [dimension][dimension]*/
	int base_dimension, /**< GR dimensionality*/
	int full_dimension, /**< Total dimension of the output fisher (ie GR_dimension + Nmod)*/
	gen_params *parameters,
	//double *parameters,
	int *amp_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	int *phase_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	double *noise
	)
{
	//populate noise and frequency
	double *internal_noise;
	bool local_noise=false;
	if (noise)
	{
		internal_noise = noise;
	}
	else{
		internal_noise = new double[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
		local_noise=true;
	}
		
	//populate derivatives

	std::complex<double> **response_deriv = new std::complex<double>*[full_dimension];
	for(int i =0 ;i<full_dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}
	bool autodiff_time_deriv=false;
	calculate_derivatives_autodiff(frequency,length, full_dimension,generation_method, parameters, response_deriv, NULL, detector,autodiff_time_deriv);
	//##########################################################
	
	//calulate fisher elements
	calculate_fisher_elements_batch(frequency, length,base_dimension, full_dimension, response_deriv, output,  internal_noise);
	//Factor of 2 for LISA's second arm
	if(detector == "LISA"){
		for(int i = 0 ; i<full_dimension;i++){
			for(int j = 0  ;j<full_dimension; j++){
				output[i][j]*=2;
			}	
		}
	}

	if(local_noise){delete [] internal_noise;}
	for(int i =0 ;i<full_dimension; i++){
		delete [] response_deriv[i];
	}
	delete [] response_deriv;
}
/*!\brief Calculates the fisher matrix for the given arguments to within numerical error using automatic differention - slower than the numerical version
 *
 * Build  around  ADOL-C -- A. Walther und A. Griewank: Getting started with ADOL-C. In U. Naumann und O. Schenk, Combinatorial Scientific Computing, Chapman-Hall CRC Computational Science, pp. 181-202 (2012).
 */
void fisher_autodiff_interp(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string detector, 
	double **output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params *parameters,
	int downsampling_factor,
	//double *parameters,
	int *amp_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	int *phase_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	double *noise
	)
{
	//populate noise and frequency
	double *internal_noise;
	bool local_noise=false;
	if (noise)
	{
		internal_noise = noise;
	}
	else{
		internal_noise = new double[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
		local_noise=true;
	}
		
	//populate derivatives

	std::complex<double> **response_deriv = new std::complex<double>*[dimension];
	for(int i =0 ;i<dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}

	//Decimate frequencies by downsampling factor
	//Ensure the final point is the end of the array
	int res = length%downsampling_factor;
	int length_ds = length/downsampling_factor; //downsampled length
	if(res != 0) length_ds+=1;
	std::complex<double> **temp_deriv = new std::complex<double>*[dimension];
	double *temp_deriv_c = new double[length_ds];
	double *temp_deriv_phase = new double[length_ds];
	double *freqs_ds = new double[length_ds];
	for(int i = 0 ; i<dimension; i++){
		temp_deriv[i] = new std::complex<double>[length_ds];
	}
	int ct = 0 ;
	for(int i = 0 ; i<length; i++){
		if(i%downsampling_factor ==0  && ct < length_ds-1){
			freqs_ds[ct] = frequency[i];
			ct++;
		}
	}
	freqs_ds[length_ds-1] = frequency[length-1];
	//calculate_derivatives_autodiff(frequency,length, dimension,generation_method, parameters, response_deriv, NULL, detector);
	bool autodiff_time_deriv=false;
	calculate_derivatives_autodiff(freqs_ds,length_ds, dimension,generation_method, parameters, temp_deriv, NULL, detector,autodiff_time_deriv);
	//Interpolate derivatives here to get back to full length
	gsl_interp_accel *my_accel_ptr= gsl_interp_accel_alloc ();
	gsl_spline *my_spline_ptr= gsl_spline_alloc (gsl_interp_linear, length_ds);
	double val;
	for(int i =0 ; i<dimension; i++){
		for(int j =0 ;  j<length_ds; j++){
			temp_deriv_c[j] = std::abs(temp_deriv[i][j]);
		}
		gsl_spline_init (my_spline_ptr, freqs_ds, temp_deriv_c, length_ds);
		for(int j=0; j<length; j++){
			response_deriv[i][j] = gsl_spline_eval (my_spline_ptr, frequency[j], my_accel_ptr);
			//y_deriv2 = gsl_spline_eval_deriv2 (my_spline_ptr, x, my_accel_ptr);
		}
		for(int j =0 ;  j<length_ds; j++){
			temp_deriv_c[j] = std::arg(temp_deriv[i][j]);
		}
		unwrap_array(temp_deriv_c,temp_deriv_phase,length_ds);
		gsl_spline_init (my_spline_ptr, freqs_ds, temp_deriv_phase, length_ds);
		for(int j=0; j<length; j++){
			//response_deriv[i][j] +=std::complex<double>(0, gsl_spline_eval (my_spline_ptr, frequency[j], my_accel_ptr));
			response_deriv[i][j] *=std::exp(std::complex<double>(0, gsl_spline_eval (my_spline_ptr, frequency[j], my_accel_ptr)));
			//y_deriv2 = gsl_spline_eval_deriv2 (my_spline_ptr, x, my_accel_ptr);
		}
	}
	gsl_spline_free (my_spline_ptr);
	gsl_interp_accel_free(my_accel_ptr);
	//##########################################################
	
	//calulate fisher elements
	std::string integration_method = "SIMPSONS";
	bool log10_f = false;
	double *weights = NULL;
	calculate_fisher_elements(frequency, length,dimension, response_deriv, output,  internal_noise,   integration_method,weights,log10_f  );

	//Factor of 2 for LISA's second arm
	if(detector == "LISA"){
		for(int i = 0 ; i<dimension;i++){
			for(int j = 0  ;j<dimension; j++){
				output[i][j]*=2;
			}	
		}
	}

	if(local_noise){delete [] internal_noise;}
	for(int i =0 ;i<dimension; i++){
		delete [] response_deriv[i];
		delete [] temp_deriv[i];
	}
	delete [] response_deriv;
	delete [] temp_deriv;
	delete [] temp_deriv_c;
	delete [] temp_deriv_phase;
	delete [] freqs_ds;
}

/*! \brief Calculates the fisher matrix for the given arguments to within numerical error using automatic differention - slower than the numerical version
 *
 * Build  around  ADOL-C -- A. Walther und A. Griewank: Getting started with ADOL-C. In U. Naumann und O. Schenk, Combinatorial Scientific Computing, Chapman-Hall CRC Computational Science, pp. 181-202 (2012).
 *
 * Note: If trying to run large number of Fishers concurrently with openmp, you may need to change some ADOL-C parameters. Copy the example .adolcrc file in data/ to the working directory and adjust as needed
 *
 * If integration method is gaussian quadrature and the frequencies are !NULL, the weights must be pre-allocated. 
 *
 * If the integration method is Newton-Cotes, the weights are not used.
 *
 * Integration schemes:
 *
 * Newton-Cotes methods:
 *
 * 	SIMPSONS
 * 	
 * Gaussian Quadrature
 *
 * 	GAUSSLEG
 *
 * log10_f used if frequencies are logarithmically (base 10) spaced. NOTE: not sure if this would work for SIMPSONS -- this hasn't been tested. For now, only gaussian quadrature routines are supported.
 */
void fisher_autodiff(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string detector, 
	double **output,/**< double [dimension][dimension]*/
	int dimension,/**<dimension of the fisher*/ 
	gen_params *parameters,/**< Injection parameters*/
	std::string integration_method,/**< Method of integration to use*/
	double * weights,/**< If using a gaussian quadrature method and the weights have been precomputed, the weights can be supplied here*/
	bool log10_f,/**< Boolean for logarithmically (base 10) spaced frequencies*/
	double *noise,/**<Precomputed PSD array*/
	//double *parameters,
	int *amp_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	int *phase_tapes/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	)
{
	//populate noise and frequency
	double *internal_noise;
	bool local_noise=false;
	if (noise)
	{
		internal_noise = noise;
	}
	else{
		internal_noise = new double[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
		local_noise=true;
	}
		
	//populate derivatives

	std::complex<double> **response_deriv = new std::complex<double>*[dimension];
	for(int i =0 ;i<dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}
	bool autodiff_time_deriv=false;
	//Gauss quad is unequal spacing, doesn't work with numerical derivatives
	if(integration_method == "GAUSSLEG"){
		autodiff_time_deriv = true;
	}
	calculate_derivatives_autodiff(frequency,length, dimension,generation_method, parameters, response_deriv, NULL, detector,autodiff_time_deriv);
	//##########################################################
	
	//calulate fisher elements
	calculate_fisher_elements(frequency, length,dimension, response_deriv, output,  internal_noise, integration_method,weights,log10_f  );

	//Factor of 2 for LISA's second arm
	if(detector == "LISA"){
		for(int i = 0 ; i<dimension;i++){
			for(int j = 0  ;j<dimension; j++){
				output[i][j]*=2;
			}	
		}
	}

	if(local_noise){delete [] internal_noise;}
	for(int i =0 ;i<dimension; i++){
		delete [] response_deriv[i];
	}
	delete [] response_deriv;
}

/*! \brief Calculates the derivatives of the detector response using automatic differentiation
 *
 * Possibly slower than the numerical derivative, but not susceptible to truncation error from finite difference
 *
 * Higher dimensional fishers actually could be faster
 *
 * NOTE: dimension parameter ALWAYS refers to the dimension of the fisher (ie the length of the source parameter vector), even though the derivatives are computed wrt dimension +1 or dimension + 2 -- the +1(+2) are for the frequency deriv(time deriv)
 *
 * autodiff_time_deriv indicates if the time derivatives for space detectors should be done, either with the full autodiff hessian (slow) or with the numerical approximation (autodiff then numerical)
 */
void calculate_derivatives_autodiff(double *frequency,
	int length,
	int dimension,
	std::string generation_method,
	gen_params *parameters,
	std::complex<double> **waveform_deriv,
	int *waveform_tapes,
	std::string detector,
	bool autodiff_time_deriv
	)
{
	//Transform gen_params to double vectors
	//double vec_parameters[dimension+1];
	int vec_param_length= dimension +1;
	if(detector == "LISA"){
		//take derivative wrt time as well, for the chain rule
		vec_param_length += 1;
	}
	double vec_parameters[vec_param_length];
	bool log_factors[dimension];
	int boundary_num= boundary_number(generation_method);
	if(boundary_num == -1){
		std::cout<<"Error -- unsupported generation method"<<std::endl;
		exit(1);
	}
	double *freq_boundaries=new double[boundary_num];
	double *grad_freqs=new double[boundary_num];
	std::string local_gen_method = local_generation_method(generation_method);
	//prep_fisher_calculation(vec_parameters,log_factors, freq_boundaries,grad_freqs,boundary_num,parameters, generation_method, dimension);
	assign_freq_boundaries(freq_boundaries, grad_freqs,boundary_num, parameters, generation_method);
	vec_parameters[0]=grad_freqs[0];
	unpack_parameters(&vec_parameters[1], parameters, generation_method,dimension, log_factors);
	double *grad_times=NULL;
	double **dt=NULL;
	double *eval_times=NULL;
	if(detector == "LISA"){
		grad_times = new double[boundary_num];
		time_phase_corrected_autodiff(grad_times, boundary_num, grad_freqs, parameters, generation_method, false);
		//for(int  i  = 0 ;  i<boundary_num;  i++){
		//	grad_times[i] =  1;
		//}
		dt = allocate_2D_array(dimension+1, length);	
		if(autodiff_time_deriv){
			time_phase_corrected_derivative_autodiff_full_hess(dt, length, frequency, parameters, generation_method, dimension, false);
		}
		else{
			time_phase_corrected_derivative_autodiff_numerical(dt, length, frequency, parameters, generation_method, dimension, false);
		}
		//time_phase_corrected_derivative_autodiff(dt, length, frequency, parameters, generation_method, dimension, false);
		eval_times = new double[length];
		time_phase_corrected_autodiff(eval_times, length, frequency, parameters, generation_method, false);
		//for(int  i  = 0 ;  i<length;  i++){
		//	eval_times[i] =  1;
		//}
			
	}
	//calculate_derivative tapes
	int tapes[boundary_num];
	for(int i =0; i<boundary_num; i++){
		tapes[i]=i;
		trace_on(tapes[i]);
		//adouble avec_parameters[dimension+1];
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,parameters,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension, parameters);
		adouble time;
		if(detector == "LISA"){
			time <<= grad_times[i];
			//map_extrinsic_angles(&a_parameters);
		}
		std::complex<adouble> a_response;
		if(!a_parameters.sky_average){
			//int status  = fourier_detector_response_equatorial(&afreq, 1, &a_response, detector, local_gen_method, &a_parameters, times);
			int status  = fourier_detector_response(&afreq, 1, &a_response, detector, local_gen_method, &a_parameters, &time);

		}
		else{
			adouble a_amp;
			adouble a_phasep;
			adouble a_phasec;
		
			int status  = fourier_amplitude(&afreq, 1, &a_amp, local_gen_method, &a_parameters);
			status  = fourier_phase(&afreq, 1, &a_phasep,  &a_phasec,local_gen_method, &a_parameters);
			a_response = a_amp * exp(std::complex<adouble>(0,a_phasep));

		}
		double response[2];
		real(a_response) >>=  response[0];	
		imag(a_response) >>=  response[1];	

		trace_off();
		//if(times){
		//	delete  [] times;
		//}
		deallocate_non_param_options(&a_parameters, parameters, generation_method);
	}
	//Evaluate derivative tapes
	int dep = 2;//Output is complex
	int indep = vec_param_length;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	double **jacob = allocate_2D_array(dep,indep);
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequency[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				if(detector == "LISA"){
					vec_parameters[vec_param_length -1] = eval_times[k];
				}
				jacobian(tapes[n], dep, indep, vec_parameters, jacob);
				for(int i =0; i<dimension; i++){
					waveform_deriv[i][k] = jacob[0][i+1] 
						+ std::complex<double>(0,1)*jacob[1][i+1];
					//correct for time deriv for LISA
					if(detector == "LISA"){
					//if(false){
						//std::complex<double> temp=waveform_deriv[i][k];
						waveform_deriv[i][k]+= 
							(jacob[0][vec_param_length-1] + std::complex<double>(0,1)*jacob[1][vec_param_length-1]) //Time derivative of WF
							* dt[i+1][k];//Derivative of time wrt source parameter
					}
				}
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		//If freq didn't fall in any boundary, set to 0
		if(!eval){
			for(int i =0; i<dimension; i++){
				waveform_deriv[i][k] = std::complex<double>(0,0);
			}	
		}
		eval = false;
	}
	//Account for Log parameters
	for(int j = 0 ; j<dimension; j++){
		if(log_factors[j]){
			for (int i = 0;i <length; i++)
			{
				//j+1 for vec_parameter because of freq in position 0
				waveform_deriv[j][i] *= (vec_parameters[j+1]);
			}
		}
	}
	deallocate_2D_array(jacob,dep,indep);
	if(freq_boundaries){
		delete [] freq_boundaries;
	}
	if(grad_freqs){
		delete [] grad_freqs;
	}
	if(grad_times){
		delete [] grad_times;	
	}
	if(dt){
		deallocate_2D_array(dt, dimension+1, length);
	}
	if(eval_times){
		delete [] eval_times;
	}

}

void num_src_params(int *N_src_params, std::string generation_method, gen_params_base<double> *params)
{
	if(generation_method.find("IMRPhenomPv2")!=std::string::npos){
		*N_src_params = 9+1;	
	}
	else if(generation_method.find("IMRPhenomD")!=std::string::npos){
		*N_src_params = 6+1;	
	}
	if(check_mod(generation_method))
	{
		*N_src_params += params->Nmod;
	}
}
void reduce_extrinsic(int *src_params, int N_src_params, std::string generation_method, gen_params_base<double>*params)
{
	int gr_dim, gr_param_dim;
	if(generation_method.find("IMRPhenomPv2")!=std::string::npos){
		src_params[0]=0;
		src_params[1]=4;
		src_params[2]=5;
		src_params[3]=6;
		src_params[4]=8;
		src_params[5]=9;
		src_params[6]=10;
		src_params[7]=11;
		src_params[8]=12;
		src_params[9]=13;
		gr_dim = 14;
		gr_param_dim = 10;
	}
	else if(generation_method.find("IMRPhenomD")!=std::string::npos){
		src_params[0]=0;
		src_params[1]=4;
		src_params[2]=5;
		src_params[3]=8;
		src_params[4]=9;
		src_params[5]=10;
		src_params[6]=11;
		gr_dim = 12;
		gr_param_dim = 7;
	}
	if(check_mod(generation_method)){
		for(int i = 0; i<params->Nmod;i++){
			src_params[gr_dim+i]=gr_param_dim+i;
		}
	}
}
/*! \brief Computes the derivative of the phase w.r.t. source parameters AS DEFINED BY FISHER FILE -- hessian of the phase
 *
 * If specific derivatives need to taken, take this routine as a template and write it yourself.
 *
 * The dt array has shape [dimension+1][length] (dimension + 1 for the frequency derivative, so dimension should only include the source parameters)
 *
 * This takes the autodiff derivative wrt source parameters, and a numerical derivative for the frequency
 *
 */
void time_phase_corrected_derivative_autodiff_numerical(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time)
{
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	int vec_param_length = dimension +1 ;//+1 for frequency 
	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	std::string local_gen_method = local_generation_method(generation_method);
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);
	double vec_parameters[vec_param_length];
	bool log_factors[dimension];
	unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
	
	for(int i = 0 ; i<vec_param_length; i++){
		for(int j = 0 ; j<length; j++){
			dt[i][j] = 0;
		}
	}
	//calculate derivative of phase
	int tapes[boundary_num];
	for(int i = 0 ; i < boundary_num ; i++){
		tapes[i] = i*12; //Random tape id 
		trace_on(tapes[i]);
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,params,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension,params);
		if(params->equatorial_orientation){
			transform_orientation_coords(&a_parameters, generation_method, "");
		}
		adouble time;
		adouble phasep, phasec;
		int status  = fourier_phase(&afreq, 1, &phasep,&phasec, local_gen_method, &a_parameters);
		double phase;
		phasep >>= phase;

		trace_off();
		deallocate_non_param_options(&a_parameters, params, generation_method);
	}
	int indep = vec_param_length;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	int dep = 1;
	double **jacob = allocate_2D_array(dep,indep);
	double **source_param_deriv= allocate_2D_array(indep, length);
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequencies[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				jacobian(tapes[n], dep, indep, vec_parameters, jacob);
				for(int i =0; i<vec_param_length; i++){
					source_param_deriv[i][k] = jacob[0][i];
				}
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		//If freq didn't fall in any boundary, set to 0
		if(!eval){
			for(int i =0; i<vec_param_length; i++){
				source_param_deriv[i][k] = 0.;
			}	
		}
		eval = false;
	}
	deallocate_2D_array(jacob, dep, indep);
	double deltaf = frequencies[1]-frequencies[0];
	for(int k = 0 ; k<vec_param_length; k++){
		//forward deriv
		dt[k][0] = (source_param_deriv[k][1] - source_param_deriv[k][0])/(2*M_PI*deltaf);
		//central difference 2nd
		dt[k][1] = (source_param_deriv[k][2] - source_param_deriv[k][0])/(4*M_PI*deltaf);
		//Central difference 4th order
		for(int i = 2 ; i<length-2; i ++){
			//dt[k][i]=(source_param_deriv[k][i+1] - source_param_deriv[k][i-1])/(4*M_PI*deltaf);
			dt[k][i]=(-source_param_deriv[k][i+2] + 8.*source_param_deriv[k][i+1]-8.*source_param_deriv[k][i-1]+source_param_deriv[k][i-2])/(24.*M_PI*deltaf);
		}
		//central difference 2nd
		dt[k][length-2] = (source_param_deriv[k][length-1] - source_param_deriv[k][length-3])/(4.*M_PI*deltaf);
		//backwards deriv
		dt[k][length - 1] = (source_param_deriv[k][length-1] - source_param_deriv[k][length-2])/(2.*M_PI*deltaf);
	}
	//divide by 2 PI
	//for(int j = 0 ; j < vec_param_length; j++){
	//	for(int i = 0 ; i<length; i++){
	//		dt[j][i]/=(2.*M_PI);
	//	}
	//}
	deallocate_2D_array(source_param_deriv, indep, length);

}
/*! \brief Computes the derivative of the phase w.r.t. source parameters AS DEFINED BY FISHER FILE -- hessian of the phase
 *
 * If specific derivatives need to taken, take this routine as a template and write it yourself.
 *
 * The dt array has shape [dimension+1][length] (dimension + 1 for the frequency derivative, so dimension should only include the source parameters)
 *
 */
void time_phase_corrected_derivative_autodiff(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time)
{
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	int vec_param_length = dimension +1 ;//+1 for frequency 
	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	std::string local_gen_method = local_generation_method(generation_method);
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);
	double vec_parameters[vec_param_length];
	bool log_factors[dimension];
	unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
	
	//int source_param[7] = {5,8,9,10,11,12,13};
	//int num_sp = 7;
	int N_src_params;
	num_src_params(&N_src_params, generation_method, params);
	int src_params[N_src_params];
	reduce_extrinsic(src_params, N_src_params, generation_method, params);
	
	for(int i = 0 ; i<vec_param_length; i++){
		for(int j = 0 ; j<length; j++){
			dt[i][j] = 0;
		}
	}
	//calculate derivative of phase
	int tapes[boundary_num];
	for(int i = 0 ; i < boundary_num ; i++){
		tapes[i] = i*12; //Random tape id 
		trace_on(tapes[i]);
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,params,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension,params);
		if(params->equatorial_orientation){
			transform_orientation_coords(&a_parameters, generation_method, "");
		}
		adouble time;
		adouble phasep, phasec;
		int status  = fourier_phase(&afreq, 1, &phasep,&phasec, local_gen_method, &a_parameters);
		double phase;
		phasep >>= phase;

		trace_off();
		deallocate_non_param_options(&a_parameters, params, generation_method);
	}
	int indep = vec_param_length;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	int dep = 1;
	int p = 2;
	int d = 2;
	double **S = new double*[indep];
        for (int k=0; k<indep; k++) {
            S[k] = new double[p];
            for (int j=0; j<p; j++)
                S[k][j] = (k==j)?1.0:0.0;
        }
	int dim = binomi(p+d,d);
	int jvec[d];
	jvec[1]=1;//Freq derivative
        double **tensor = new double*[dep];
        for(int k = 0 ; k<dep; k++)
             tensor[k] = new double[dim];
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequencies[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				//tensor_eval(tapes[n], dep,indep, d,p,vec_parameters,tensor,S );
				for(int i =0; i<indep; i++){
					if(check_list(i, src_params, N_src_params)){
						for(int l = 0 ; l<indep; l++){
							S[l][1]=(i==l)? 1.0:0.0;
						}
						tensor_eval(tapes[n], dep,indep, d,p,vec_parameters,tensor,S );
						//jvec[0]=i+1;	
						jvec[0]=2;	
						dt[i][k] = tensor[0][tensor_address(d,jvec)] ;
					}
				}
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		//If freq didn't fall in any boundary, set to 0
		if(!eval){
			for(int i =0; i<vec_param_length; i++){
				dt[i][k] = 0.;
			}	
		}
		eval = false;
	}
	//deallocate_2D_array(hess,indep,indep);
	//divide by 2 PI
	for(int j = 0 ; j < vec_param_length; j++){
		for(int i = 0 ; i<length; i++){
			dt[j][i]/=(2.*M_PI);
		}
	}
	for(int k = 0; k<indep; k++){
		delete [] S[k];
	}
	for(int k = 0; k<dep; k++){
		delete [] tensor[k];
	}
	delete [] S;
	delete [] tensor;

}
/*! \brief Computes the derivative of the phase w.r.t. source parameters AS DEFINED BY FISHER FILE -- hessian of the phase
 *
 * If specific derivatives need to taken, take this routine as a template and write it yourself.
 *
 * The dt array has shape [dimension+1][length] (dimension + 1 for the frequency derivative, so dimension should only include the source parameters)
 *
 */
void time_phase_corrected_derivative_autodiff_sparse(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time)
{
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	int vec_param_length = dimension +1 ;//+1 for frequency 
	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	std::string local_gen_method = local_generation_method(generation_method);
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);
	double vec_parameters[vec_param_length];
	bool log_factors[dimension];
	unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
	
	for(int i = 0 ; i<vec_param_length; i++){
		for(int j = 0 ; j<length; j++){
			dt[i][j] = 0;
		}
	}
	//calculate derivative of phase
	int tapes[boundary_num];
	for(int i = 0 ; i < boundary_num ; i++){
		tapes[i] = i*12; //Random tape id 
		trace_on(tapes[i]);
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,params,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension,params);
		if(params->equatorial_orientation){
			transform_orientation_coords(&a_parameters, generation_method, "");
		}
		adouble time;
		adouble phasep, phasec;
		int status  = fourier_phase(&afreq, 1, &phasep,&phasec, local_gen_method, &a_parameters);
		double phase;
		phasep >>= phase;

		trace_off();
		deallocate_non_param_options(&a_parameters, params, generation_method);
	}
	int indep = vec_param_length;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	//double **hess = allocate_2D_array(indep,indep);
	unsigned int *rind=NULL;
	unsigned int *cind=NULL;
	double *values=NULL;
	int options[2]; options[0]=0; options[1]=0;
	int nnz;
	vec_parameters[0]=frequencies[0];
	sparse_hess(tapes[0], indep, 0,vec_parameters,&nnz, &rind, &cind, &values,options );
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequencies[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				sparse_hess(tapes[n], indep, 1,vec_parameters,&nnz, &rind, &cind, &values,options );
				for(int i =0; i<nnz; i++){
					if(rind[i]==0){
						dt[cind[i]][k] = values[i] ;
					}
					//for(int j = 0 ; j<nnz; j++){
					//	std::cout<<rind[j]<<" "<<cind[j]<<" "<<values[j]<<std::endl;
					//}
					//std::cout<<std::endl;
				}
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		//If freq didn't fall in any boundary, set to 0
		if(!eval){
			for(int i =0; i<vec_param_length; i++){
				dt[i][k] = 0.;
			}	
		}
		eval = false;
	}
	//deallocate_2D_array(hess,indep,indep);
	//divide by 2 PI
	for(int j = 0 ; j < vec_param_length; j++){
		for(int i = 0 ; i<length; i++){
			dt[j][i]/=(2.*M_PI);
		}
	}
	free(rind);free(cind);free(values);

}
/*! \brief Computes the derivative of the phase w.r.t. source parameters AS DEFINED BY FISHER FILE -- hessian of the phase
 *
 * If specific derivatives need to taken, take this routine as a template and write it yourself.
 *
 * The dt array has shape [dimension+1][length] (dimension + 1 for the frequency derivative, so dimension should only include the source parameters)
 *
 */
void time_phase_corrected_derivative_autodiff_full_hess(double **dt, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, int dimension, bool correct_time)
{
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	int vec_param_length = dimension +1 ;//+1 for frequency 
	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	std::string local_gen_method = local_generation_method(generation_method);
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);
	double vec_parameters[vec_param_length];
	bool log_factors[dimension];
	unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
	
	//calculate derivative of phase
	int tapes[boundary_num];
	for(int i = 0 ; i < boundary_num ; i++){
		tapes[i] = i*12; //Random tape id 
		trace_on(tapes[i]);
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,params,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension,params);
		if(params->equatorial_orientation){
			transform_orientation_coords(&a_parameters, generation_method, "");
		}
		adouble time;
		adouble phasep, phasec;
		int status  = fourier_phase(&afreq, 1, &phasep,&phasec, local_gen_method, &a_parameters);
		double phase;
		phasep >>= phase;

		trace_off();
		deallocate_non_param_options(&a_parameters, params, generation_method);
	}
	int indep = vec_param_length;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	double **hess = allocate_2D_array(indep,indep);
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequencies[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				hessian(tapes[n], indep, vec_parameters, hess);
				for(int i =0; i<vec_param_length; i++){
					dt[i][k] = hess[i][0] ;
					//std::cout<<hess[i][0]<<" ";
				}
				//std::cout<<std::endl;
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		//If freq didn't fall in any boundary, set to 0
		if(!eval){
			for(int i =0; i<vec_param_length; i++){
				dt[i][k] = 0.;
			}	
		}
		eval = false;
	}
	deallocate_2D_array(hess,indep,indep);
	//divide by 2 PI
	for(int j = 0 ; j < vec_param_length; j++){
		for(int i = 0 ; i<length; i++){
			dt[j][i]/=(2.*M_PI);
		}
	}

}
/*! \brief Computes the derivative of the phase w.r.t. source parameters AS DEFINED BY FISHER FILE -- hessian of the phase -- numerical 
 *
 * IN PROGRESS -- DO NOT USE
 *
 * If specific derivatives need to taken, take this routine as a template and write it yourself.
 *
 * The dt array has shape [dimension+1][length] (dimension + 1 for the frequency derivative, so dimension should only include the (full) source parameters)
 *
 */
template<class T>
void time_phase_corrected_derivative_numerical(T **dt, int length, T *frequencies,gen_params_base<T> *params, std::string generation_method, int dimension, bool correct_time)
{
	//bool save_shift_time = params->shift_time;
	//params->shift_time = false;
	//std::string local_gen = "IMRPhenomD";
	std::string local_gen=generation_method;
	//################################################
	T *phase_pp = new T[length];//source param plus
	T *phase_cross = new T[length];
	T *phase_pm = new T[length];//source param minus
	bool log_factors[dimension];
	T *parameters_vec = new T[dimension];
	T *param_p = new T[dimension];
	T *param_m = new T[dimension];
	double epsilon = 1e-6;
	T deltaf = frequencies[1]-frequencies[0];
	lambda_parameters<T> lambda;
	source_parameters<T> s_param;
	//################################################
	//
	//################################################
	//std::string local_gen_method = local_generation_method(gen_method);
	unpack_parameters(parameters_vec, params, generation_method, dimension, log_factors);
	//##########################################################
	for (int k =0; k<dimension; k++){
		for( int j =0;j<dimension;j++){
			param_p[j] = parameters_vec[j] ;
			param_m[j] = parameters_vec[j] ;
		}
		param_p[k] = parameters_vec[k] + epsilon;
		param_m[k] = parameters_vec[k] - epsilon;
		gen_params waveform_params;
		repack_non_parameter_options(&waveform_params,params, generation_method);
		repack_parameters(param_p, &waveform_params, generation_method, dimension, params);
		fourier_phase(frequencies, length, phase_pp, phase_cross, local_gen, &waveform_params);

		repack_parameters(param_m, &waveform_params, generation_method, dimension, params);
		fourier_phase(frequencies, length, phase_pm, phase_cross, local_gen, &waveform_params);
		//################################################
		T fRD, fdamp,fpeak;
		if(local_gen.find("IMRPhenomPv2")!=std::string::npos){
			IMRPhenomPv2<T> modelp;
			s_param = source_parameters<T>::populate_source_parameters(params);
			s_param.spin1z = params->spin1[2];
			s_param.spin2z = params->spin2[2];
			s_param.chip = params->chip;
			s_param.phip = params->phip;
			s_param.phiRef = params->phiRef;
			s_param.f_ref = params->f_ref;
			s_param.incl_angle = params->incl_angle;
			modelp.PhenomPv2_Param_Transform_reduced(&s_param);
			s_param.sky_average = params->sky_average;
			s_param.cosmology=params->cosmology;
			modelp.assign_lambda_param(&s_param,&lambda);	
			modelp.post_merger_variables(&s_param);
			fRD = s_param.fRD;
			fdamp = s_param.fdamp;
			fpeak = modelp.fpeak(&s_param , &lambda);
		}
		else if(local_gen.find("IMRPhenomD")!=std::string::npos){
			IMRPhenomD<T> model;
			s_param = source_parameters<T>::populate_source_parameters(params);
			s_param.sky_average = params->sky_average;
			s_param.f_ref = params->f_ref;
			s_param.phiRef = params->phiRef;
			s_param.cosmology=params->cosmology;
			s_param.incl_angle=params->incl_angle;
			model.assign_lambda_param(&s_param,&lambda);	
			model.post_merger_variables(&s_param);
			fRD = s_param.fRD;
			fdamp = s_param.fdamp;
			fpeak = model.fpeak(&s_param , &lambda);
		}
		//################################################
		//Factor of 2 pi for the definition of time from frequency
		//if(local_gen == "IMRPhenomD"){
		if(local_gen.find("IMRPhenom")!=std::string::npos){
			//Currently using Nico's fix
			if(correct_time){
				//T f = frequencies[0];
				//bool check = true, check2=true;
				//T pt1, pt2, f1,f2;
				//int i = 0 ;
				////One sided, to start it off
				//dt[k][0] = (phase_plus[1]-phase_plus[0])/(2.*M_PI*deltaf);
				//i++;
				//while(f < .95*fRD && i<length-1)
				//{
				//	f = frequencies[i];
				//	//central difference for the rest of the steps
				//	dt[i] = (phase_plus[i+1]-phase_plus[i-1])/(4.*M_PI*deltaf);
				//	if(check){
				//		if(f>fpeak){
				//			pt1 = dt[k][i];
				//			f1 = f;
				//			check=false;
				//		}
				//	}
				//	else{
				//		if(check2){
				//			if(f>.9*fRD){
				//				pt2 = dt[k][i];
				//				f2 = f;
				//				check2=false;
				//			}
				//		}
				//	}
				//	i++;
				//}	
				//T f_intercept = dt[k][i-1];
				//T f_mr = f;
				//T freq_slope_pm = (pt2-pt1)/(f2-f1);
				//while(f<1.5*fRD && i<length-1)
				//{
				//	f = frequencies[i];
				//	dt[k][i] =f_intercept+ (f-f_mr)*freq_slope_pm;
				//	//Stop if observation goes past 20 years
				//	i++;
				//	if(dt[k][i-1]>(params->tc+630720000)){
				//		break;
				//	}

				//}
				//T time_transition = dt[k][i-1];
				//T f_transition = f;
				//while(i<length-1)
				//{
				//	f = frequencies[i];
				//	dt[k][i] = time_transition +pow_int(f-f_transition,2);
				//	//Stop if observation goes past 20 years
				//	i++;
				//	if(dt[k][i-1]>(params->tc+630720000)){
				//		break;
				//	}
				//	
				//}
				////if stopped at 20 years, fill out with frozen time
				//if(i != length){
				//	while(i<length){
				//		dt[k][i] = dt[k][i-1];
				//		i++;
				//	}
				//}
				//else{
				//	dt[k][length-1] = (phase_plus[length-1] - phase_plus[length-2])/(2*M_PI*deltaf);
				//}
			}
			else{
				//dt[k][0] = (phase_pp[1]-phase_pm[1]-phase_pp[0] + phase_pm[0])/(8*M_PI*deltaf * epsilon);
				if(length>2){
					for(int i = 1  ;i<length-1; i++){
						dt[k][i] = (phase_pp[i+1] - phase_pm[i+1] -phase_pp[i-1]+phase_pm[i-1])/(8*M_PI*deltaf*epsilon);
					}
				//dt[k][length-1] = (phase_pp[length-1]-phase_pm[length-1]-phase_pp[length-2]+phase_pm[length-2])/(4*M_PI*deltaf*epsilon);
					//YES THIS IS WRONG
					dt[k][length-1] = dt[k][length-2];
					dt[k][0] = dt[k][1];
				
				}
			}
		}
	}
	//################################################
	delete [] phase_pp;
	delete [] phase_pm;
	delete [] phase_cross;
	delete [] param_p;
	delete [] param_m;
	delete [] parameters_vec;
	//params->shift_time = save_shift_time;
}
template void time_phase_corrected_derivative_numerical<double>(double **, int, double *, gen_params_base<double> *, std::string, int, bool);
/*! \brief Utility for mapping generation method string to one accepted by the waveform_generation routines
 *
 * Certain combinations of parameters are labeled by generation method strings not under the waveform_generation routines, so a transformation is needed
 */
std::string local_generation_method(std::string generation_method)
{
	std::string local_gen_method = generation_method;
	if(generation_method.find("MCMC") != std::string::npos && generation_method.find("Full") != std::string::npos)
	{
		local_gen_method.erase(0,5);
		local_gen_method.erase(local_gen_method.length()-5,5);
	}
	else if(generation_method.find("MCMC") != std::string::npos)
	{
		local_gen_method.erase(0,5);
	}
	return local_gen_method;
		
}
/*! \brief Adjust parameters for detector specific configurations (namely, LISA introduces extra transitions that needs to be accounted for)
 *
 * This is kept separate to improve the modularity of the code. Waveform specific parameters are taken care of in prep_fisher_calculation and detector specific parameters are taken care of here. 
 */
void detect_adjust_parameters( double *freq_boundaries,double *grad_freqs, int *boundary_num,gen_params_base<double> *input_params, std::string generation_method, std::string detector,int dim)
{
	if(detector == "LISA"){
		if(generation_method.find("IMRPhenom") != std::string::npos){
			gen_params_base<adouble> internal_params;
			transform_parameters(input_params, &internal_params);
			source_parameters<adouble> s_param;
			s_param = source_parameters<adouble>::populate_source_parameters(&internal_params);
			s_param.sky_average = internal_params.sky_average;
			s_param.f_ref = internal_params.f_ref;
			s_param.phiRef = internal_params.phiRef;
			s_param.shift_time = false;
			s_param.cosmology=internal_params.cosmology;
			s_param.incl_angle=internal_params.incl_angle;
			lambda_parameters<adouble> lambda;
			double M, fRD, fpeak;
			if(generation_method.find("IMRPhenomPv2") != std::string::npos){
				IMRPhenomPv2<adouble> modelp;
				s_param.spin1z = internal_params.spin1[2];
				s_param.spin2z = internal_params.spin2[2];
				s_param.chip = internal_params.chip;
				s_param.phip = internal_params.phip;
				modelp.PhenomPv2_Param_Transform_reduced(&s_param);
				modelp.assign_lambda_param(&s_param, &lambda);
				modelp.post_merger_variables(&s_param);
				M = s_param.M.value();
				fRD = s_param.fRD.value();
				fpeak = modelp.fpeak(&s_param, &lambda).value();
			}
			else if(generation_method.find("IMRPhenomD")!=std::string::npos){
				IMRPhenomD<adouble> modeld;
				modeld.assign_lambda_param(&s_param, &lambda);
				modeld.post_merger_variables(&s_param);
				M = s_param.M.value();
				fRD = s_param.fRD.value();
				fpeak = modeld.fpeak(&s_param, &lambda).value();
			}
		}
	}

}
/*! \brief Unpacks the input gen_params object into a double array for use with the fisher routines
 */
void unpack_parameters(double *parameters, gen_params_base<double> *input_params, std::string generation_method, int dimension, bool *log_factors)
{
	if(!input_params->sky_average){
		if(generation_method.find("IMRPhenomPv2") != std::string::npos){
			if(generation_method.find("MCMC") != std::string::npos){
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				parameters[0]=input_params->RA;
				parameters[1]=sin(input_params->DEC);
				if(input_params->equatorial_orientation){
					parameters[2]=input_params->theta_l;
					parameters[3]=input_params->phi_l;
				}
				else{
					parameters[2]=input_params->psi;
					parameters[3]=cos(input_params->incl_angle);
				}
				parameters[4]=input_params->phiRef;
				parameters[5]=input_params->tc;
				parameters[6]=log(input_params->Luminosity_Distance);
				parameters[7]=log(calculate_chirpmass(input_params->mass1, 
					input_params->mass2));
				parameters[8]=calculate_eta(input_params->mass1, 
					input_params->mass2);
				//parameters[9]=input_params->spin1[2];
				//parameters[10]=input_params->spin2[2];
				//parameters[11]=input_params->chip;
				//parameters[12]=input_params->phip;
				double spin1sph[3];
				double spin2sph[3];
				transform_cart_sph(input_params->spin1,spin1sph);
				transform_cart_sph(input_params->spin2,spin2sph);
				parameters[9]=spin1sph[0];
				parameters[10]=spin2sph[0];
				parameters[11]=cos(spin1sph[1]);
				parameters[12]=cos(spin2sph[1]);
				parameters[13]=spin1sph[2]-spin2sph[2];

			}
			else{
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				log_factors[6] = true;//Distance
				log_factors[7] = true;//chirpmass

				parameters[0]=input_params->RA;
				parameters[1]=input_params->DEC;
				if(input_params->equatorial_orientation){
					parameters[2]=input_params->theta_l;
					parameters[3]=input_params->phi_l;
				}
				else{
					parameters[2]=input_params->psi;
					parameters[3]=input_params->incl_angle;
				}
				parameters[4]=input_params->phiRef;
				parameters[5]=input_params->tc;
				parameters[6]=input_params->Luminosity_Distance;
				parameters[7]=calculate_chirpmass(input_params->mass1, input_params->mass2);
				parameters[8]=calculate_eta(input_params->mass1, input_params->mass2);
				parameters[9]=input_params->spin1[2];
				parameters[10]=input_params->spin2[2];
				parameters[11]=input_params->chip;
				parameters[12]=input_params->phip;

			}

		}
		else if(generation_method.find("IMRPhenomD") != std::string::npos){
			if(generation_method.find("MCMC") != std::string::npos){
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				//double spin1spher[3];
				//double spin2spher[3];
				//parameters[0]=input_params->RA;
				//parameters[1]=input_params->DEC;
				//if(input_params->equatorial_orientation){
				//	parameters[2]=input_params->theta_l;
				//	parameters[3]=input_params->phi_l;
				//}
				//else{
				//	parameters[2]=input_params->psi;
				//	parameters[3]=cos(input_params->incl_angle);
				//}
				//parameters[4]=log(input_params->Luminosity_Distance);
				//parameters[5]=log(calculate_chirpmass(input_params->mass1, 
				//	input_params->mass2));
				//parameters[6]=calculate_eta(input_params->mass1, 
				//	input_params->mass2);
				//parameters[7]=input_params->spin1[2];
				//parameters[8]=input_params->spin2[2];
				
				double spin1spher[3];
				double spin2spher[3];
				parameters[0]=input_params->RA;
				parameters[1]=sin(input_params->DEC);
				if(input_params->equatorial_orientation){
					parameters[2]=input_params->theta_l;
					parameters[3]=input_params->phi_l;
				}
				else{
					parameters[2]=input_params->psi;
					parameters[3]=cos(input_params->incl_angle);
				}
				parameters[4] = input_params->phiRef;
				parameters[5] = input_params->tc;
				parameters[6]=log(input_params->Luminosity_Distance);
				parameters[7]=log(calculate_chirpmass(input_params->mass1, 
					input_params->mass2));
				parameters[8]=calculate_eta(input_params->mass1, 
					input_params->mass2);
				parameters[9]=input_params->spin1[2];
				parameters[10]=input_params->spin2[2];

			}
			else{
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				log_factors[6] = true;//Distance
				log_factors[7] = true;//chirpmass

				parameters[0]=input_params->RA;
				parameters[1]=input_params->DEC;
				if(input_params->equatorial_orientation){
					parameters[2]=input_params->theta_l;
					parameters[3]=input_params->phi_l;
				}
				else{
					parameters[2]=input_params->psi;
					parameters[3]=input_params->incl_angle;
				}
				parameters[4]=input_params->phic;
				parameters[5]=input_params->tc;
				parameters[6]=input_params->Luminosity_Distance;
				parameters[7]=calculate_chirpmass(input_params->mass1, 
					input_params->mass2);
				parameters[8]=calculate_eta(input_params->mass1, 
					input_params->mass2);
				parameters[9]=input_params->spin1[2];
				parameters[10]=input_params->spin2[2];

			}
		
		}

	}
	else{
		if(generation_method.find("IMRPhenomPv2") != std::string::npos){
			//Need to populate
			if(generation_method.find("MCMC") != std::string::npos){
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				parameters[0]=log(calculate_chirpmass(input_params->mass1, 
					input_params->mass2));
				parameters[1]=calculate_eta(input_params->mass1, 
					input_params->mass2);
				double spin1sph[3];
				double spin2sph[3];
				transform_cart_sph(input_params->spin1,spin1sph);
				transform_cart_sph(input_params->spin2,spin2sph);
				parameters[2]=spin1sph[0];
				parameters[3]=spin2sph[0];
				parameters[4]=cos(spin1sph[1]);
				parameters[5]=cos(spin2sph[1]);
				parameters[6]=spin1sph[2]-spin2sph[2];


			}
			else{

			}
	
		}
		else if(generation_method.find("IMRPhenomD") != std::string::npos){
			if(generation_method.find("MCMC") != std::string::npos){
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				//log_factors[0] = true;//chirpmass

				parameters[0]=log(calculate_chirpmass(input_params->mass1, 
					input_params->mass2));
				parameters[1]=calculate_eta(input_params->mass1, 
					input_params->mass2);
				parameters[2]=input_params->spin1[2];
				parameters[3]=input_params->spin2[2];

			}
			else{
				for(int i = 0 ; i<dimension; i++){
					log_factors[i] = false;
				}
				log_factors[0] = true;//A0
				log_factors[3] = true;//chirpmass
				log_factors[4] = true;//eta

				parameters[3] = calculate_chirpmass(input_params->mass1, input_params->mass2);
				parameters[0] = A0_from_DL(parameters[3]*MSOL_SEC,
					input_params->Luminosity_Distance*MPC_SEC,
					input_params->sky_average);
				parameters[1] = input_params->phic;
				parameters[2] = input_params->tc;
				parameters[4] = calculate_eta(input_params->mass1, 
					input_params->mass2);
				parameters[5] = (input_params->spin1[2] + input_params->spin2[2])/2.;
				parameters[6] = (input_params->spin1[2] - input_params->spin2[2])/2.;

			}
		
		}
	}
	if( check_mod(generation_method)){
		int base = dimension-input_params->Nmod;
		for(int i = 0 ;i<input_params->Nmod; i++){
			parameters[base+i] = input_params->betappe[i];
		}
	}

}
/*! \brief Repack the parameters from an adouble vector to a gen_params_base<adouble> object and freqeuncy 
 *
 * This is one of the places where the generation-method/dimension/sky_average specific modifications should go
 */
template<class T>
void repack_parameters(T *avec_parameters, gen_params_base<T> *a_params, std::string generation_method, int dim, gen_params_base<double> *original_params)
{
	if(!a_params->sky_average){
		if(generation_method.find("IMRPhenomPv2") != std::string::npos){
			if(generation_method.find("MCMC")!=std::string::npos){
				a_params->mass1 = calculate_mass1(exp(avec_parameters[7]),
					avec_parameters[8]);
				a_params->mass2 = calculate_mass2(exp(avec_parameters[7]),
					avec_parameters[8]);
				a_params->Luminosity_Distance = exp(avec_parameters[6]);
				a_params->RA = avec_parameters[0];
				a_params->DEC = asin(avec_parameters[1]);
				if(a_params->equatorial_orientation){
					a_params->theta_l = avec_parameters[2];
					a_params->phi_l=avec_parameters[3];
				}
				else{
					a_params->psi = avec_parameters[2];
					a_params->incl_angle=acos(avec_parameters[3]);
				}
				a_params->phiRef = avec_parameters[4];
				//a_params->chi1_l = avec_parameters[8];
				//a_params->chi2_l = avec_parameters[9];
				//a_params->spin1[2] = avec_parameters[9];
				//a_params->spin2[2] = avec_parameters[10];
				//a_params->chip = avec_parameters[11];
				//a_params->phip = avec_parameters[12];

				//###################################
				//Fishers don't necessarily respect the bounds of 
				//acos, so bump it back into range if we're on the edge
				//##################################
				T local_theta1;
				if(avec_parameters[11] > 1){
					local_theta1 = 0;
				}
				else if(avec_parameters[11]<-1){
					local_theta1 = M_PI;
				}
				else{
					local_theta1=acos(avec_parameters[11]);
				}
				T local_theta2;
				if(avec_parameters[12] > 1){
					local_theta2 = 0;
				}
				else if(avec_parameters[12]<-1){
					local_theta2 = M_PI;
				}
				else{
					local_theta2=acos(avec_parameters[12]);
				}
				a_params->tc=avec_parameters[5];
				//##################################
				//##################################
				T spin1cart[3];
				T spin1sph[3] = {avec_parameters[9],local_theta1,
					avec_parameters[13]};
				T spin2cart[3];
				T spin2sph[3] = {avec_parameters[10],local_theta2,0};
				transform_sph_cart(spin1sph,spin1cart);
				transform_sph_cart(spin2sph,spin2cart);
				a_params->spin1[0] = spin1cart[0];
				a_params->spin1[1] = spin1cart[1];
				a_params->spin1[2] = spin1cart[2];
				a_params->spin2[0] = spin2cart[0];
				a_params->spin2[1] = spin2cart[1];
				a_params->spin2[2] = spin2cart[2];
				


			}
			else{
				a_params->mass1 = calculate_mass1(avec_parameters[7],
					avec_parameters[8]);
				a_params->mass2 = calculate_mass2(avec_parameters[7],
					avec_parameters[8]);
				a_params->Luminosity_Distance = avec_parameters[6];
				a_params->RA = avec_parameters[0];
				a_params->DEC = avec_parameters[1];
				if(a_params->equatorial_orientation){
					a_params->theta_l = avec_parameters[2];
					a_params->phi_l=avec_parameters[3];
				}
				else{
					a_params->psi = avec_parameters[2];
					a_params->incl_angle=avec_parameters[3];
				}
				a_params->phiRef = avec_parameters[4];
				a_params->tc = avec_parameters[5];
				//a_params->chi1_l = avec_parameters[9];
				//a_params->chi2_l = avec_parameters[10];
				a_params->spin1[2] = avec_parameters[9];
				a_params->spin2[2] = avec_parameters[10];
				a_params->chip = avec_parameters[11];
				a_params->phip = avec_parameters[12];

			}	
		}	
		else if(generation_method.find("IMRPhenomD") != std::string::npos){
			if(generation_method.find("MCMC")!=std::string::npos){
				a_params->mass1 = calculate_mass1(exp(avec_parameters[7]),
					avec_parameters[8]);
				a_params->mass2 = calculate_mass2(exp(avec_parameters[7]),
					avec_parameters[8]);
				a_params->Luminosity_Distance = exp(avec_parameters[6]);
				a_params->RA = avec_parameters[0];
				a_params->DEC = asin(avec_parameters[1]);
				if(a_params->equatorial_orientation){
					a_params->theta_l = avec_parameters[2];
					a_params->phi_l=avec_parameters[3];
				}
				else{
					a_params->psi = avec_parameters[2];
					a_params->incl_angle=acos(avec_parameters[3]);
				}
				//T spin1sph[3] = {avec_parameters[9],0,0};
				//T spin2sph[3] = {avec_parameters[10],0,0};
				a_params->spin1[2] = avec_parameters[9];
				a_params->spin2[2] = avec_parameters[10];
				//transform_sph_cart(spin1sph,a_params->spin1);
				//transform_sph_cart(spin2sph,a_params->spin2);
				//maximized out
				a_params->phiRef = avec_parameters[4];
				a_params->phic = 0;
				a_params->tc = avec_parameters[5];

				//a_params->mass1 = calculate_mass1(exp(avec_parameters[5]),
				//	avec_parameters[6]);
				//a_params->mass2 = calculate_mass2(exp(avec_parameters[5]),
				//	avec_parameters[6]);
				//a_params->Luminosity_Distance = exp(avec_parameters[4]);
				//a_params->RA = avec_parameters[0];
				//a_params->DEC = avec_parameters[1];
				//if(a_params->equatorial_orientation){
				//	a_params->theta_l = avec_parameters[2];
				//	a_params->phi_l=avec_parameters[3];
				//}
				//else{
				//	a_params->psi = avec_parameters[2];
				//	a_params->incl_angle=(avec_parameters[3]);
				//}
				//a_params->spin1[2] = avec_parameters[7];
				//a_params->spin2[2] = avec_parameters[8];
				////maximized out
				//a_params->phiRef = 1;
				//a_params->phic = 0;
				//a_params->tc = 1;

			}
			else{
				a_params->mass1 = calculate_mass1(avec_parameters[7],
					avec_parameters[8]);
				a_params->mass2 = calculate_mass2(avec_parameters[7],
					avec_parameters[8]);
				a_params->Luminosity_Distance = avec_parameters[6];
				a_params->RA = avec_parameters[0];
				a_params->DEC = avec_parameters[1];
				if(a_params->equatorial_orientation){
					a_params->theta_l = avec_parameters[2];
					a_params->phi_l=avec_parameters[3];
				}
				else{
					a_params->psi = avec_parameters[2];
					a_params->incl_angle=avec_parameters[3];
				}
				a_params->phic = avec_parameters[4];
				a_params->tc = avec_parameters[5];
				//T spin1sph[3] = {avec_parameters[9],0,0};
				//T spin2sph[3] = {avec_parameters[10],0,0};
				//transform_sph_cart(spin1sph,a_params->spin1);
				//transform_sph_cart(spin2sph,a_params->spin2);
				a_params->spin1[2] = avec_parameters[9];
				a_params->spin2[2] = avec_parameters[10];

			}	
		}	
	}
	else{
		if(generation_method.find("IMRPhenomPv2") != std::string::npos){
			if(generation_method.find("MCMC")!=std::string::npos){

				a_params->mass1 = calculate_mass1(exp(avec_parameters[0]),avec_parameters[1]);
				a_params->mass2 = calculate_mass2(exp(avec_parameters[0]),avec_parameters[1]);

				T local_theta1;
				if(avec_parameters[4] > 1){
					local_theta1 = 0;
				}
				else if(avec_parameters[4]<-1){
					local_theta1 = M_PI;
				}
				else{
					local_theta1=acos(avec_parameters[4]);
				}
				T local_theta2;
				if(avec_parameters[5] > 1){
					local_theta2 = 0;
				}
				else if(avec_parameters[5]<-1){
					local_theta2 = M_PI;
				}
				else{
					local_theta2=acos(avec_parameters[5]);
				}
				//##################################
				//##################################
				T spin1cart[3];
				T spin1sph[3] = {avec_parameters[2],local_theta1,
					avec_parameters[6]};
				T spin2cart[3];
				T spin2sph[3] = {avec_parameters[3],local_theta2,0};
				transform_sph_cart(spin1sph,spin1cart);
				transform_sph_cart(spin2sph,spin2cart);
				a_params->spin1[0] = spin1cart[0];
				a_params->spin1[1] = spin1cart[1];
				a_params->spin1[2] = spin1cart[2];
				a_params->spin2[0] = spin2cart[0];
				a_params->spin2[1] = spin2cart[1];
				a_params->spin2[2] = spin2cart[2];

				a_params->tc=0;
				a_params->phiRef=0;
				a_params->RA=0;
				a_params->DEC=0;
				a_params->psi=0;
				a_params->Luminosity_Distance=100;

			}
			else{

			}	
		}	
		else if(generation_method.find("IMRPhenomD") != std::string::npos){
			if(generation_method.find("MCMC")!=std::string::npos){
				a_params->mass1 = calculate_mass1(exp(avec_parameters[0]),
					avec_parameters[1]);
				a_params->mass2 = calculate_mass2(exp(avec_parameters[0]),
					avec_parameters[1]);
				a_params->Luminosity_Distance = 1000;
				//T spin1sph[3] = {avec_parameters[2],0,0};
				//T spin2sph[3] = {avec_parameters[3],0,0};
				//transform_sph_cart(spin1sph,a_params->spin1);
				//transform_sph_cart(spin2sph,a_params->spin2);
				a_params->spin1[2] = avec_parameters[2];
				a_params->spin2[2] = avec_parameters[3];
				a_params->phic=0;
				a_params->phiRef=0;
				a_params->tc=0;
				//a_params->incl_angle=0;
				//a_params->RA=1.;
				//a_params->DEC=1.;
				//a_params->psi=1.;
				//a_params->gmst=1185389807.3;

			}
			else{
				a_params->mass1 = calculate_mass1(avec_parameters[3],
					avec_parameters[4]);
				a_params->mass2 = calculate_mass2(avec_parameters[3],
					avec_parameters[4]);
				a_params->Luminosity_Distance = 
					DL_from_A0((T)(avec_parameters[3]*MSOL_SEC),
					avec_parameters[0],a_params->sky_average)/MPC_SEC;
				a_params->tc = avec_parameters[2];
				a_params->phic = avec_parameters[1];
				T chi1 = avec_parameters[5] + avec_parameters[6];
				T chi2 = avec_parameters[5] - avec_parameters[6];
				//T spin1sph[3] = {chi1,0,0};
				//T spin2sph[3] = {chi2,0,0};
				//transform_sph_cart(spin1sph,a_params->spin1);
				//transform_sph_cart(spin2sph,a_params->spin2);
				a_params->spin1[2]  = chi1;
				a_params->spin2[2]  = chi2;
				//a_params->incl_angle=0;
				//a_params->RA=1.;
				//a_params->DEC=1.;
				//a_params->psi=1.;
				//a_params->gmst=1185389807.3;
			}	
		}	

	}
	if( check_mod(generation_method)){
		int base = dim - a_params->Nmod;
		for(int i = 0 ;i<a_params->Nmod; i++){
			a_params->betappe[i] = avec_parameters[base+i];
		}
	}

}

/*! \brief Utilitiy to transfer non-parameter options from one gen_params structure to another
 *
 * If ppE waveform ALLOCATES MEMORY -- MUST be deallocated
 */
template<class T>
void repack_non_parameter_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method)
{
	waveform_params->sky_average = input_params->sky_average;
	waveform_params->f_ref = input_params->f_ref;
	waveform_params->gmst = input_params->gmst;
	waveform_params->horizon_coord = input_params->horizon_coord;
	waveform_params->equatorial_orientation = input_params->equatorial_orientation;
	waveform_params->NSflag1 = input_params->NSflag1;
	waveform_params->NSflag2 = input_params->NSflag2;
	waveform_params->shift_time = false;
	waveform_params->shift_phase = input_params->shift_phase;
	waveform_params->LISA_alpha0 = input_params->LISA_alpha0;
	waveform_params->LISA_phi0 = input_params->LISA_phi0;
	if( check_mod(gen_method)){
		waveform_params->bppe = input_params->bppe;
		waveform_params->Nmod = input_params->Nmod;
		waveform_params->betappe = new T[waveform_params->Nmod];
	}

}
template void repack_non_parameter_options<double>(gen_params_base<double> *, gen_params_base<double> *, std::string);
template void repack_non_parameter_options<adouble>(gen_params_base<adouble> *, gen_params_base<double> *, std::string);

template<class T>
void deallocate_non_param_options(gen_params_base<T> *waveform_params, gen_params_base<double> *input_params, std::string gen_method)
{
	if( check_mod(gen_method)){
		delete [] waveform_params->betappe	;
	}
}
template void deallocate_non_param_options<double>(gen_params_base<double> *, gen_params_base<double> *, std::string);
template void deallocate_non_param_options<adouble>(gen_params_base<adouble> *, gen_params_base<double> *, std::string);


/*! \brief Subroutine to calculate fisher elements for a subset of the fisher
 *
 * Skips elements that have dimensions (i,j) for i!=j && i>base_dim && j>base_dim
 *
 * Sets non-computed elements to zero
 *
 */
void calculate_fisher_elements_batch(double *frequency, 
	int length, 
	int base_dimension, 
	int full_dimension, 
	std::complex<double> **response_deriv, 
	double **output,
	double *psd)
{
	//list of modifications
	int mod_list[full_dimension-base_dimension];
	for(int i = base_dimension; i<full_dimension; i++){
		mod_list[i-base_dimension]=i;
	}

	double *integrand=new double[length];
	for (int j=0;j<full_dimension; j++)
	{
		
		for (int k = 0; k<j; k++)
		{
			//Mod mod element of fisher, set to zero
			if(check_list(k, mod_list, full_dimension-base_dimension) &&check_list(k, mod_list, full_dimension-base_dimension)){
				output[j][k] = 0;	
			}
			//GR GR or GR mod element of the fisher
			else{
				for (int i =0;i<length;i++)
				{
					integrand[i] = 
						real( 
						(response_deriv[j][i]*
						std::conj(response_deriv[k][i]))
						/psd[i]);
				}
				
				output[j][k] = 4*simpsons_sum(
							frequency[1]-frequency[0], length, integrand);	
			}
			output[k][j] = output[j][k];
		}

	}

	//All diagonal terms are calculated
	for (int j = 0; j<full_dimension; j ++)
	{

		for (int i =0;i<length;i++){
				integrand[i] = 
					real( (response_deriv[j][i]*std::conj(response_deriv[j][i]))
					/psd[i]);
		}
		
		output[j][j] = 4*simpsons_sum(
					frequency[1]-frequency[0], length, integrand);	
	}
	delete [] integrand;	

}
void calculate_fisher_elements(double *frequency, 
	int length, 
	int dimension, 
	std::complex<double> **response_deriv, 
	double **output,
	double *psd,
	std::string integration_method, 
	double *weights, 
	bool log10_f)
{
	double *integrand=new double[length];
	for (int j=0;j<dimension; j++)
	{
		for (int k = 0; k<j; k++)
		{
			for (int i =0;i<length;i++)
			{
				integrand[i] = 
					real( 
					(response_deriv[j][i]*
					std::conj(response_deriv[k][i]))
					/psd[i]);
				//Jacobian for integrating in logspace
				if(log10_f){
					integrand[i]*=(frequency[i])*LOG10;
				}
			}
			if(integration_method =="GAUSSLEG"){
				double sum = 0;
				for(int i = 0 ; i<length; i++){
					sum+= ( weights[i] * integrand[i]);
	
				}
				output[j][k] = 4*sum;	
			}
			else if(integration_method =="SIMPSONS"){
			
				output[j][k] = 4*simpsons_sum(
							frequency[1]-frequency[0], length, integrand);	
			}
			output[k][j] = output[j][k];
		}

	}

	for (int j = 0; j<dimension; j ++)
	{

		for (int i =0;i<length;i++){
				integrand[i] = 
					real( (response_deriv[j][i]*std::conj(response_deriv[j][i]))
					/psd[i]);
				//Jacobian for integrating in logspace
				if(log10_f){
					integrand[i]*=(frequency[i])*LOG10;
				}
		}
			if(integration_method =="GAUSSLEG"){
				double sum = 0;
				for(int i = 0 ; i<length; i++){
					sum+= ( weights[i] * integrand[i]);
	
				}
				output[j][j] = 4*sum;	
			}
			else if(integration_method =="SIMPSONS"){
				output[j][j] = 4*simpsons_sum(
					frequency[1]-frequency[0], length, integrand);	
			}
	}
	delete [] integrand;	

}
//#################################################################
template void repack_parameters<adouble>(adouble *, gen_params_base<adouble> *, std::string, int, gen_params_base<double> *);
template void repack_parameters<double>(double *, gen_params_base<double> *, std::string, int, gen_params_base<double> *);

void prep_gsl_subroutine(gsl_subroutine *params_packed)
{
	int dimension = params_packed->dim;
	std::string generation_method = params_packed->generation_method;
	std::string detector = params_packed->detector;

	params_packed->boundary_num = boundary_number(generation_method);

	int boundary_num = params_packed->boundary_num;
	params_packed->freq_boundaries = new double[boundary_num];
	params_packed->grad_freqs = new double[boundary_num];
	params_packed->log_factors = new bool[dimension];
	params_packed->waveform_tapes = new int[boundary_num];
	if(detector=="LISA"){
		params_packed->time_tapes = new int[boundary_num];
		params_packed->phase_tapes = new int[boundary_num];
	}
	double *freq_boundaries = params_packed->freq_boundaries;
	double *grad_freqs = params_packed->grad_freqs;
	bool *log_factors = params_packed->log_factors;
	gen_params *params = params_packed->gen_params_in;
	if(boundary_num == -1){
		std::cout<<"Error -- unsupported generation method"<<std::endl;
		exit(1);
	}
	int vec_param_length = dimension+1;
	double vec_parameters[vec_param_length];
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);
	//unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
}
void tape_phase_gsl_subroutine(gsl_subroutine * params_packed)
{
	int boundary_num = params_packed->boundary_num;
	double *freq_boundaries = params_packed->freq_boundaries;
	double *grad_freqs = params_packed->grad_freqs;
	std::string detector = params_packed->detector;
	std::string generation_method = params_packed->generation_method;
	gen_params *params = params_packed->gen_params_in;
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	std::string local_gen_method = local_generation_method(generation_method);
	
	//calculate derivative of phase
	int *tapes = params_packed->phase_tapes;
	gen_params_base<adouble> aparams;
	transform_parameters(params, &aparams);
	for(int i = 0 ; i<boundary_num ; i++){
		tapes[i]=(i+1)*31;	
		trace_on(tapes[i]);
		adouble freq;
		freq <<= grad_freqs[i];
		adouble phasep, phasec;
		fourier_phase(&freq, 1, &phasep, &phasec, generation_method, &aparams);
		double phaseout;
		phasep>>=phaseout;
		trace_off();
	}
	if(check_mod(generation_method)){
		delete [] aparams.betappe;
		delete [] aparams.bppe;
	}

}
void tape_time_gsl_subroutine(gsl_subroutine * params_packed)
{
	int id1 = params_packed->id1;
	int id2 = params_packed->id2;
	int dimension = params_packed->dim;
	int boundary_num = params_packed->boundary_num;
	double *freq_boundaries = params_packed->freq_boundaries;
	double *grad_freqs = params_packed->grad_freqs;
	bool *log_factors = params_packed->log_factors;
	std::string detector = params_packed->detector;
	std::string generation_method = params_packed->generation_method;
	gen_params *params = params_packed->gen_params_in;
	//calculate hessian of phase, take [0][j] components to get the derivative of time
	int vec_param_length = dimension +1 ;//+1 for frequency 
	std::string local_gen_method = local_generation_method(generation_method);
	double vec_parameters[vec_param_length];
	unpack_parameters(&vec_parameters[1], params, generation_method, dimension, log_factors);
	
	//calculate derivative of phase
	int *tapes = params_packed->time_tapes;
	for(int i = 0 ; i < boundary_num ; i++){
		tapes[i] = (i+1)*17; //Random tape id 
		trace_on(tapes[i]);
		adouble avec_parameters[vec_param_length];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			if(j == id1+1 || j == id2+1){
				avec_parameters[j]<<=vec_parameters[j];	
			}
			else{
				avec_parameters[j]=vec_parameters[j];	
			}
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,params,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension,params);
		adouble phasep, phasec;
		int status  = fourier_phase(&afreq, 1, &phasep,&phasec, local_gen_method, &a_parameters);
		double phase;
		phasep >>= phase;

		trace_off();
		deallocate_non_param_options(&a_parameters, params, generation_method);
	}

}
void tape_waveform_gsl_subroutine(gsl_subroutine * params_packed)
{
	int id1 = params_packed->id1;
	int id2 = params_packed->id2;
	int dimension = params_packed->dim;
	int boundary_num = params_packed->boundary_num;
	double *freq_boundaries = params_packed->freq_boundaries;
	double *grad_freqs = params_packed->grad_freqs;
	bool *log_factors = params_packed->log_factors;
	std::string generation_method = params_packed->generation_method;
	std::string detector = params_packed->detector;
	gen_params *parameters = params_packed->gen_params_in;

	int vec_param_length= dimension +1;
	if(detector == "LISA"){
		//take derivative wrt time as well, for the chain rule
		vec_param_length += 1;
	}
	double vec_parameters[vec_param_length];
	std::string local_gen_method = local_generation_method(generation_method);
	vec_parameters[0]=grad_freqs[0];
	unpack_parameters(&vec_parameters[1], parameters, generation_method,dimension, log_factors);
	double *grad_times=NULL;
	if(detector == "LISA"){
		grad_times = new double[boundary_num];
		time_phase_corrected_autodiff(grad_times, boundary_num, grad_freqs, parameters, generation_method, false);
			
	}
	//calculate_derivative tapes
	int *tapes = params_packed->waveform_tapes;
	for(int i =0; i<boundary_num; i++){
		tapes[i]= (i+1)*7;
		trace_on(tapes[i]);
		//adouble avec_parameters[dimension+1];
		adouble avec_parameters[vec_param_length];
		//avec_parameters[0] =grad_freqs[i];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			if(j == id1+1 || j == id2+1){
				avec_parameters[j]<<=vec_parameters[j];	
			}
			else{
				avec_parameters[j]=vec_parameters[j];	
			}
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		adouble afreq;
		afreq = avec_parameters[0];
		//############################################
		//Non variable parameters
		repack_non_parameter_options(&a_parameters,parameters,generation_method);
		//############################################
		repack_parameters(&avec_parameters[1],&a_parameters,generation_method, dimension, parameters);
		adouble time;
		if(detector == "LISA"){
			time <<= grad_times[i];
			//map_extrinsic_angles(&a_parameters);
		}
		std::complex<adouble> a_response;
		if(!a_parameters.sky_average){
			int status  = fourier_detector_response(&afreq, 1, &a_response, detector, local_gen_method, &a_parameters, &time);

		}
		else{
			adouble a_amp;
			adouble a_phasep;
			adouble a_phasec;
		
			int status  = fourier_amplitude(&afreq, 1, &a_amp, local_gen_method, &a_parameters);
			status  = fourier_phase(&afreq, 1, &a_phasep,  &a_phasec,local_gen_method, &a_parameters);
			a_response = a_amp * exp(std::complex<adouble>(0,a_phasep));

		}
		double response[2];
		real(a_response) >>=  response[0];	
		imag(a_response) >>=  response[1];	

		trace_off();
		deallocate_non_param_options(&a_parameters, parameters, generation_method);
	}
	if(detector=="LISA"){
		delete [] grad_times;
	}
}
/*! \brief Routine that implements GSL numerical integration to calculate the Fishers
 *
 * This can be faster than brute force calculations in fisher_autodiff, but that depends
 *
 * Trade offs: 
 *
 * Every element is calculated independently, so no information is retained between elements. In the brute force calculation, there is information reused.
 *
 * However, time can be saved by spending less time on trivial elements (identically 0 elements, etc) and better spent on complicated elements
 *
 * Does not have a direct interpretation in terms of integration time, as the scheme is adaptative. Sampling frequency and integration time are ``as good as they need to be'' to calculate the fisher
 *
 * Implements (GSL_INTEG_GAUSS15)
 */
void fisher_autodiff_gsl_integration(double *frequency_bounds, /**<Bounds of integration in fourier space*/
	string generation_method, /**<Method of waveform generation*/
	string sensitivity_curve, /**<Sensitivity curve to be used for the PSD -- MUST BE ANALYTIC*/
	string detector, /**< Detector to use for the response function*/
	double **output,/**<[out] Output Fisher -- must be preallocated -- shape [dimension][dimension]*/
	double **error,/**<[out] Estimated error, as specified by GSL's integration -- must be preallocated -- shape [dimension][dimension]*/
	int dimension, /**<Dimension of the Fisher */
	gen_params *parameters,/**< Generation parameters specifying source parameters and waveform options*/
	double abserr,/**<Target absolute error (0 if this should be ignored -- ONE type of error must be specified)*/
	double relerr/**<Target relative error (0 if this should be ignored -- ONE type of error must be specified)*/
	)
{
	fisher_autodiff_gsl_integration(frequency_bounds, generation_method, sensitivity_curve, detector, output,error, dimension, parameters, abserr, relerr, "", false);
}

/*! \brief Routine that implements GSL numerical integration to calculate the Fishers
 *
 * This can be faster than brute force calculations in fisher_autodiff, but that depends
 *
 * Trade offs: 
 *
 * Every element is calculated independently, so no information is retained between elements. In the brute force calculation, there is information reused.
 *
 * However, time can be saved by spending less time on trivial elements (identically 0 elements, etc) and better spent on complicated elements
 *
 * Does not have a direct interpretation in terms of integration time, as the scheme is adaptative. Sampling frequency and integration time are ``as good as they need to be'' to calculate the fisher
 *
 * Implements (GSL_INTEG_GAUSS15)
 *
 * Now includes option to log error instead of ending program for certain types of errors
 */
void fisher_autodiff_gsl_integration(double *frequency_bounds, /**<Bounds of integration in fourier space*/
	string generation_method, /**<Method of waveform generation*/
	string sensitivity_curve, /**<Sensitivity curve to be used for the PSD -- MUST BE ANALYTIC*/
	string detector, /**< Detector to use for the response function*/
	double **output,/**<[out] Output Fisher -- must be preallocated -- shape [dimension][dimension]*/
	double **error,/**<[out] Estimated error, as specified by GSL's integration -- must be preallocated -- shape [dimension][dimension]*/
	int dimension, /**<Dimension of the Fisher */
	gen_params *parameters,/**< Generation parameters specifying source parameters and waveform options*/
	double abserr,/**<Target absolute error (0 if this should be ignored -- ONE type of error must be specified)*/
	double relerr,/**<Target relative error (0 if this should be ignored -- ONE type of error must be specified)*/
	std::string error_log,/**< File to write non-critical error codes to (roundoff error)*/
	bool logerr/**<Whether or not to end program with certain error codes, or to log them and continue*/
	)
{
	gsl_subroutine params_packed ;
	params_packed.detector = detector;
	params_packed.sensitivity_curve = sensitivity_curve;
	params_packed.generation_method = generation_method;
	params_packed.gen_params_in = parameters;
	params_packed.dim = dimension;
	gsl_function F;
	//double abserr = 0, relerr = 1e-3;
	double result;	
	double err;	
	int id1,id2;
	size_t np = 1e5;	//Max number of division
	std::ofstream log_file;
	//int *err_rows;
	//int *err_cols;
	bool no_err=true;
	if(logerr){
		log_file.open(error_log);	
		//err_rows = new int [dimension];
		//err_cols = new int [dimension];
		
	}
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);
	for(int i = 0 ; i<dimension; i++){
		for(int j = 0 ; j<=i ; j++){
			//i>j
			params_packed.id1 = j;
			params_packed.id2 = i;
			//Find boundary_num, boundary_freqs, grad_freqs
			prep_gsl_subroutine(&params_packed);
			//Tape functions
			tape_waveform_gsl_subroutine( &params_packed);
			if(detector == "LISA"){
				tape_time_gsl_subroutine(&params_packed);
				tape_phase_gsl_subroutine(&params_packed);
			}


			F.function = &calculate_integrand_autodiff_gsl_subroutine;
			F.params = &params_packed;

			//std::cout<<"Integrating "<<i<<" "<<j<<std::endl;
			if(logerr){
				gsl_set_error_handler_off();
			}
			int errcode = gsl_integration_qag(&F, frequency_bounds[0],frequency_bounds[1], abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result, &err);
			if(logerr && errcode){
				if(errcode == 18){
					log_file<<i<<" "<<j<<" : "<<gsl_strerror(errcode);
					no_err = false;
					//err_rows[
				}
				else{
					std::cout<<"Error -- not roundoff"<<std::endl;;
					exit(-1);
				}
			}	

			output[i][j] = 4*result;
			output[j][i] = output[i][j];
			error[i][j] = err;
			error[j][i] = error[i][j];
		
			delete [] params_packed.waveform_tapes;
			if(detector=="LISA"){
				delete [] params_packed.time_tapes;
				delete [] params_packed.phase_tapes;
			}
			delete [] params_packed.freq_boundaries ;
			delete [] params_packed.grad_freqs ;
			delete [] params_packed.log_factors ;

			//std::cout<<"2*Result "<<2*output[i][j]<<std::endl;
			//std::cout<<"Error "<<err<<std::endl;
			//std::cout<<"intervals "<<w->size<<std::endl;
		}
	}
	if(logerr){
		log_file.close();
	}
	if(logerr && no_err){
		std::remove(error_log.c_str());
	}

	gsl_integration_workspace_free (w);

	if(detector == "LISA"){
		for(int i = 0 ; i<dimension;i++){
			for(int j = 0  ;j<dimension; j++){
				output[i][j]*=2;
			}	
		}
	}
}
/*! \brief Routine that implements GSL numerical integration to calculate the Fishers -- batch modifications version
 *
 * Calculates Fisher for multiple modifications at a time, neglecting covariance between modifications (set to 0 in Fisher)
 *
 * Modifications MUST BE evaluated at 0 for this routine to calculate correct results
 *
 * This can be faster than brute force calculations in fisher_autodiff, but that depends
 *
 * Trade offs: 
 *
 * Every element is calculated independently, so no information is retained between elements. In the brute force calculation, there is information reused.
 *
 * However, time can be saved by spending less time on trivial elements (identically 0 elements, etc) and better spent on complicated elements
 *
 * Does not have a direct interpretation in terms of integration time, as the scheme is adaptative. Sampling frequency and integration time are ``as good as they need to be'' to calculate the fisher
 *
 * Implements (GSL_INTEG_GAUSS15)
 */
void fisher_autodiff_gsl_integration_batch_mod(double *frequency_bounds, /**<Bounds of integration in fourier space*/
	string generation_method, /**<Method of waveform generation*/
	string sensitivity_curve, /**<Sensitivity curve to be used for the PSD -- MUST BE ANALYTIC*/
	string detector, /**< Detector to use for the response function*/
	double **output,/**<[out] Output Fisher -- must be preallocated -- shape [full_dimension][full_dimension]*/
	double **error,/**<[out] Estimated error, as specified by GSL's integration -- must be preallocated -- shape [full_dimension][full_dimension]*/
	int base_dimension, /**< Dimension of base model (ie GR dimension)*/
	int full_dimension, /**< Full dimension (GR dimension + Nmod)*/
	gen_params *parameters,/**< Generation parameters specifying source parameters and waveform options*/
	double abserr,/**<Target absolute error (0 if this should be ignored -- ONE type of error must be specified)*/
	double relerr/**<Target relative error (0 if this should be ignored -- ONE type of error must be specified)*/
	)
{
	fisher_autodiff_gsl_integration_batch_mod(frequency_bounds, generation_method,
		sensitivity_curve, detector, output, error, base_dimension, full_dimension, parameters, abserr, relerr, "",false);
}
/*! \brief Routine that implements GSL numerical integration to calculate the Fishers -- batch modifications version
 *
 * Calculates Fisher for multiple modifications at a time, neglecting covariance between modifications (set to 0 in Fisher)
 *
 * Modifications MUST BE evaluated at 0 for this routine to calculate correct results
 *
 * This can be faster than brute force calculations in fisher_autodiff, but that depends
 *
 * Trade offs: 
 *
 * Every element is calculated independently, so no information is retained between elements. In the brute force calculation, there is information reused.
 *
 * However, time can be saved by spending less time on trivial elements (identically 0 elements, etc) and better spent on complicated elements
 *
 * Does not have a direct interpretation in terms of integration time, as the scheme is adaptative. Sampling frequency and integration time are ``as good as they need to be'' to calculate the fisher
 *
 * Implements (GSL_INTEG_GAUSS15)
 *
 * Now includes option to log error instead of ending program for certain types of errors
 */
void fisher_autodiff_gsl_integration_batch_mod(double *frequency_bounds, /**<Bounds of integration in fourier space*/
	string generation_method, /**<Method of waveform generation*/
	string sensitivity_curve, /**<Sensitivity curve to be used for the PSD -- MUST BE ANALYTIC*/
	string detector, /**< Detector to use for the response function*/
	double **output,/**<[out] Output Fisher -- must be preallocated -- shape [full_dimension][full_dimension]*/
	double **error,/**<[out] Estimated error, as specified by GSL's integration -- must be preallocated -- shape [full_dimension][full_dimension]*/
	int base_dimension, /**< Dimension of base model (ie GR dimension)*/
	int full_dimension, /**< Full dimension (GR dimension + Nmod)*/
	gen_params *parameters,/**< Generation parameters specifying source parameters and waveform options*/
	double abserr,/**<Target absolute error (0 if this should be ignored -- ONE type of error must be specified)*/
	double relerr,/**<Target relative error (0 if this should be ignored -- ONE type of error must be specified)*/
	std::string error_log,/**< File to write non-critical error codes to (roundoff error)*/
	bool logerr/**<Whether or not to end program with certain error codes, or to log them and continue*/
	)
{
	int mod_list[full_dimension-base_dimension];
	for(int i = base_dimension; i<full_dimension; i++){
		mod_list[i-base_dimension]=i;
	}
	gsl_subroutine params_packed ;
	params_packed.detector = detector;
	params_packed.sensitivity_curve = sensitivity_curve;
	params_packed.generation_method = generation_method;
	params_packed.gen_params_in = parameters;
	params_packed.dim = full_dimension;
	gsl_function F;
	//double abserr = 0, relerr = 1e-3;
	double result;	
	double err;	
	int id1,id2;
	size_t np = 1e5;	//Max number of division
	std::ofstream log_file;
	bool no_err = true;
	if(logerr){
		log_file.open(error_log);	
	}
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);
	for(int i = 0 ; i<full_dimension; i++){
		for(int j = 0 ; j<=i ; j++){
			//i>j
			if(check_list(i, mod_list, full_dimension-base_dimension) && check_list(j, mod_list, full_dimension-base_dimension) && i!=j){
				output[i][j] = 0;	
				output[j][i] = 0;	
				error[i][j] = 0;	
				error[j][i] = 0;	
			}
			else{
				params_packed.id1 = j;
				params_packed.id2 = i;
				//Find boundary_num, boundary_freqs, grad_freqs
				prep_gsl_subroutine(&params_packed);
				//Tape functions
				tape_waveform_gsl_subroutine( &params_packed);
				if(detector == "LISA"){
					tape_time_gsl_subroutine(&params_packed);
					tape_phase_gsl_subroutine(&params_packed);
				}


				F.function = &calculate_integrand_autodiff_gsl_subroutine;
				F.params = &params_packed;

				//std::cout<<"Integrating "<<i<<" "<<j<<std::endl;
				if(logerr){
					gsl_set_error_handler_off();
				}
				int errcode = gsl_integration_qag(&F, frequency_bounds[0],frequency_bounds[1], abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result, &err);
				if(logerr && errcode){
					if(errcode == 18){
						log_file<<i<<" "<<j<<" : "<<gsl_strerror(errcode)<<std::endl;
						no_err = false;
					}
					else{
						std::cout<<"Error -- not roundoff"<<std::endl;;
						exit(-1);
					}
				}	
				output[i][j] = 4*result;
				output[j][i] = output[i][j];
				error[i][j] = err;
				error[j][i] = error[i][j];
		
				delete [] params_packed.waveform_tapes;
				if(detector=="LISA"){
					delete [] params_packed.time_tapes;
					delete [] params_packed.phase_tapes;
				}
				delete [] params_packed.freq_boundaries ;
				delete [] params_packed.grad_freqs ;
				delete [] params_packed.log_factors ;
				//std::cout<<"2*Result "<<2*output[i][j]<<std::endl;
				//std::cout<<"Error "<<err<<std::endl;
				//std::cout<<"intervals "<<w->size<<std::endl;

			}
		}
	}
	if(logerr){
		log_file.close();
	}
	if(logerr && no_err){
		std::remove(error_log.c_str());
	}

	gsl_integration_workspace_free (w);

	if(detector == "LISA"){
		//std::cout<<"Doubling"<<std::endl;
		for(int i = 0 ; i<full_dimension;i++){
			for(int j = 0  ;j<full_dimension; j++){
				output[i][j]*=2;
			}	
		}
	}


}
/*! \brief Calculates the derivatives of the detector response using automatic differentiation -- one frequency for gsl_integration
 *
 * Possibly slower than the numerical derivative, but not susceptible to truncation error from finite difference
 *
 * Higher dimensional fishers actually could be faster
 *
 * NOTE: dimension parameter ALWAYS refers to the dimension of the fisher (ie the length of the source parameter vector), even though the derivatives are computed wrt dimension +1 or dimension + 2 -- the +1(+2) are for the frequency deriv(time deriv)
 */
double calculate_integrand_autodiff_gsl_subroutine(double frequency, void *params_in)
{
	//double start = clock();
	gsl_subroutine params_packed = *(gsl_subroutine *)params_in;
	std::string detector=  params_packed.detector;
	std::string generation_method=  params_packed.generation_method;
	gen_params *parameters = params_packed.gen_params_in;
	std::string sensitivity_curve = params_packed.sensitivity_curve;
	int * wf_tapes = params_packed.waveform_tapes;
	int * time_tapes = params_packed.time_tapes;
	int boundary_num = params_packed.boundary_num;
	double *freq_boundaries = params_packed.freq_boundaries;
	bool *log_factors = params_packed.log_factors;
	int id1 = params_packed.id1;
	int id2 = params_packed.id2;
	int dimension = params_packed.dim;
	//std::cout<<id1<<" "<<id2<<std::endl;

	//Transform gen_params to double vectors
	int vec_param_length= dimension +1;
	int indep ;
	if(id1!=id2){
		indep= 3;//2 plus frequency
	}
	else{
		indep = 2;//1 plus frequency
	}
	if(detector == "LISA"){
		//take derivative wrt time as well, for the chain rule
		vec_param_length += 1;
		indep +=1;
	}
	std::complex<double> waveform_deriv[2];
	double vec_parameters[vec_param_length];
	unpack_parameters(&vec_parameters[1], parameters, generation_method,dimension, log_factors);
	double eval_time;
	//std::cout<<"Prep time: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	//start = clock();
	if(detector == "LISA"){
		time_phase_corrected_autodiff(&eval_time, 1, &frequency, parameters, generation_method, false, params_packed.phase_tapes);
	}
	//std::cout<<"Time calc time: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	//start = clock();
	//Evaluate derivative tapes
	int dep = 2;//Output is complex
	bool eval = false;//Keep track of when a boundary is hit
	double indep_vec[indep];
	double **jacob = allocate_2D_array(dep,indep);
	std::complex<double> dt;
	vec_parameters[0]=frequency;
	indep_vec[0] = frequency;
	indep_vec[1] = vec_parameters[id1+1];
	if(id1!=id2){
		indep_vec[2] = vec_parameters[id2+1];
	}
	if(detector == "LISA"){
		indep_vec[indep -1] = eval_time;
	}
	//std::cout<<indep_vec[0]<<" "<<indep_vec[1]<<" "<<indep_vec[2]<<" "<<indep_vec[3]<<" "<<std::endl;
	double **dt_hess;
	for(int n = 0 ; n<boundary_num; n++){
		if(indep_vec[0]<freq_boundaries[n]){
			jacobian(wf_tapes[n], dep, indep, indep_vec, jacob);
			
			if(detector=="LISA"){
			//if(false){
				dt_hess = new double*[indep-1];
				for(int j = 0 ; j<indep-1; j++){
					dt_hess[j]=new double[indep-1];
				}
				double time_eval[indep-1];
				time_eval[0]=indep_vec[0];
				time_eval[1]=indep_vec[1];
				if(id1!=id2){
					time_eval[2]=indep_vec[2];
				}
				hessian(time_tapes[n], indep-1, time_eval, dt_hess);
				//if(id1!=id2)	
				//std::cout<<dt_hess[1][0]<<" "<<dt_hess[2][0]<<std::endl;
				//else
				//std::cout<<dt_hess[1][0]<<std::endl;
			}
			int derivs;
			if(id1!=id2){derivs =2;}
			else{derivs =1;}
			for(int i =0; i<derivs; i++){
				waveform_deriv[i] = jacob[0][i+1] 
					+ std::complex<double>(0,1)*jacob[1][i+1];
				//correct for time deriv for LISA
				if(detector == "LISA" ){
				//if(false ){
					waveform_deriv[i]+= 
						(jacob[0][indep-1] + std::complex<double>(0,1)*jacob[1][indep-1]) //Time derivative of WF
						* dt_hess[i+1][0]/(2*M_PI);//Derivative of time wrt source parameter -- 2 pi from definition of time 
				}
			}
			if(detector=="LISA"){
			//if(false){
				for(int j = 0 ; j<indep-1; j++){
					delete [] dt_hess[j];
				}
				delete [] dt_hess;
			}
			//Mark successful derivative
			eval = true;
			if(id1==id2){waveform_deriv[1]=waveform_deriv[0];}
			//Skip the rest of the bins
			break;
		}
	}
	//std::cout<<"Eval calc time: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	//start = clock();
	//If freq didn't fall in any boundary, set to 0
	if(!eval){
		for(int i =0; i<indep; i++){
			waveform_deriv[i] = std::complex<double>(0,0);
		}	
	}
	//Account for Log parameters
	if(log_factors[id1]){
		waveform_deriv[0] *= vec_parameters[id1+1];
	}
	//if(id1!=id2){
	//This needs to be here even if id1==id2
	if(log_factors[id2]){
		waveform_deriv[1] *= vec_parameters[id2+1];
	}
	deallocate_2D_array(jacob,dep,indep);
	double psdroot;
	populate_noise(&frequency,sensitivity_curve, &psdroot, 1);
	//std::cout<<std::real((waveform_deriv[0]*conj(waveform_deriv[1])))/(psdroot*psdroot)<<std::endl;
	return std::real((waveform_deriv[0]*conj(waveform_deriv[1])))/(psdroot*psdroot);
}
