#include "waveform_util.h"
#include "util.h"
#include "waveform_generator.h"
#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomP.h"
#include "detector_util.h"
#include <fftw3.h>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <adolc/taping.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
/*!\file 
 * Utilities for waveforms - SNR calculation and detector response
 * 	
 * includes snr and detector response
 *
 * Includes some utilities useful for MCMC and fisher calculations, as well as time-frequency methods for detectors like LISA
 */


/*! \brief Utility to calculate the snr of a fourier transformed data stream while maximizing over the coalescence parameters phic and tc
 *
 * The gen_params structure holds the parameters for the template to be used (the maximimum likelihood parameters)
 */
double data_snr_maximized_extrinsic(double *frequencies, /**< Frequencies used by data*/
				int length,/**< length of the data*/
				std::complex<double> *data,/**< input data in the fourier domain*/
				double *psd,/**< PSD for the detector that created the data*/
				std::string detector,/**< Name of the detector --See noise_util for options */
				std::string generation_method,/**< Generation method for the template -- See waveform_generation.cpp for options*/
				gen_params *param/**< gen_params structure for the template*/
				)
{
	
	/*produce noise curve*/
	double *noise = (double *)malloc(sizeof(double)*length);	
	//populate_noise(frequencies,detector, noise, length);
	for(int i = 0; i<length; i++)
		//noise[i] = noise[i]*noise[i];
		noise[i ] =psd[i];

	std::complex<double> *detector_response
			 = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	double *integrand
			 = (double *)malloc(sizeof(double)*length);
	/*produce the waveform*/
	fourier_detector_response(frequencies, length, detector_response, detector,generation_method, param);

	/*Calculate the template snr integrand 4*Re(h* h /S(f)) - factor of q for the plus, cross modes 
 * 	effect on the detector*/
	for (int i = 0; i<length;i++)
		integrand[i] = 4.*real(conj(detector_response[i])*detector_response[i])/noise[i]; 
	double delta_f = frequencies[1]-frequencies[0];
	double snr_template;
	snr_template = sqrt(simpsons_sum(delta_f,length, integrand));

	fftw_complex *in, *out; 
	fftw_plan p;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        p = fftw_plan_dft_1d(length, in, out,FFTW_FORWARD, FFTW_MEASURE);
	std::complex<double> g_tilde;
        for (int i=0;i<length; i++)
        {
                g_tilde = 4.*conj(data[i]) * detector_response[i] / noise[i];
                in[i][0] = real(g_tilde);
                in[i][1] = imag(g_tilde);
        }

        double *g = (double *)malloc(sizeof(double) * length);

        fftw_execute(p);

        for (int i=0;i<length; i++)
        {
                g[i] = std::abs(std::complex<double>(out[i][0],out[i][1])) ;
        }

        double max = *std::max_element(g, g+length)*delta_f;
	


	fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
	fftw_cleanup();

	free(g);
	free(noise);
	free(detector_response);
	free(integrand);
	return max/(snr_template);
}
/*! \brief Light wrapper for the data_snr_maximized_extrinsic method
 *
 * Splits the data into real and imaginary, so all the arguments are C-safe
 */
double data_snr_maximized_extrinsic(double *frequencies, /**< Frequencies used by data*/
				int length,/**< length of the data*/
				double *data_real,/**< input data in the fourier domain -- real part*/
				double *data_imag,/**< input data in the fourier domain -- imaginary part*/
				double *psd,/**< PSD for the detector that created the data*/
				std::string detector,/**< Name of the detector --See noise_util for options */
				std::string generation_method,/**< Generation method for the template -- See waveform_generation.cpp for options*/
				gen_params *param/**< gen_params structure for the template*/
				)
{
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>) * length);
        for (int i =0; i<length; i++)
                data[i] = std::complex<double>(data_real[i],data_imag[i]);	
	double snr;
	snr = data_snr_maximized_extrinsic(frequencies,
				length,
				data,
				psd,
				detector,
				generation_method,
				param);
	free(data);
	return snr;
}
/*! \brief Caclulates the snr given a detector and waveform (complex) and frequencies
 *      
 * This function computes the un-normalized snr: \sqrt( ( H | H ) )
 */     
double calculate_snr(std::string detector, /**< detector name - must match the string of populate_noise precisely*/
                        std::complex<double> *waveform,/**< complex waveform */
                        double *frequencies,/**< double array of frequencies that the waveform is evaluated at*/
                        int length/**< length of the above two arrays*/
                        )
{
        double *noise = (double *)malloc(sizeof(double)*length);
        populate_noise(frequencies,detector, noise,  length);
        for (int i = 0; i< length; i++)
                noise[i] = noise[i]*noise[i];
        double *integrand = (double *) malloc(sizeof(double)*length);
        for (int i = 0; i<length; i++)
                integrand[i] = 4.* real(conj(waveform[i])*waveform[i]/noise[i]);
        double integral = trapezoidal_sum(frequencies, length, integrand);
	free(integrand);
	free(noise);
        return sqrt(integral);
}

/* \brief calculates the detector response for a given waveform and detector -- polarization angle =0
 */
template<class T>
int fourier_detector_response(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			std::complex<T> *detector_response, /**< [out] detector response*/
			T theta, /**< polar angle (rad) theta in detector frame*/
			T phi, /**< azimuthal angle (rad) phi in detector frame*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	return fourier_detector_response(frequencies, length, hplus, hcross, detector_response, theta, phi, (T)0., detector);
	
}
/* \brief calculates the detector response for a given waveform and detector
 */
template<class T>
int fourier_detector_response(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			std::complex<T> *detector_response, /**< [out] detector response*/
			T theta, /**< polar angle (rad) theta in detector frame*/
			T phi, /**< azimuthal angle (rad) phi in detector frame*/ 
			T psi, /**< polarization angle (rad) phi in detector frame*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	T fplus, fcross, Fplus, Fcross, c2psi, s2psi;
	
	if(	detector == "LIGO" || 
		detector == "Livingston" || 
		detector == "LIVINGSTON" || 
		detector == "livingston" || 
		detector == "Hanford" || 
		detector == "HANFORD" || 
		detector == "hanford" || 
		detector == "VIRGO" ||
		detector == "Virgo" ||
		detector == "virgo")
	{
		fplus = right_interferometer_plus(theta,phi);
		fcross = right_interferometer_cross(theta,phi);	
		c2psi = cos(2*psi);
		s2psi = sin(2*psi);
		Fplus = fplus*c2psi- fcross*s2psi;
		Fcross = fplus*s2psi+ fcross*s2psi;
	}
	for (int i =0; i <length; i++)
	{
		detector_response[i] = Fplus * hplus[i] 
					+ (Fcross )*hcross[i];
	}	
	return status;
	
}

/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
	int length,/**<integer length of all the arrays*/
	std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
	std::string detector,
	std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
	gen_params_base<T> *parameters/**<structure containing all the source parameters*/
	)
{
	int status = 1;
	//generate waveform
	std::complex<T> *waveform_plus =
		(std::complex<T> *)malloc(sizeof(std::complex<T>) * length);
	std::complex<T> *waveform_cross=
		(std::complex<T> *)malloc(sizeof(std::complex<T>) * length);
	status = fourier_waveform(frequencies, 
			length,
			waveform_plus, 
			waveform_cross, 
			generation_method,
			parameters
			);
	status = fourier_detector_response(frequencies, 
			length, 
			waveform_plus, 
			waveform_cross,
			response, 
			parameters->theta, 
			parameters->phi, 
			parameters->psi, 
			detector 
			) ;
	
		
	//Deallocate memory
	free(waveform_plus);
	free(waveform_cross);

	return status;
}
/* \brief calculates the detector response for a given waveform and detector -- using equatorial coordinates and greenwich mean sidereal time
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			std::complex<T> *detector_response, /**< [out] detector response*/
			T ra, /**< Right Ascension in rad*/
			T dec, /**< Declination in rad*/ 
			T psi, /**< polarization angle (rad) */ 
			double gmst, /**< greenwich mean sidereal time*/ 
			T *times,
			T LISA_alpha0,
			T LISA_phi0,
			T LISA_thetal,
			T LISA_phil,
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	//Not an elegant solution, but should work..
	T fplus ;
	T fcross ;
	T *Fplus;
	T *Fcross;
	if(detector=="LISA"){
		Fplus = new T[length];
		Fcross = new T[length];
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,LISA_thetal,LISA_phil, Fplus, Fcross);
	}
	else{
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,LISA_thetal,LISA_phil, &fplus, &fcross);
	}

	
	if(detector == "LISA"){
		for (int i =0; i <length; i++)
		{
			detector_response[i] = Fplus[i] * hplus[i] 
						+ (Fcross[i] )*hcross[i];
		}	
		delete [] Fplus; 
		delete [] Fcross;

	}
	else{
		for (int i =0; i <length; i++)
		{
			detector_response[i] = fplus * hplus[i] 
						+ (fcross )*hcross[i];
		}	
	}
	return status;
	
}
/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary for equatorial coordinates for TERRESTIAL detectors, where the earth's rotation during the signal is minor
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ 
			int length,/**<integer length of all the arrays*/
			std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters/**<structure containing all the source parameters*/
			)
{
	T *times=NULL;
	return fourier_detector_response_equatorial(frequencies, length,response, detector, generation_method, parameters, times);	

}
/*!\brief Function to produce the detector response caused by impinging gravitational waves from a quasi-circular binary for equatorial coordinates
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 * Detector options include classic interferometers like LIGO/VIRGO (coming soon: ET and LISA)
 * 	
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually 
 */
template<class T>
int fourier_detector_response_equatorial(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ 
			int length,/**<integer length of all the arrays*/
			std::complex<T> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters,/**<structure containing all the source parameters*/
			T *times
			)
{
	int status = 1;
	//generate waveform
	std::complex<T> *waveform_plus = new std::complex<T>[length];
	std::complex<T> *waveform_cross = new std::complex<T>[length];
	status = fourier_waveform(frequencies, 
			length,
			waveform_plus, 
			waveform_cross, 
			generation_method,
			parameters
			);
	status = fourier_detector_response_equatorial(frequencies, 
			length, 
			waveform_plus, 
			waveform_cross,
			response, 
			parameters->RA, 
			parameters->DEC, 
			parameters->psi, 
			parameters->gmst, 
			times,
			parameters->LISA_alpha0,
			parameters->LISA_phi0,
			parameters->LISA_thetal,
			parameters->LISA_phil,
			detector 
			) ;
	
		
	//Deallocate memory
	delete [] waveform_plus;
	delete [] waveform_cross;

	return status;
}
/*! \brief Calculates the amplitude (magnitude) and phase (argument) of the response of a given detector
 *
 * This is for general waveforms, and will work for precessing waveforms
 *
 * Not as fast as non-precessing, but that can't be helped. MUST include plus/cross polarizations
 */
int fourier_detector_amplitude_phase(double *frequencies, 
			int length,
			double *amplitude, 
			double *phase, 
			std::string detector,
			std::string generation_method,
			gen_params *parameters
			)
{
	std::complex<double> *response = (std::complex<double>*)malloc(
			sizeof(std::complex<double> ) * length);	
	fourier_detector_response(frequencies, length,
				response, detector, generation_method, parameters);
	for(int i =0 ; i< length; i++){
		amplitude[i] = std::abs(response[i]);
		phase[i] = std::arg(response[i]);
	}

	free(response);
}

/*! \brief Mapping from phase to time AUTODIFFERENTIATION using ADOL-C
 *
 * This is NOT autodiff safe. ADOL-C does not support wrapping a section of code as active twice. To find the derivative of this function, you must use the hessian function of ADOL-C
 *
 * Made for use with detectors like LISA
 *
 * Using https://arxiv.org/abs/1809.04799 t = (1/2PI) d phi/ d f
 *
 * This breaks down near merger,so at fRD, the relationship between frequency and time is extrapolated as a line
 *
 * Currently, just uses IMRPhenomD as a proxy regardless of what method is being used, as this is the analytically known function
 *
 * For IMRPhenomPv2, the phase has to be wrapped, because arctan is taken of the waveform because of the euler rotations. This might make the numerical derivative unpredictable
 *
 */
void time_phase_corrected_autodiff(double *times, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, bool correct_time)
{
	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);	
	gen_params_base<adouble> aparams;
	transform_parameters(params, &aparams);
	int tapes[boundary_num];
	for(int i = 0 ; i<boundary_num ; i++){
		tapes[i]=i*8;	
		trace_on(tapes[i]);
		adouble freq;
		freq <<= grad_freqs[i];
		adouble phasep, phasec;
		fourier_phase(&freq, 1, &phasep, &phasec, generation_method, &aparams);
		double phaseout;
		phasep>>=phaseout;
		trace_off();
	}
	
	bool eval = false;
	double freq;	
	for(int k = 0; k<length; k++){
		freq = frequencies[k];
		for(int n = 0; n<boundary_num ; n++){
			if(freq < freq_boundaries[n]){
				gradient(tapes[n], 1, &freq, &times[k]);
				//Mark successful derivative
				eval = true;
				//Skip the rest of the bins
				break;
			}
		}
		if(!eval){
			times[k]=0;
		}
		eval = false;
	}

	//divide by 2 PI
	for(int i = 0 ; i<length; i++){
		times[i]/=(2.*M_PI);
	}
	if(check_ppE(generation_method)){
		delete [] aparams.betappe;
		delete [] aparams.bppe;
	}
	
}
/*! \brief Utility to inform the fisher routine how many logical boundaries should be expected
 *
 * The automatic derivative code requires a new tape for each logical branch of the program, so each waveform_generation method needs to add the number of branches here
 */
int boundary_number(std::string method)
{
	if(method.find("IMRPhenomP") != std::string::npos || 
		method.find("IMRPhenomD")!=std::string::npos){
		return 5;
	}
	return -1;
}
/*! \brief Mapping from phase to time NUMERICAL
 *
 * Made for use with detectors like LISA
 *
 * Using https://arxiv.org/abs/1809.04799 t = (1/2PI) d phi/ d f
 *
 * This breaks down near merger,so at fRD, the relationship between frequency and time is extrapolated as a line
 *
 * Currently, just uses IMRPhenomD as a proxy regardless of what method is being used, as this is the analytically known function
 *
 * For IMRPhenomPv2, the phase has to be wrapped, because arctan is taken of the waveform because of the euler rotations. This might make the numerical derivative unpredictable
 *
 * Just uses a second order numerical derivative for now
 */
template<class T>
void time_phase_corrected(T *times, int length, T *frequencies,gen_params_base<T> *params, std::string generation_method, bool correct_time)
{
	//bool save_shift_time = params->shift_time;
	//params->shift_time = false;
	std::string local_gen = "IMRPhenomD";
	//################################################
	T *phase_plus = new T[length];
	T *phase_cross = new T[length];
	//################################################
	//################################################
	fourier_phase(frequencies, length, phase_plus, phase_cross, local_gen, params);
	//################################################
	source_parameters<T> s_param;
	s_param = source_parameters<T>::populate_source_parameters(params);
	s_param.sky_average = params->sky_average;
	s_param.f_ref = params->f_ref;
	s_param.phiRef = params->phiRef;
	s_param.cosmology=params->cosmology;
	s_param.incl_angle=params->incl_angle;
	IMRPhenomD<T> model;
	lambda_parameters<T> lambda;
	model.assign_lambda_param(&s_param,&lambda);	
	model.post_merger_variables(&s_param);
	T fRD = s_param.fRD;
	T fdamp = s_param.fdamp;
	T fpeak = model.fpeak(&s_param , &lambda);
	T deltaf = frequencies[1]-frequencies[0];
	//################################################
	//Factor of 2 pi for the definition of time from frequency
	//IMRPhenomD returns (-) phase
	if(local_gen == "IMRPhenomD"){
		//Currently using Nico's fix
		if(correct_time){
			T f = frequencies[0];
			bool check = true, check2=true;
			T pt1, pt2, f1,f2;
			int i = 0 ;
			//One sided, to start it off
			times[0] = (phase_plus[1]-phase_plus[0])/(2.*M_PI*deltaf);
			i++;
			while(f < .95*fRD && i<length-1)
			{
				f = frequencies[i];
				//central difference for the rest of the steps
				times[i] = (phase_plus[i+1]-phase_plus[i-1])/(4.*M_PI*deltaf);
				if(check){
					if(f>fpeak){
						pt1 = times[i];
						f1 = f;
						check=false;
					}
				}
				else{
					if(check2){
						if(f>.9*fRD){
							pt2 = times[i];
							f2 = f;
							check2=false;
						}
					}
				}
				i++;
			}	
			T f_intercept = times[i-1];
			T f_mr = f;
			T freq_slope_pm = (pt2-pt1)/(f2-f1);
			while(f<1.5*fRD && i<length-1)
			{
				f = frequencies[i];
				times[i] =f_intercept+ (f-f_mr)*freq_slope_pm;
				//Stop if observation goes past 20 years
				i++;
				if(times[i-1]>(params->tc+630720000)){
					break;
				}

			}
			T time_transition = times[i-1];
			T f_transition = f;
			while(i<length-1)
			{
				f = frequencies[i];
				times[i] = time_transition +pow_int(f-f_transition,2);
				//Stop if observation goes past 20 years
				i++;
				if(times[i-1]>(params->tc+630720000)){
					break;
				}
				
			}
			//if stopped at 20 years, fill out with frozen time
			if(i != length){
				while(i<length){
					times[i] = times[i-1];
					i++;
				}
			}
			else{
				times[length-1] = (phase_plus[length-1] - phase_plus[length-2])/(2*M_PI*deltaf);
			}
		}
		else{
			times[0] = (phase_plus[1]-phase_plus[0])/(2*M_PI*deltaf);
			if(length>2){
				for(int i = 1  ;i<length-1; i++){
					times[i] = (phase_plus[i+1]-phase_plus[i-1])/(4*M_PI*deltaf);
				}
				times[length-1] = (phase_plus[length-1]-phase_plus[length-2])/(2*M_PI*deltaf);
			}
		}
	}
	//################################################
	delete [] phase_plus;
	delete [] phase_cross;
	//params->shift_time = save_shift_time;
}

/*! \brief map between inclination angle and polarization angle to the spherical coordinates of the binary's total angular momentum in the SS frame
 *
 * CURRENTLY NOT RIGHT MUST FIX
 */
template<class T>
void map_extrinsic_angles(gen_params_base<T> *params)
{
	params->LISA_thetal = params->incl_angle;
	params->LISA_phil = params->psi;
}

void assign_freq_boundaries(double *freq_boundaries, 
	double *intermediate_freqs, 
	int boundary_num, 
	gen_params_base<double> *input_params, 
	std::string generation_method)
{
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
	if(	(
		generation_method =="IMRPhenomPv2" ||
		generation_method =="MCMC_IMRPhenomPv2" ||
		generation_method =="ppE_IMRPhenomPv2_Inspiral" ||
		generation_method =="ppE_IMRPhenomPv2_IMR" ||
		generation_method =="MCMC_ppE_IMRPhenomPv2_Inspiral" ||
		generation_method =="MCMC_ppE_IMRPhenomPv2_IMR" ||
		generation_method =="dCS_IMRPhenomPv2" ||
		generation_method =="EdGB_IMRPhenomPv2" ||
		generation_method =="MCMC_dCS_IMRPhenomPv2" ||
		generation_method =="MCMC_EdGB_IMRPhenomPv2" 
		)
		&& !input_params->sky_average){

		IMRPhenomPv2<adouble> modelp;
		s_param.spin1z = internal_params.spin1[2];
		s_param.spin2z = internal_params.spin2[2];
		s_param.chip = internal_params.chip;
		s_param.phip = internal_params.phip;
		modelp.PhenomPv2_Param_Transform_reduced(&s_param);
		modelp.assign_lambda_param(&s_param, &lambda);
		modelp.post_merger_variables(&s_param);
		double M = s_param.M.value();
		double fRD = s_param.fRD.value();
		double fpeak = modelp.fpeak(&s_param, &lambda).value();
		//###########################################
		freq_boundaries[0] = .014/M;
		freq_boundaries[1] = .018/M;
		if(fRD/2. < fpeak){
			freq_boundaries[2] = fRD/2.;
			freq_boundaries[3] = fpeak;
		}
		else{
			freq_boundaries[3] = fRD/2.;
			freq_boundaries[2] = fpeak;
		}
		freq_boundaries[4] = .2/M;//End waveform
		//###########################################
		intermediate_freqs[0] = freq_boundaries[0]*.9;
		for(int i = 1 ; i<boundary_num; i++){
			intermediate_freqs[i] = freq_boundaries[i-1]+(double)(freq_boundaries[i]-freq_boundaries[i-1])/2.;
		}
	}
	else if(
		(generation_method =="IMRPhenomD" || 
		generation_method == "ppE_IMRPhenomD_Inspiral"|| 
		generation_method == "ppE_IMRPhenomD_IMR") 
		|| 
		(generation_method =="MCMC_IMRPhenomD_Full" || 
		generation_method == "MCMC_ppE_IMRPhenomD_Inspiral_Full"|| 
		generation_method == "MCMC_ppE_IMRPhenomD_IMR_Full")
		|| 
		(generation_method =="MCMC_dCS_IMRPhenomD_Full" || 
		generation_method == "MCMC_EdGB_IMRPhenomD_Full")
		|| 
		(generation_method =="dCS_IMRPhenomD" || 
		generation_method == "EdGB_IMRPhenomD")
		|| 
		(generation_method =="MCMC_IMRPhenomD" || 
		generation_method == "MCMC_ppE_IMRPhenomD_Inspiral" ||
		generation_method == "MCMC_ppE_IMRPhenomD_IMR" )
		){

		IMRPhenomD<adouble> modeld;
		modeld.assign_lambda_param(&s_param, &lambda);
		modeld.post_merger_variables(&s_param);
		double M = s_param.M.value();
		double fRD = s_param.fRD.value();
		double fpeak = modeld.fpeak(&s_param, &lambda).value();
		//###########################################
		freq_boundaries[0] = .014/M;
		freq_boundaries[1] = .018/M;
		if(fRD/2. < fpeak){
			freq_boundaries[2] = fRD/2.;
			freq_boundaries[3] = fpeak;
		}
		else{
			freq_boundaries[3] = fRD/2.;
			freq_boundaries[2] = fpeak;
		}
		freq_boundaries[4] = .2/M;//End waveform
		//###########################################
		intermediate_freqs[0] = freq_boundaries[0]*.9;
		for(int i = 1 ; i<boundary_num; i++){
			intermediate_freqs[i] = freq_boundaries[i-1]+(double)(freq_boundaries[i]-freq_boundaries[i-1])/2.;
		}
	}
	if(check_ppE(generation_method)){
		delete [] internal_params.betappe;
		delete [] internal_params.bppe;
	}
	
	

}
bool check_ppE(std::string generation_method)
{
	if(generation_method.find("ppE") != std::string::npos || 
		generation_method.find("dCS") !=std::string::npos ||
		generation_method.find("EdGB") !=std::string::npos 
		)
	{
		return true;
		
	}
	return false;
}

/*! \brief Utility to find the integration bounds for Fisher matrices for increasing speed of Fisher evaluation
 *
 * Numerically finds the frequencies at which the Fisher should be evaluated at. 
 * 
 * Uses the bisection search algorithm for the cases where the waveform enters/leaves the band at SNR>1
 *
 * Uses a 100 pt grid search (logarithmically spaced) if the signal has SNR<1 when entering and leaving
 *
 * integrand_bounds[0] ~ frequency at which |h|/(sqrt S) ~signal_to_noise +/- tol
 *
 * integrand_bounds[1] ~ frequency at which |h|/(sqrt S) ~signal_to_noise +/- tol
 */
void integration_bounds(gen_params_base<double> *params, /**< Parameters of the waveform*/
	std::string generation_method, /*<< Generation method to use for the waveform*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use (must be one of the analytic curves in the detector_utilitiy file*/
	double fmin, /**< minimum frequency to use (specific to the detector)*/
	double fmax, /**< max frequency to use (specific to the detector)*/
	double signal_to_noise,/**< Target ratio of |h|/ sqrt(S) (typically ~0.1)*/
	double tol,/**< This is a numerical algorithm, so the tolerance must be specified*/
	double *integration_bounds/**< [out] bounds fo the integral shape -- [2] -- (fmin,fmax)*/
	) 
{
	double integration_time = 12;
	bool lower=false, upper=false;
	std::complex<double> response;
	double eval_freq;
	double time;
	double psd;
	double ratio;
	
	//Check lowest frequency
	eval_freq = fmin;	
	time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
		generation_method, true);
	fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
		generation_method, params, &time);
	populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
	ratio = std::abs(response)/psd ;
	if( ratio>(signal_to_noise-tol))
	{
		integration_bounds[0] = eval_freq;
		lower = true;
	}

	//Check highest frequency
	eval_freq = fmax;
	time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
		generation_method, true);
	fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
		generation_method, params, &time);
	populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
	ratio = std::abs(response)/psd ;
	if( ratio>(signal_to_noise-tol))
	{
		integration_bounds[1] = eval_freq;
		upper = true;
	}

	//Search for bounds
	double fmin_search = fmin;	
	double fmax_search = fmax;	
	bool continue_search=true;
	if(upper && !lower)
	{
		while(continue_search)
		{
			//Bisection in log freq space
			eval_freq=sqrt(fmax_search*fmin_search);
			time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
				generation_method, true);
			fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
				generation_method, params, &time);
			populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
			ratio = std::abs(response)/psd ;
			
			//The function can be so steep the algorithm cannot determine a valid 
			//frequency with the required tolerance because of floating point error
			//check to see if the difference is at all meaningful
			if(2*(fmax_search - fmin_search)/(fmax_search+fmin_search) <1e-12){
				continue_search =false;
				integration_bounds[1]=eval_freq;

			}
			if(ratio > (signal_to_noise + tol)){
				fmax_search = eval_freq;
			}
			else if(ratio < (signal_to_noise - tol)){
				fmin_search = eval_freq;
			}
			else{
				continue_search =false;
				integration_bounds[0]=eval_freq;
			}
		}
	}
	else if(lower && !upper){
		while(continue_search)
		{
			//Bisection in log freq space
			eval_freq=sqrt(fmax_search*fmin_search);
			time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
				generation_method, true);
			fourier_detector_response_equatorial(&eval_freq, 1, &response, detector, 
				generation_method, params, &time);
			populate_noise(&eval_freq, sensitivity_curve, &psd, 1,integration_time);
			ratio = std::abs(response)/psd ;

			//The function can be so steep the algorithm cannot determine a valid 
			//frequency with the required tolerance because of floating point error
			//check to see if the difference is at all meaningful
			if(2*(fmax_search - fmin_search)/(fmax_search+fmin_search) <1e-12){
				continue_search =false;
				integration_bounds[1]=eval_freq;

			}
			if(ratio > (signal_to_noise + tol)){
				fmin_search = eval_freq;
			}
			else if(ratio < (signal_to_noise - tol)){
				fmax_search = eval_freq;
			}
			else{
				continue_search =false;
				integration_bounds[1]=eval_freq;
			}
		}

	}
	//For now, just doing a sparse grid search. This can be refined later 
	//After doing grid search, it would be best to do a bisection search like above
	//with the new bounds after the grid search
	else if(!lower && !upper){
		//Grid search to see if any point is above the noise curve
		int vec_length= 100;
		std::complex<double> response_vec[vec_length];
		double time_vec[vec_length];
		double freq_vec[vec_length];
		double psd_vec[vec_length];
		double ratio_vec[vec_length];

		//Logarithmically space frequencies
		double delta_f_factor = pow(fmax/fmin,1./(vec_length-1));
		for(int i = 0 ; i<vec_length; i++){
			freq_vec[i]= pow_int(delta_f_factor,i)*fmin;
		}
		time_phase_corrected_autodiff(time_vec, vec_length, freq_vec, params, 
			generation_method, true);
		fourier_detector_response_equatorial(freq_vec, vec_length, response_vec, detector, 
			generation_method, params, time_vec);
		populate_noise(freq_vec, sensitivity_curve, psd_vec, vec_length,integration_time);
		for(int i = 0 ; i<vec_length; i++){
			ratio_vec[i]=std::abs(response_vec[i])/psd_vec[i];
		}
		bool search_lower = true, search_upper = true;
		for(int i = 0 ; i<vec_length; i++){
			if(search_lower)
			{
				if(ratio_vec[i] >(signal_to_noise -tol) )
				{
					fmin_search = freq_vec[i];
					search_lower=false;
				}		
			}
			else if(search_upper){
				if(ratio_vec[i] <(signal_to_noise -tol) )
				{
					fmax_search = freq_vec[i];
					search_upper=false;
					break;
				}
			}
		}
		if(!search_lower && !search_upper){
			integration_bounds[0]=fmin_search;
			integration_bounds[1]=fmax_search;
		}
		else{
			integration_bounds[0]=fmin;
			integration_bounds[1]=fmax;
		}
	}
}

//###########################################################################
template void map_extrinsic_angles<double>(gen_params_base<double> *);
template void map_extrinsic_angles<adouble>(gen_params_base<adouble> *);

template void  time_phase_corrected<double>(double *, int, double *, gen_params_base<double> *, std::string, bool );
template void  time_phase_corrected<adouble>(adouble *, int, adouble *, gen_params_base<adouble> *, std::string, bool);

template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::complex<double> *,std::complex<double> *, double, double, std::string);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *,std::complex<adouble> *, adouble, adouble, std::string);
//
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double, double, std::string);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble, adouble, std::string);
//
//
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *);
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *, double *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *, adouble *);
//
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double>*);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble>*);
//
template int fourier_detector_response_equatorial<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double , double, double,double *, double, double, double, double, std::string);
template int fourier_detector_response_equatorial<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble , adouble, double,adouble*, adouble, adouble, adouble, adouble,  std::string);
