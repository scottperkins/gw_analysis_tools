#include "waveform_util.h"
#include "util.h"
#include "GWATConfig.h"
#include "waveform_generator.h"
#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomP.h"
#include "detector_util.h"
#include "pn_waveform_util.h"
#include <fftw3.h>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <adolc/taping.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <gsl/gsl_integration.h>
/*!\file 
 * Utilities for waveforms - SNR calculation and detector response
 * 	
 * includes snr and detector response
 *
 * Includes some utilities useful for MCMC and fisher calculations, as well as time-frequency methods for detectors like LISA
 */

struct gsl_snr_struct
{
	gen_params_base<double> *params;
	std::string SN;
	std::string generation_method;
	std::string detector;
};

double data_snr(double *frequencies, 
	int length,
	std::complex<double> *data,
	std::complex<double> *response,
	double *psd
	)	
{
	double *integrand
			 = (double *)malloc(sizeof(double)*length);

	double delta_f = frequencies[1]-frequencies[0];
	//for (int i = 0; i<length;i++)
	//	integrand[i] = 4.*real(conj(response[i])*response[i])/psd[i]; 
	//double snr_template;
	//snr_template = sqrt(simpsons_sum(delta_f,length, integrand));

	
	for (int i = 0; i<length;i++)
		integrand[i] = 4.*real(conj(data[i])*response[i])/psd[i]; 
	//double inner_prod = sqrt(simpsons_sum(delta_f,length, integrand));
	double inner_prod = (simpsons_sum(delta_f,length, integrand));
	//std::cout<<"WU: "<<inner_prod<<std::endl;
	free(integrand);
	//return inner_prod/(snr_template);
	return sqrt(inner_prod);

}

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

double calculate_snr(std::string sensitivity_curve,
	std::string detector,
	std::string generation_method,
	gen_params_base<double> *params,
	double *frequencies,
	int length,
	std::string integration_method,
	double *weights,
	bool log10_freq)
{
	double *times;
	if(params->equatorial_orientation){
		transform_orientation_coords(params, generation_method,detector);
	}
	if(detector == "LISA" && !params->sky_average){
		times = new double[length];
		if(integration_method == "GAUSSLEG"){
			time_phase_corrected_autodiff(times, length, frequencies, params, generation_method, false, NULL);
		}
		else{
			time_phase_corrected(times, length, frequencies, params, generation_method, false);
		}
	}
	std::complex<double> *response = new std::complex<double>[length];
	fourier_detector_response(frequencies, length, response, detector, generation_method, params,times);
	double snr = calculate_snr(sensitivity_curve, response, frequencies, length);
	if(detector == "LISA"){
		//snr+=calculate_snr(sensitivity_curve, response, frequencies, length);
		if(!params->sky_average){
			delete [] times;
		}
		//snr*=sqrt(2); //Two detectors
	}
	delete [] response;	
	return snr;

}
/**< \brief Routine to calculate the SNR of a template with GSL quadrature integration
 *
 * Sometimes, this is faster than the ``grid'' style integration
 *
 * Supports sky-averaged templates, but this should only be used with non-precessing waveforms
 *
 */
int calculate_snr_gsl(double *snr,/**< [out] SNR*/
	std::string sensitivity_curve,/**< Noise curve */
	std::string detector,/**<Detector to compute response -- can be empty is SA*/
	std::string generation_method,/**<Generation method */
	gen_params_base<double> *params,/**< Source Parameters*/
	double f_min,/**< Lower frequency bound*/
	double f_max,/**< Upper frequency bound*/
	double relative_error/**< Relative error threshold*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params, generation_method,detector);
	}
	int np=100;
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(np) ;
	int err =  calculate_snr_gsl(snr,sensitivity_curve, detector, generation_method, params, f_min, f_max,relative_error,w,np);
	gsl_integration_workspace_free(w);
	return err;
}
/**< \brief Routine to calculate the SNR of a template with GSL quadrature integration
 *
 * Sometimes, this is faster than the ``grid'' style integration
 *
 * Supports sky-averaged templates, but this should only be used with non-precessing waveforms
 *
 */
int calculate_snr_gsl(double *snr,
	std::string sensitivity_curve,/**< Noise curve */
	std::string detector,/**<Detector to compute response -- can be empty is SA*/
	std::string generation_method,/**<Generation method */
	gen_params_base<double> *params,/**< Source Parameters*/
	double f_min,/**< Lower frequency bound*/
	double f_max,/**< Upper frequency bound*/
	double relative_error,/**< Relative error threshold*/
	gsl_integration_workspace *w, /**< User-allocated gsl_integration_workspace*/
	int np/**<Size of gsl_integration_workspace allocation*/
	)
{
	bool SA_save=params->sky_average;
	if(sensitivity_curve.find("SADC") != std::string::npos && params->sky_average){
		params->sky_average=false;	
	}
	gsl_snr_struct helper_params;
	helper_params.params = params;
	helper_params.generation_method = generation_method;
	helper_params.SN = sensitivity_curve;
	helper_params.detector = detector;
	gsl_function F;
	if(SA_save){
		F.function = [](double f, void *param){return integrand_snr_SA_subroutine(f,param);};
	}
	else{
		F.function = [](double f, void *param){return integrand_snr_subroutine(f,param);};
	}
	F.params = (void *)&helper_params;
	double result, err;
	gsl_set_error_handler_off();
	int errcode = gsl_integration_qag(&F,f_min, f_max, 0,relative_error, np, GSL_INTEG_GAUSS15,w, &result, &err);
	*snr = sqrt(result);
	//sqrt 2 for second LISA detector
	//if(detector=="LISA"){
	//	*snr *=sqrt(2.);	
	//}
	params->sky_average = SA_save;
	return errcode;
}

/*! \brief Internal function to calculate the SNR integrand for sky-averaged waveforms
 */
double integrand_snr_SA_subroutine(double f, void *subroutine_params)
{
	gsl_snr_struct cast_params = *(gsl_snr_struct *)subroutine_params;
	std::complex<double> wfp, wfc;
	fourier_waveform(&f, 1,&wfp, &wfc, cast_params.generation_method, cast_params.params);
	double SN;
	populate_noise(&f, cast_params.SN, &SN,1);
	SN*=SN;
	//std::cout<<f<<" "<<wfp<<" "<<SN<<std::endl;
	return 4*std::real(std::conj(wfp)*wfp)/SN;
}
/*! \brief Internal function to calculate the SNR integrand for full waveforms
 */
double integrand_snr_subroutine(double f, void *subroutine_params)
{
	gsl_snr_struct cast_params = *(gsl_snr_struct *)subroutine_params;
	double times[2];
	if(cast_params.detector == "LISA" ){
		//PROBLEM
		double temp_f[2] = {f, f+1.e-5};
		//times = new double[length];
		//time_phase_corrected_autodiff(&times[0], 1, &f, cast_params.params, cast_params.generation_method, false, NULL);
		time_phase_corrected(times, 2, temp_f, cast_params.params, cast_params.generation_method, false);
	}
	double time = times[0];

	std::complex<double> response;
	fourier_detector_response(&f, 1,&response, cast_params.detector,cast_params.generation_method, cast_params.params, &time);
	double SN;
	populate_noise(&f, cast_params.SN, &SN,1);
	SN*=SN;
	//std::cout<<4*std::real(std::conj(response)*response)/SN<<std::endl;
	//std::cout<<SN<<" "<<std::real(std::conj(response)*response)<<std::endl;
	return 4*std::real(std::conj(response)*response)/SN;
}
/*! \brief Caclulates the snr given a detector and waveform (complex) and frequencies
 *      
 * This function computes the un-normalized snr: \sqrt( ( H | H ) )
 */     
double calculate_snr(std::string sensitivity_curve, /**< detector name - must match the string of populate_noise precisely*/
                        std::complex<double> *waveform,/**< complex waveform */
                        double *frequencies,/**< double array of frequencies that the waveform is evaluated at*/
                        int length/**< length of the above two arrays*/
                        )
{
        double *noise = (double *)malloc(sizeof(double)*length);
        populate_noise(frequencies,sensitivity_curve, noise,  length);
        for (int i = 0; i< length; i++){
                noise[i] = noise[i]*noise[i];
	}
        double *integrand = (double *) malloc(sizeof(double)*length);
        for (int i = 0; i<length; i++)
                integrand[i] = 4.* real(conj(waveform[i])*waveform[i]/noise[i]);
        double integral = trapezoidal_sum(frequencies, length, integrand);
	free(integrand);
	free(noise);
        return sqrt(integral);
}

double calculate_snr_internal(double *psd, 
	std::complex<double>*waveform,
	double *frequencies, 
	int length,
	std::string integration_method,
	double *weights,
	bool log10_freq)
{
        double *integrand = (double *) malloc(sizeof(double)*length);
        for (int i = 0; i<length; i++){
                integrand[i] = 4.* real(conj(waveform[i])*waveform[i]/psd[i]);
		if(log10_freq){
			integrand[i]*= (frequencies[i]*LOG10);	
		}
	}
        double integral=0;
	if(integration_method == "SIMPSONS"){
		//integral = trapezoidal_sum(frequencies, length, integrand);
		integral = simpsons_sum(1./(frequencies[1]-frequencies[0]), length, integrand);
	}
	else if(integration_method == "GAUSSLEG"){
        	for (int i = 0; i<length; i++){
			integral+= weights[i]*integrand[i];
		
		}
	
	}
	else{
		std::cout<<"UNSUPPORTED INTEGRATION METHOD"<<std::endl;	
	}
	free(integrand);
        return sqrt(integral);
}

/* \brief calculates the detector response for a given waveform and detector -- polarization angle =0
 */
template<class T>
int fourier_detector_response_horizon(T *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			std::complex<T> *hplus, /*<precomputed plus polarization of the waveform*/ 
			std::complex<T> *hcross, /**<precomputed cross polarization of the waveform*/ 
			std::complex<T> *detector_response, /**< [out] detector response*/
			T theta, /**< polar angle (rad) theta in detector frame*/
			T phi, /**< azimuthal angle (rad) phi in detector frame*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	return fourier_detector_response_horizon(frequencies, length, hplus, hcross, detector_response, theta, phi, (T)0., detector);
	
}
/* \brief calculates the detector response for a given waveform and detector
 */
template<class T>
int fourier_detector_response_horizon(T *frequencies, /**<array of frequencies corresponding to waveform*/
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
		Fcross = fplus*s2psi+ fcross*c2psi;
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
int fourier_detector_response_horizon(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
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
	status = fourier_detector_response_horizon(frequencies, 
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
			T theta_j_ecl,
			T phi_j_ecl,
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
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,theta_j_ecl,phi_j_ecl, Fplus, Fcross);
	}
	else{
		detector_response_functions_equatorial(detector, ra, dec, psi, gmst,times, length,LISA_alpha0, LISA_phi0,theta_j_ecl,phi_j_ecl, &fplus, &fcross);
	}
	
	if(detector == "LISA"){
		for (int i =0; i <length; i++)
		{
			//Doppler phase shift
			T theta_s, phi_s, phi_t;
			ecl_from_eq((T)(M_PI/2. - dec), ra, &theta_s, &phi_s);
			phi_t = LISA_phi0 + 2. * M_PI * times[i] / T_year;
			T doppler_phase_shift = 2*M_PI * frequencies[i] * AU_SEC *sin(theta_s) *cos(phi_t - phi_s);

			detector_response[i] = (Fplus[i] * hplus[i] 
						+ (Fcross[i] )*hcross[i]) * std::exp(std::complex<T>(0, - doppler_phase_shift ) );
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
	if(parameters->equatorial_orientation){
		transform_orientation_coords(parameters, generation_method,detector);
	}
	else{
		if(detector=="LISA"){
			std::cout<<"ERROR -- fourier_detector_response_equatorial -- LISA currently only accepts equatorial_orientation parameters"<<std::endl;
		}
	}
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
			parameters->theta_j_ecl,
			parameters->phi_j_ecl,
			detector 
			) ;
	
		
	//Deallocate memory
	delete [] waveform_plus;
	delete [] waveform_cross;

	return status;
}

/*! \brief Wrapper to handle all detector_response calls -- horizon and equatorial
 *
 */
template<class T>
int fourier_detector_response(T *frequencies,
	int length,
	std::complex<T> *response,
	std::string detector,
	std::string generation_method,
	gen_params_base<T> *parameters,
	T *times)
{
	int status;
	if(parameters->horizon_coord){
		status=fourier_detector_response_horizon(frequencies, length, response, detector, generation_method, parameters);
	}
	else{
		status=fourier_detector_response_equatorial(frequencies, length, response, detector, generation_method, parameters, times);

	}
	return status;

}
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double> *, double *);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble> *, adouble *);

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
	fourier_detector_response_horizon(frequencies, length,
				response, detector, generation_method, parameters);
	for(int i =0 ; i< length; i++){
		amplitude[i] = std::abs(response[i]);
		phase[i] = std::arg(response[i]);
	}

	free(response);
	return 0;
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
void time_phase_corrected_autodiff(double *times, int length, double *frequencies,gen_params_base<double> *params, std::string generation_method, bool correct_time, int *tapes_in)
{
	bool save_dep = params->dep_postmerger;
	params->dep_postmerger=true;
	bool save_shift_time = params->shift_time;
	params->shift_time=false;
	//if(generation_method.find("Pv2") != std::string::npos && params->chip == -1){
	//	IMRPhenomPv2<double> model;
	//	params->chip = model.PhenomPv2_inplane_spin(params);
	//	params->phip = 2*M_PI;
	//}

	int boundary_num = boundary_number(generation_method);
	double freq_boundaries[boundary_num];
	double grad_freqs[boundary_num];
	assign_freq_boundaries(freq_boundaries, grad_freqs, boundary_num, params, generation_method);	
	int *tapes;
	gen_params_base<adouble> aparams;
	if(tapes_in){ 
		tapes = tapes_in;
	}
	else{
		tapes = new int[boundary_num];
		transform_parameters(params, &aparams);
		if(params->equatorial_orientation){
			transform_orientation_coords(&aparams,generation_method,"");
		}
		for(int i = 0 ; i<boundary_num ; i++){
			tapes[i]=(i+1)*8;	
			trace_on(tapes[i]);
			adouble freq;
			freq <<= grad_freqs[i];
			adouble phasep, phasec;
			fourier_phase(&freq, 1, &phasep, &phasec, generation_method, &aparams);
			double phaseout;
			phasep>>=phaseout;
			trace_off();
		}
	}
	bool eval = false;
	double freq;	
	for(int k = 0; k<length; k++){
		freq = frequencies[k];
		for(int n = 0; n<boundary_num ; n++){
			if(freq < freq_boundaries[n]){
				//std::cout<<freq<<std::endl;
				//std::cout<<n<<std::endl;
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
	if(!tapes_in){
		delete [] tapes;
		if(check_mod(generation_method)){
			delete [] aparams.betappe;
			delete [] aparams.bppe;
		}
	}
	params->dep_postmerger=save_dep;
	params->shift_time=save_shift_time;
	
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
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	std::string local_gen = generation_method;
	if(length ==1){
		T temp_f[2];
		T epsilon =1e-6;
		temp_f[0] = frequencies[0]-epsilon;
		temp_f[1] = frequencies[0]+epsilon;
		T temp_deltaf = 2*epsilon;
		T temp_phase_plus[2];
		T temp_phase_cross[2];
		fourier_phase(temp_f, 2, temp_phase_plus, temp_phase_cross, local_gen, params);
		times[0] = (temp_phase_plus[1]-temp_phase_plus[0])/(4*M_PI*temp_deltaf);
		return ;
	}
	//bool save_shift_time = params->shift_time;
	//params->shift_time = false;
	//std::string local_gen = "IMRPhenomD";
	//################################################
	T *phase_plus = new T[length];
	T *phase_cross = new T[length];
	//################################################
	//################################################
	fourier_phase(frequencies, length, phase_plus, phase_cross, local_gen, params);
	//################################################
	T fdamp, fRD, fpeak, deltaf;
	if(local_gen.find("IMRPhenomD")!=std::string::npos){
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
		fRD = s_param.fRD;
		fdamp = s_param.fdamp;
		fpeak = model.fpeak(&s_param , &lambda);
		deltaf = frequencies[1]-frequencies[0];
	}
	else if(local_gen.find("IMRPhenomPv2")!=std::string::npos){
		source_parameters<T> s_param;
		s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		s_param.chip = params->chip;
		s_param.phip = params->phip;
		IMRPhenomPv2<T> model;
		lambda_parameters<T> lambda;
		model.assign_lambda_param(&s_param,&lambda);	
		model.post_merger_variables(&s_param);
		fRD = s_param.fRD;
		fdamp = s_param.fdamp;
		fpeak = model.fpeak(&s_param , &lambda);
		deltaf = frequencies[1]-frequencies[0];
	}
	//################################################
	//Factor of 2 pi for the definition of time from frequency
	//IMRPhenomD returns (-) phase
	if(local_gen.find("IMRPhenom") !=std::string::npos){
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
			else{
				//only option with length 2 vector
				times[1] = times[0];
			}
		}
	}
	//################################################
	delete [] phase_plus;
	delete [] phase_cross;
	//params->shift_time = save_shift_time;
}

template<class T>
void transform_orientation_coords(gen_params_base<T> *parameters,std::string generation_method,std::string detector)
{
	if(generation_method.find("IMRPhenomP") != std::string::npos){
		T theta_s = M_PI/2. - parameters->DEC;
		T phi_s = parameters->RA;
		//N should point source to detector
		T Neq[3] = {-sin(theta_s)*cos(phi_s), -sin(theta_s)*sin(phi_s), -cos(theta_s)};
		T Leq[3] = {sin(parameters->theta_l)*cos(parameters->phi_l), 
			sin(parameters->theta_l)*sin(parameters->phi_l), 
			cos(parameters->theta_l)};
		parameters->incl_angle = acos(Neq[0]*Leq[0]+Neq[1]*Leq[1]+Neq[2]*Leq[2]);
		if(parameters->incl_angle == 0 || parameters->incl_angle == M_PI){
			std::cout<<"ERROR -- L == N || L == -N : The source frame is ill defined, and the conversion between theta_l and phi_l to theta_j and phi_j should not be trusted"<<std::endl; 
		}
		//Populate JSF here
		T JSF[3];
		if(generation_method.find("Pv2")!=std::string::npos){
			IMRPhenomPv2<T> model;
			model.PhenomPv2_JSF_from_params(parameters, JSF);
		}
		else{
			std::cout<<"ERROR -- Only Pv2 is supported right now"<<std::endl;	
		}
		T Jeq[3];
		equatorial_from_SF(JSF, parameters->theta_l, parameters->phi_l,(T) (theta_s), phi_s, parameters->incl_angle, parameters->phiRef, Jeq);
		T Jeqsph[3];
		transform_cart_sph(Jeq,Jeqsph);
		T theta_j = Jeqsph[1];
		T phi_j= Jeqsph[2];

		T iota_j;
		if(detector==""){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC, theta_j, phi_j,&parameters->psi, &iota_j);
			ecl_from_eq(theta_j, phi_j, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
		else if(detector!="LISA"){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC, theta_j, phi_j,&parameters->psi, &iota_j);
		}
		else{
			ecl_from_eq(theta_j, phi_j, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
	}
	else{
		T theta_s = M_PI/2. - parameters->DEC;
		T phi_s = parameters->RA;
		//N should point source to detector
		T Neq[3] = {-sin(theta_s)*cos(phi_s), -sin(theta_s)*sin(phi_s), -cos(theta_s)};
		T Leq[3] = {sin(parameters->theta_l)*cos(parameters->phi_l), 
			sin(parameters->theta_l)*sin(parameters->phi_l), 
			cos(parameters->theta_l)};
		parameters->incl_angle = acos(Neq[0]*Leq[0]+Neq[1]*Leq[1]+Neq[2]*Leq[2]);
		if(detector==""){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC,parameters->theta_l, parameters->phi_l, &parameters->psi, &parameters->incl_angle);
			ecl_from_eq(parameters->theta_l, parameters->phi_l, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
		if(detector!="LISA"){
			terr_pol_iota_from_equat_sph(parameters->RA, parameters->DEC,parameters->theta_l, parameters->phi_l, &parameters->psi, &parameters->incl_angle);
		}
		else{
			ecl_from_eq(parameters->theta_l, parameters->phi_l, &parameters->theta_j_ecl, &parameters->phi_j_ecl);
		}
	
	}
}
template void transform_orientation_coords<double>(gen_params_base<double> *, std::string,std::string);
template void transform_orientation_coords<adouble>(gen_params_base<adouble> *, std::string,std::string);

void assign_freq_boundaries(double *freq_boundaries, 
	double *intermediate_freqs, 
	int boundary_num, 
	gen_params_base<double> *input_params, 
	std::string generation_method)
{
	if(input_params->equatorial_orientation){
		transform_orientation_coords(input_params,generation_method,"");
	}
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
	//IMRPhenomPv2<double> model;
	//model.PhenomPv2_inplane_spin(input_params);
	if(	(
		generation_method.find("IMRPhenomPv2") != std::string::npos
		)
		&& !input_params->sky_average){

		IMRPhenomPv2<adouble> modelp;
		s_param.spin1z = internal_params.spin1[2];
		s_param.spin2z = internal_params.spin2[2];
		if((internal_params.chip + 1)>DOUBLE_COMP_THRESH){
			s_param.chip = internal_params.chip;
			s_param.phip = internal_params.phip;
			modelp.PhenomPv2_Param_Transform_reduced(&s_param);
		}
		else{
			s_param.spin1x = internal_params.spin1[0];
			s_param.spin1y = internal_params.spin1[1];
			s_param.spin1z = internal_params.spin1[2];
			s_param.spin2x = internal_params.spin2[0];
			s_param.spin2y = internal_params.spin2[1];
			s_param.spin2z = internal_params.spin2[2];
			modelp.PhenomPv2_Param_Transform(&s_param);

		}
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
		generation_method.find("IMRPhenomD") != std::string::npos
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
	if(check_mod(generation_method)){
		delete [] internal_params.betappe;
		delete [] internal_params.bppe;
	}
	
	

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
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	double integration_time = 48;//Just for sensitivity curve
	bool lower=false, upper=false;
	std::complex<double> response;
	double eval_freq;
	double time;
	double psd;
	double ratio;
	
	//Check lowest frequency
	eval_freq = fmin;	
	//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
	//	generation_method, true);
	time_phase_corrected(&time, 1, &eval_freq, params, 
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
	//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
	//	generation_method, true);
	time_phase_corrected(&time, 1, &eval_freq, params, 
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
			//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			//	generation_method, true);
			time_phase_corrected(&time, 1, &eval_freq, params, 
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
			//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			//	generation_method, true);
			time_phase_corrected(&time, 1, &eval_freq, params, 
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
		//time_phase_corrected_autodiff(time_vec, vec_length, freq_vec, params, 
		//	generation_method, true);
		time_phase_corrected(time_vec, vec_length, freq_vec, params, 
			generation_method, true);
		fourier_detector_response_equatorial(freq_vec, vec_length, response_vec, detector, 
			generation_method, params, time_vec);
		populate_noise(freq_vec, sensitivity_curve, psd_vec, vec_length,integration_time);
		for(int i = 0 ; i<vec_length; i++){
			ratio_vec[i]=std::abs(response_vec[i])/psd_vec[i];
		}
		bool search_lower = true, search_upper = true;
		for(int i = 0 ;i<vec_length; i++){
			if(ratio_vec[i] >(signal_to_noise -tol) )
			{
				fmin_search = freq_vec[i];
				search_lower=false;
				break;
			}		

		}
		for(int i = vec_length-1 ;i>0; i--){
			if(ratio_vec[i] >(signal_to_noise -tol) )
			{
				fmax_search = freq_vec[i];
				search_upper=false;
				break;
			}		

		}
		//for(int i = 0 ; i<vec_length; i++){
		//	std::cout<<i<<" "<<ratio_vec[i]<<" "<<freq_vec[i]<<std::endl;
		//	if(search_lower)
		//	{
		//		if(ratio_vec[i] >(signal_to_noise -tol) )
		//		{
		//			fmin_search = freq_vec[i];
		//			search_lower=false;
		//		}		
		//	}
		//	else if(search_upper){
		//		if(ratio_vec[i] <(signal_to_noise -tol) )
		//		{
		//			fmax_search = freq_vec[i];
		//			search_upper=false;
		//			break;
		//		}
		//	}
		//}
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
/*! \brief Determines the integration bounds for the log likelihood or fisher given some observation time, sampling frequency, detector, and sensitivity curve
 *
 * Sensitivity curve has to be one of the options in detector_util analytic options
 *
 * The current scheme is to use the frequency bounds determined by the SNR if the binary spends less than the integration time in band. If the merger spends more time in band than the integration time, the frequencies are determined to be (f_integration_time, f_high_band)
 *
 * Accounts for the integration time, unlike the integration_bounds routine
 *
 * returns 0 if successful, returns 1 if bounds could not be found due to roundoff, and returns 2 if entirely unsuccessful
 */
int integration_interval(double sampling_freq, /**< Frequency at which the detector operates*/
	double integration_time, /**< Time of observation in seconds*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use -- must match analytic choices in detector_util*/
	std::string generation_method,/**< method to use for the waveform generation*/
	gen_params_base<double> *params,/**< parameters of the source*/
	double *freq_bounds/**< [out] Output bounds*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,detector);
	}
	double fmax = sampling_freq /2.; //Nyquist limit
	//double fmin= 0; //DC component
	double fmin= 1e-6; //DC component
	double delta_f =  1./integration_time;
	
	double bounds_from_band[2];
	integration_bounds( params, generation_method, detector, sensitivity_curve, fmin, fmax, .1, .01, bounds_from_band);
	//std::cout<<"Integration bounds from band "<<bounds_from_band[0]<<" "<<bounds_from_band[1]<<std::endl;
	double times[2];
	//time_phase_corrected_autodiff(times, 2, bounds_from_band, params, generation_method, true);
	time_phase_corrected(&times[0], 1, &bounds_from_band[0], params, generation_method, false);
	time_phase_corrected(&times[1], 1, &bounds_from_band[1], params, generation_method, false);
	double T_band = -times[1]+times[0];
	//std::cout<<T_band<<std::endl;

	//Conform the output of the band calculation to the pre-set frequency grid
	bool max_found=false, min_found=false;

	if(T_band < integration_time){
		freq_bounds[0]=bounds_from_band[0];	
		freq_bounds[1]=bounds_from_band[1];	
		return 0;
	}
	else{
		freq_bounds[1] = bounds_from_band[1];
		bool continue_search=true;
		double eval_freq, time;
		double max_search=bounds_from_band[1], min_search=bounds_from_band[0];
		double tolerance = .01*integration_time;
		
		while(continue_search)
		{
			//eval_id = (max_id_search +min_id_search)/2;
			eval_freq = (max_search +min_search)/2;
			//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			//	generation_method, true);
			time_phase_corrected(&time, 1, &eval_freq, params, 
				generation_method, false);
			//std::cout<<(-times[1]+time)/T_year<<std::endl;
			if( 	( ( (-times[1] + time) > integration_time-tolerance )  && 
				( (-times[1] + time) < integration_time+tolerance) ) ) 
			{
				continue_search = false;
				freq_bounds[0]=eval_freq;
				return 0;
				//std::cout<<(-times[1]+time)/T_year<<std::endl;
			}
			else if( fabs(max_search-min_search)<DOUBLE_COMP_THRESH){
				continue_search = false;
				freq_bounds[0] = eval_freq;
				return 1;

			}
			else if(( (-times[1] + time) < (integration_time-tolerance) )){
				max_search = eval_freq;
			}
			else if(( (-times[1] + time) > (integration_time+tolerance) )){
				min_search = eval_freq;
			}

		}
	}
	return 2;
	
}

/*! \brief Determines the integration bounds for the log likelihood or fisher given some observation time, sampling frequency, detector, and sensitivity curve
 *
 * Sensitivity curve has to be one of the options in detector_util analytic options
 *
 * The current scheme is to use the frequency bounds determined by the SNR if the binary spends less than the integration time in band. If the merger spends more time in band than the integration time, the frequencies are determined to be (f_integration_time, f_high_band)
 */
void integration_interval_discrete(double sampling_freq, /**< Frequency at which the detector operates*/
	double integration_time, /**< Time of observation in seconds*/
	std::string detector, /**< Detector to use for the response function*/
	std::string sensitivity_curve, /**< Sensitivity curve to use -- must match analytic choices in detector_util*/
	std::string generation_method,/**< method to use for the waveform generation*/
	gen_params_base<double> *params,/**< parameters of the source*/
	double *freq_bounds/**< [out] Output bounds*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,detector);
	}
	double fmax = sampling_freq /2.; //Nyquist limit
	//double fmin= 0; //DC component
	double fmin= 1e-6; //DC component
	double delta_f =  1./integration_time;
	int N = (fmax-fmin)/(delta_f)+1;
	double *frequencies = new double[N];
	for(int i = 0 ; i<N; i++){
		frequencies[i] = fmin + i*delta_f;	
	}
	
	double bounds_from_band[2];
	integration_bounds( params, generation_method, detector, sensitivity_curve, fmin, fmax, .1, .01, bounds_from_band);
	//std::cout<<"Integration bounds from band "<<bounds_from_band[0]<<" "<<bounds_from_band[1]<<std::endl;
	double times[2];
	//time_phase_corrected_autodiff(times, 2, bounds_from_band, params, generation_method, true);
	time_phase_corrected(times, 2, bounds_from_band, params, generation_method, true);
	double T_band = -times[1]+times[0];
	//std::cout<<T_band<<std::endl;

	//Conform the output of the band calculation to the pre-set frequency grid
	bool max_found=false, min_found=false;
	int max_id , min_id ;
	int i=0 ;

	while((!max_found || !min_found) && i<N){
		if(bounds_from_band[1] >= (frequencies[i]-delta_f/2. ) &&
			bounds_from_band[1] <= (frequencies[i]+delta_f/2. )){
			max_id = i;
			max_found = true;
		}
		if(bounds_from_band[0] >= (frequencies[i]-delta_f/2. ) &&
			bounds_from_band[0] <= (frequencies[i]+delta_f/2. )){
			min_id = i;
			min_found = true;
		}
		i++;
	}

	if(T_band < integration_time){
		freq_bounds[0]=frequencies[min_id];	
		freq_bounds[1]=frequencies[max_id];	
	}
	else{
		freq_bounds[1] = frequencies[max_id];
		bool continue_search=true;
		double eval_freq, time;
		int min_id_search=min_id, max_id_search=max_id, eval_id;
		double tolerance = .1*integration_time;
		
		while(continue_search)
		{
			eval_id = (max_id_search +min_id_search)/2;
			eval_freq = frequencies[eval_id];
			//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
			//	generation_method, true);
			time_phase_corrected(&time, 1, &eval_freq, params, 
				generation_method, true);
			if( 	( ( (-times[1] + time) > integration_time-tolerance )  && 
				( (-times[1] + time) < integration_time+tolerance) )
				||
				(freq_bounds[1]-eval_freq < 100*delta_f) ) 
			{

				continue_search = false;
				freq_bounds[0]=eval_freq;
				//std::cout<<(times[1]-time)/T_year<<std::endl;
			}
			else if(( (-times[1] + time) < (integration_time-tolerance) )){
				max_id_search = eval_id;
			}
			else if(( (-times[1] + time) > (integration_time+tolerance) )){
				min_id_search = eval_id;
			}

		}
	}
	
	delete [] frequencies;
}
template<class T>
void postmerger_params(gen_params_base<T>*params,
	std::string generation_method,
	T *fpeak,
	T *fdamp,
	T *fRD
	)
{
	if(generation_method.find("IMRPhenomD")!=std::string::npos){
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
		*fRD = s_param.fRD;
		*fdamp = s_param.fdamp;
		*fpeak = model.fpeak(&s_param , &lambda);
	}
	else if(generation_method.find("IMRPhenomPv2")!=std::string::npos){
		source_parameters<T> s_param;
		s_param = source_parameters<T>::populate_source_parameters(params);
		s_param.sky_average = params->sky_average;
		s_param.f_ref = params->f_ref;
		s_param.phiRef = params->phiRef;
		s_param.cosmology=params->cosmology;
		s_param.incl_angle=params->incl_angle;
		s_param.chip = params->chip;
		s_param.phip = params->phip;
		IMRPhenomPv2<T> model;
		lambda_parameters<T> lambda;
		model.assign_lambda_param(&s_param,&lambda);	
		model.post_merger_variables(&s_param);
		*fRD = s_param.fRD;
		*fdamp = s_param.fdamp;
		*fpeak = model.fpeak(&s_param , &lambda);
	}

}
template void postmerger_params<double>(gen_params_base<double> *,std::string, double *, double *, double*);
template void postmerger_params<adouble>(gen_params_base<adouble> *,std::string, adouble *, adouble *, adouble*);

/*! \brief Convenience function to Calculate the time before merger using numerical methods
 *
 * Uses numerical -- omp safe and thread safe
 */
void Tbm_to_freq(gen_params_base<double> *params,/**< Generation parameters of the source*/
	std::string generation_method,/**< Generation method for the waveform*/
	double Tbm,/**< target time before merger*/
	double *freq,/**< Frequency at the input time before merger*/
	double tol /**< Tolerance for the scheme*/
	)
{
	if(params->equatorial_orientation){
		transform_orientation_coords(params,generation_method,"");
	}
	double fpeak,fRD,fdamp ;
	postmerger_params(params, generation_method, &fpeak, &fdamp, &fRD);
	//std::cout<<"f fpeak "<<fpeak<<std::endl;	
	double time_peak=1, time_Tbm;
	time_phase_corrected_autodiff(&time_peak, 1, &fpeak, params, 
		generation_method, false);
	//std::cout<<"Time fpeak "<<time_peak<<std::endl;	
	bool continue_search = true;
	double Tbmp = (1+tol)*Tbm;
	double Tbmm = (1-tol)*Tbm;
	double eval_freq, time;
	double fmax_search = fpeak;
	double fmin_search = fpeak;
	double T=0;
	while(T<Tbm){
		fmin_search = .1*fmin_search;
		//time_phase_corrected_autodiff(&time, 1, &fmin_search, params, 
		//	generation_method, false);
		time_phase_corrected(&time, 1, &fmin_search, params, 
			generation_method, false);
		//T = time_peak- time;
		T = -time_peak+ time;
		//std::cout<<"Time  "<<time<<std::endl;	
	
	}
	while(continue_search)
	{
		//Bisection in log freq space
		eval_freq=sqrt(fmax_search*fmin_search);
		//time_phase_corrected_autodiff(&time, 1, &eval_freq, params, 
		//	generation_method, true);
		time_phase_corrected(&time, 1, &eval_freq, params, 
			generation_method, true);
		//T = time_peak-time;	
		T = -time_peak+time;	
		//The function can be so steep the algorithm cannot determine a valid 
		//frequency with the required tolerance because of floating point error
		//check to see if the difference is at all meaningful
		if(2*(fmax_search - fmin_search)/(fmax_search+fmin_search) <1e-12){
			continue_search =false;
			*freq=eval_freq;

		}
		if(T > (Tbmp )){
			fmin_search = eval_freq;
		}
		else if(T< (Tbmm)){
			fmax_search = eval_freq;
		}
		else{
			continue_search =false;
			*freq=eval_freq;
			//std::cout<<T/T_year<<std::endl;
		}
	}
}
/*! \brief Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now
 *
 * If no time before merger satisfies the requirements, both are set to -1
 */
void threshold_times(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double f_lower,/**<Lower bound of search*/
	double f_upper,/**<upper bound of search*/
	std::string SN,/**< Noise curve name*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output frequencies */
	double tolerance /**< Percent tolerance on SNR search*/
	)
{
	int length = (f_upper-f_lower)/T_obs;
	double deltaf = 1./T_obs;
	double *freqs = new double[length];
	for(int i = 0 ; i<length; i++){
		freqs[i]=f_lower +i*deltaf;
	}
	double *SN_curve = new double[length];
	populate_noise(freqs, SN, SN_curve, length);
	for(int i = 0 ; i<length; i++){
		SN_curve[i] = SN_curve[i]*SN_curve[i];	
	}
	threshold_times(params, generation_method, T_obs, T_wait, freqs, SN_curve, length, SNR_thresh, threshold_times_out,tolerance);
	delete [] freqs;
	delete [] SN_curve;

}

/*! \brief Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now -- No angular dependence used ( only uses plus polarization -- assumes iota = psi = 0 )
 *
 * Assumes this is for multiband -- ie stellar mass BHs -- Only uses pn approximation of time frequency relation
 *
 * If no time before merger satisfies the requirements, both are set to -1
 */
void threshold_times(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double *freqs,/**<Maximum frequency array*/
	double *SN,/**< Noise curve array, should be prepopulated from f_lower to f_upper with spacing 1./T_obs*/
	int length,/**< Length of maximum frequency array*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output frequencies */
	double tolerance /**< Percent tolerance on SNR search*/
	)
{
	if(!params->sky_average){ std::cout<<"NOT sky averaged -- This is not supported by threshold_freqs"<<std::endl;}
	
	params->sky_average = false;
	double bounds[2];
	
	//Max number of iterations -- safety net
	int max_iter = 100;
	int ct = 0;
	double deltaf = freqs[1]-freqs[0];
	double chirpmass = calculate_chirpmass(params->mass1, params->mass2)*MSOL_SEC;
	std::complex<double> *hplus= new std::complex<double>[length];
	std::complex<double> *hcross= new std::complex<double>[length];
	bool not_found = true;
	int bound_id_lower = 0, bound_id_upper = length-1;//Current ids
	int bound_id_lower_prev = 0, bound_id_upper_prev = length-1;//Current ids
	double t_mer=0; //Time before merger
	double snr, snr_prev;
	int precalc_wf_id;
	//Determine if any t_mer between [0,T_wait] allows for an SNR>SNR_thresh 
	
	//Frequency T_obs before it leaves band
	//Using PN f(t) instead of local, because its faster and simpler -- maybe upgrade later
	//Shouldn't matter this far from merger
	t_mer= t_0PN(freqs[bound_id_upper], chirpmass)+T_obs;
	bound_id_lower = (f_0PN(t_mer, chirpmass)- freqs[0])/deltaf;

	fourier_waveform(&freqs[bound_id_lower], bound_id_upper-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
	precalc_wf_id = bound_id_lower;
	snr = std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
	snr_prev=snr;
	bound_id_lower_prev = bound_id_lower;
	bound_id_upper_prev = bound_id_lower;
	double t1 = t_mer, t2=t_mer;
	if(snr>SNR_thresh){not_found= false;}
	else{
		bool bound_search=true, t1_moved=true,t2_moved=true;
		while(bound_search){
			t_mer *=2.;
			bound_id_lower = (f_0PN(t_mer, chirpmass)- freqs[0])/deltaf;
			bound_id_upper = (f_0PN(t_mer-T_obs, chirpmass)- freqs[0])/deltaf;

			fourier_waveform(&freqs[bound_id_lower], bound_id_lower_prev-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
			precalc_wf_id = bound_id_lower;
			snr = std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			if(snr<snr_prev){bound_search=false;}	
		}	
		snr_prev=snr;
		t2 = t_mer;
		while(not_found &&ct < max_iter && t_mer < T_wait && (t2-t1)>T_day){
			t_mer = (t1+t2)/2.;
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(snr>SNR_thresh){not_found=false;}
			else{
				if(t1_moved){
					if(snr>snr_prev){t1=t_mer;t1_moved=true; t2_moved=false;}
					else{t2=t_mer;t2_moved=true; t1_moved=false;;}
				}
				if(t2_moved){
					if(snr>snr_prev){t2=t_mer;t2_moved=true; t1_moved=false;}
					else{t1=t_mer;t1_moved=true; t2_moved=false;;}

				}
			}
			snr_prev = snr;

		}
	}	
	
	//If no SNR is larger than threshold, return 
	if(not_found){
		threshold_times_out[0] = -1;
		threshold_times_out[1] = -1;
	}
	//Find roots
	else{
		double t_save = t_mer;
		ct=0;
		bool found_lower_root=false;	
		bool found_upper_root=false;	
		t1=t_save, t2=t_save;//We know t_mer is over the threshold
		//Find a lower bound for bisection search
		while(snr>SNR_thresh){
			t1/=2.;
			bound_id_lower = ( f_0PN(t1 ,chirpmass) - freqs[0])/deltaf;
			if(t1-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t1 - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
		}
		while(!found_lower_root){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_lower_root=true;threshold_times_out[0]=t_mer;}
			else{
				if(snr>SNR_thresh){ t2 = t_mer;	}
				else{ t1=t_mer;}
			}
			snr_prev=snr;
			
		}
		ct=0;
		t1=t_save; t2=t_save;
		do{
			t2*=2.;
			if(t2>T_wait){t2=T_wait;}	
			bound_id_lower = ( f_0PN(t2 ,chirpmass) - freqs[0])/deltaf;
			if(t2-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t2 - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			if(bound_id_lower < precalc_wf_id){
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			if(t2==T_wait && snr>SNR_thresh){ found_upper_root=true; threshold_times_out[1]=T_wait;break;}
		}while(snr>SNR_thresh  );
		while(!found_upper_root){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			bound_id_lower = ( f_0PN(t_mer ,chirpmass) - freqs[0])/deltaf;
			if(t_mer-T_obs > 0){
				bound_id_upper = ( std::min( f_0PN(  t_mer - T_obs,chirpmass) , freqs[length-1] ) - freqs[0])/deltaf;
			}
			else{
				bound_id_upper = (  freqs[length-1]  - freqs[0])/deltaf;

			}
			//Update waveform if using new frequencies
			if(bound_id_lower < precalc_wf_id){
				fourier_waveform(&freqs[bound_id_lower], precalc_wf_id-bound_id_lower, &hplus[bound_id_lower], &hcross[bound_id_lower], generation_method,params);
				precalc_wf_id = bound_id_lower;
			}
			bound_id_lower_prev= bound_id_lower;
			bound_id_upper_prev= bound_id_upper;
			snr =std::sqrt(1.)*calculate_snr_internal(&SN[bound_id_lower], &hplus[bound_id_lower],&freqs[bound_id_lower],bound_id_upper-bound_id_lower);
			ct++;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_upper_root=true;threshold_times_out[1]=t_mer;}
			else{
				if(snr>SNR_thresh){ t1 = t_mer;	}
				else{ t2=t_mer;}
			}
			snr_prev=snr;
			
		}
	}
	//Cleanup
	delete [] hplus;
	delete [] hcross;
}


/*! \brief Utility for calculating the threshold times before merger that result in an SNR>SNR_thresh --GSL quad integration implementation
 *
 * See arXiv 1902.00021
 *
 * Binary must merge within time T_wait
 *
 * SNR is calculated with frequencies [f(t_mer),f(t_mer-T_obs)] or [f(t_mer),0] depending on whether the binary has merged or not
 *
 * Assumes sky average -- Only supports PhenomD for now -- No angular dependence used ( only uses plus polarization -- assumes iota = psi = 0 )
 *
 * Assumes this is for multiband -- ie stellar mass BHs -- Only uses pn approximation of time frequency relation
 *
 * If no time before merger satisfies the requirements, both are set to -1
 *
 * ALL temporal quantities in seconds or Hz
 *
 * Return values: 
 * 	
 * 	0 -- success
 * 	
 * 	11-- Failure: SNR was 0 in lower bound
 *
 * 	12-- Failure: SNR was 0 in upper bound
 *
 * 	13 -- partial success: Closest values output, but roundoff error prevented the routine from reaching the desired accuracy
 */
int threshold_times_gsl(gen_params_base<double> *params,
	std::string generation_method, /**<Generation method to use for the waveform*/
	double T_obs, /**<Observation time -- also specifies the frequency spacing (\delta f = 1./T_obs)*/
	double T_wait, /**<Wait time -- Maximum time for binaries to coalesce */
	double fmin,/**<Maximum frequency array*/
	double fmax,/**<Maximum frequency array*/
	std::string SN,/**< Noise curve array, should be prepopulated from f_lower to f_upper with spacing 1./T_obs*/
	double SNR_thresh, /**< Threshold SNR */
	double *threshold_times_out,/**<[out] Output times */
	double tolerance, /**< Percent tolerance on SNR search*/
	gsl_integration_workspace *w,
	int np
	)
{
	int fail_ct=0;
	int fail_THRESH = 10; //Attempts to fix an error due to roundoff
	double SHIFT_UP=1.0 + 1e-2;
	double SHIFT_DOWN=1.0 - 1e-2;

	bool round_off_error=false;	
	
	bool save_SA = params->sky_average;
	if(!params->sky_average){ std::cout<<"NOT sky averaged -- This is not supported by threshold_freqs"<<std::endl;}
	params->sky_average = false;

	bool stellar_mass = true;	
	double fpeak, fdamp,fRD;
	postmerger_params(params,generation_method,&fpeak,&fdamp, &fRD);
	
	if(fpeak < fmax){
		stellar_mass=false;
	}
	stellar_mass=true;

	double bounds[2];
	
	//Max number of iterations -- safety net
	double chirpmass = calculate_chirpmass(params->mass1, params->mass2)*MSOL_SEC;
	bool not_found = true;
	double f_lower = fmin, f_upper = fmax;//Current ids
	double f_lower_prev = fmin, f_upper_prev = fmax;//Current ids
	double t_mer=0; //Time before merger
	double snr, snr_prev;
	double rel_err = tolerance;
	//Determine if any t_mer between [0,T_wait] allows for an SNR>SNR_thresh 
	
	//Frequency T_obs before it leaves band
	//Using PN f(t) instead of local, because its faster and simpler -- maybe upgrade later
	//Shouldn't matter this far from merger
	if(stellar_mass){
		t_mer= t_0PN(f_upper, chirpmass)+T_obs;
		f_lower = f_0PN(t_mer,chirpmass);	
	}
	else{
		f_upper = fpeak*1.1;	
		f_lower = f_0PN(T_obs,chirpmass);	
	}

	snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
	snr_prev=snr;
	double t1 = t_mer, t2=t_mer;
	if(snr>SNR_thresh){not_found= false;}
	else{
		bool bound_search=true, t1_moved=true,t2_moved=true;
		while(bound_search){
			t_mer *=2.;
			if(stellar_mass){
				f_lower = (f_0PN(t_mer, chirpmass));
				f_upper = (f_0PN(t_mer-T_obs, chirpmass));
			}
			else{
				//Much much slower
				Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance);
				Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance);
			}
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			if(snr<snr_prev){bound_search=false;}	
			else if(f_lower < fmin){ 
				t_mer = t_0PN(fmin, chirpmass); 
				bound_search=false;
			}
		}	
		snr_prev=snr;
		t2 = t_mer;
		while(not_found && t_mer < T_wait && (t2-t1)>T_day){
			t_mer = (t1+t2)/2.;
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				else{
					f_upper =   fmax;

				}
			}
			else{
				Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance);
				if(t_mer-T_obs > 0){
					Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance);
				}
				else{
					f_upper =   fpeak*1.1;

				}

			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			if(snr>SNR_thresh){not_found=false;}
			else{
				if(t1_moved){
					if(snr>snr_prev){t1=t_mer;t1_moved=true; t2_moved=false;}
					else{t2=t_mer;t2_moved=true; t1_moved=false;}
				}
				if(t2_moved){
					if(snr>snr_prev){t2=t_mer;t2_moved=true; t1_moved=false;}
					else{t1=t_mer;t1_moved=true; t2_moved=false;}

				}
			}
			snr_prev = snr;

		}
	}	
	//If no SNR is larger than threshold, return 
	if(not_found){
		threshold_times_out[0] = -1;
		threshold_times_out[1] = -1;
	}
	//Find roots
	else{
		double t_save = t_mer;
		bool found_lower_root=false;	
		bool found_upper_root=false;	
		double t1_prev=t_save, t2_prev=t_save;
		t1=t_save, t2=t_save;//We know t_mer is over the threshold
		//Find a lower bound for bisection search
		while(snr>SNR_thresh){
			t1/=2.;
			if(stellar_mass){
				f_lower =  f_0PN(t1 ,chirpmass) ;
				if(t1-T_obs > 0){
					f_upper =  std::min( f_0PN(  t1 - T_obs,chirpmass) , fmax ) ;
				}
				else{
					f_upper =   fmax;

				}
			}
			else{
				Tbm_to_freq(params, generation_method, t1, &f_lower,tolerance);
				if(t1-T_obs > 0){
					Tbm_to_freq(params, generation_method, t1-T_obs, &f_upper,tolerance);
				}
				else{
					f_upper =   1.1*fpeak;

				}
			}
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			//std::cout<<f_lower<<" "<<f_upper<<std::endl;
		}
		while(!found_lower_root){
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax);
				}
				else{
					f_upper = fmax;

				}
			}
			else{
				Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance);
				if(t_mer-T_obs > 0){
					Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance);
				}
				else{
					f_upper =   1.1*fpeak;

				}

			}
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			//std::cout<<snr<<std::endl;
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_lower_root=true;threshold_times_out[0]=t_mer;}
			else if( (fabs(t1-t2)*2/fabs(t1+t2)<1e-14) ){
				found_lower_root=true;
				threshold_times_out[0]=t_mer;	
				round_off_error = true;
			}
			else{
				//Sometimes, integration messes up
				if(snr==0){
					//If the SNR hits 0 for the lower root, its typically because the 
					//Waveform is sharply dropping off post merger, and this algorithm can't work
					//with such a steep function. But if this is the case, the error is on the order of 
					//seconds, which is completely negligible for most things this routine is used for
					found_lower_root=true; threshold_times_out[0]=(t1+t2 )/2.;
					//if(fail_ct>fail_THRESH){
					//	threshold_times_out[0] = -1;
					//	threshold_times_out[1] = -1;
					//	return 11;
					//}
					//else{
					//	if(f_lower_prev == f_lower){
					//		t2 = t2_prev*SHIFT_UP;
					//	}
					//	else{
					//		t1 = t1_prev*SHIFT_DOWN;
					//	}
					//	fail_ct++;
					//}
					
				}
				else{
					if(snr>SNR_thresh){ t2 = t_mer;	}
					else{ t1=t_mer;}
				}
			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			t1_prev= t1;
			t2_prev= t2;
			snr_prev=snr;
			//std::cout<<"LOWER: "<<f_lower<<" "<<f_upper<<" "<<t1<<" "<<t2<<" "<<snr<<std::endl;
		}
		fail_ct =0;
		t1=t_save; t2=t_save;
		do{
			t2*=2.;
			if(stellar_mass){
				if(t2>T_wait){t2=T_wait;}	
				f_lower =  f_0PN(t2 ,chirpmass) ;
				if(t2-T_obs > 0){
					f_upper =  std::min( f_0PN(  t2 - T_obs,chirpmass) , fmax ) ;
				}
				else{
					f_upper =fmax;

				}
			}
			else{
				Tbm_to_freq(params, generation_method, t2, &f_lower,tolerance);
				if(t2-T_obs > 0){
					Tbm_to_freq(params, generation_method, t2-T_obs, &f_upper,tolerance);
				}
				else{
					f_upper =   1.1*fpeak;

				}
			}
			//std::cout<<"UPPER SEARCH: "<<f_upper<<" "<<f_lower<<" "<<t1<<" "<<t2<<T_wait<<std::endl;
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			if(t2==T_wait && snr>SNR_thresh){ found_upper_root=true; threshold_times_out[1]=T_wait;break;}
			//std::cout<<"UPPER SEARCH: "<<snr<<std::endl;
		}while(snr>SNR_thresh  );
		while(!found_upper_root){
			
			t_mer = (t1+t2)/2.;
			//Find new frequency bound ids
			if(stellar_mass){
				f_lower =  f_0PN(t_mer ,chirpmass) ;
				if(t_mer-T_obs > 0){
					f_upper =  std::min( f_0PN(  t_mer - T_obs,chirpmass) , fmax );
				}
				else{
					f_upper = fmax;

				}
			}
			else{
				Tbm_to_freq(params, generation_method, t_mer, &f_lower,tolerance);
				if(t_mer-T_obs > 0){
					Tbm_to_freq(params, generation_method, t_mer-T_obs, &f_upper,tolerance);
				}
				else{
					f_upper =   1.1*fpeak;

				}
			}
			f_lower_prev= f_lower;
			f_upper_prev= f_upper;
			snr = snr_threshold_subroutine(	f_lower, f_upper, rel_err,params, generation_method,SN, w,np);
			if(std::abs(snr-SNR_thresh)/SNR_thresh<tolerance ){found_upper_root=true;threshold_times_out[1]=t_mer;}
			else if( (fabs(t1-t2)*2/fabs(t1+t2)<1e-14) ){
				found_upper_root=true;
				threshold_times_out[1]=t_mer;	
				round_off_error = true;
			}
			else{
				if(snr==0){
					threshold_times_out[0] = -1;
					threshold_times_out[1] = -1;
					params->sky_average = save_SA;
					return 12;
					//if(fail_ct>fail_THRESH){
					//	threshold_times_out[0] = -1;
					//	threshold_times_out[1] = -1;
					//	return 12;
					//}
					//else{
					//	if(f_lower_prev == f_lower){
					//		t2 = t2_prev*SHIFT_UP;
					//	}
					//	else{
					//		t1 = t1_prev*SHIFT_DOWN;
					//	}
					//	fail_ct++;
					//}
					
				}
				else{
					if(snr>SNR_thresh){ t1 = t_mer;	}
					else{ t2=t_mer;}
				}
			}
			snr_prev=snr;
			//std::cout<<"Upper: "<<f_lower<<" "<<f_upper<<" "<<snr<<std::endl;
			
		}
	}
	if(round_off_error){
		params->sky_average = save_SA;
		return 13;
	}
	params->sky_average = save_SA;
	//std::cout<<std::endl;
	return 0;
}
double snr_threshold_subroutine(double fmin, double fmax, double rel_err, gen_params_base<double> *params, std::string generation_method,std::string SN, gsl_integration_workspace *w, int np)
{
	gsl_snr_struct helper_params;
	helper_params.params = params;
	helper_params.generation_method = generation_method;
	helper_params.SN = SN;
	gsl_function F;
	F.function = [](double f, void *param){return integrand_snr_SA_subroutine(f,param);};
	F.params = (void *)&helper_params;
	double result, err;
	int errcode = gsl_integration_qag(&F,fmin, fmax, 0,rel_err, np, GSL_INTEG_GAUSS15,w, &result, &err);
	return sqrt(result);
}

//###########################################################################
//template void map_extrinsic_angles<double>(gen_params_base<double> *);
//template void map_extrinsic_angles<adouble>(gen_params_base<adouble> *);

template void  time_phase_corrected<double>(double *, int, double *, gen_params_base<double> *, std::string, bool );
template void  time_phase_corrected<adouble>(adouble *, int, adouble *, gen_params_base<adouble> *, std::string, bool);

template int fourier_detector_response_horizon<double>(double *, int, std::complex<double> *, std::complex<double> *,std::complex<double> *, double, double, std::string);
template int fourier_detector_response_horizon<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *,std::complex<adouble> *, adouble, adouble, std::string);
//
template int fourier_detector_response_horizon<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double, double, std::string);
template int fourier_detector_response_horizon<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble, adouble, std::string);
//
//
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *);
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *, double *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *, adouble *);
//
template int fourier_detector_response_horizon<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double>*);
template int fourier_detector_response_horizon<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble>*);
//
template int fourier_detector_response_equatorial<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double , double, double,double *, double, double, double, double, std::string);
template int fourier_detector_response_equatorial<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble , adouble, double,adouble*, adouble, adouble, adouble, adouble,  std::string);
