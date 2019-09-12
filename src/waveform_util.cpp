#include "waveform_util.h"
#include "util.h"
#include "waveform_generator.h"
#include "detector_util.h"
#include <fftw3.h>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
/*!\file 
 * Utilities for waveforms - SNR calculation and detector response
 * 	
 * includes snr and detector response
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
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	T Fplus, Fcross;
	
	detector_response_functions_equatorial(detector, ra, dec, psi, gmst, &Fplus, &Fcross);
	
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
int fourier_detector_response_equatorial(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ int length,/**<integer length of all the arrays*/
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
	status = fourier_detector_response_equatorial(frequencies, 
			length, 
			waveform_plus, 
			waveform_cross,
			response, 
			parameters->RA, 
			parameters->DEC, 
			parameters->psi, 
			parameters->gmst, 
			detector 
			) ;
	
		
	//Deallocate memory
	free(waveform_plus);
	free(waveform_cross);

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
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::complex<double> *,std::complex<double> *, double, double, std::string);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *,std::complex<adouble> *, adouble, adouble, std::string);
//
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double, double, std::string);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble, adouble, std::string);
//
//
template int fourier_detector_response_equatorial<double>(double *, int , std::complex<double> *,std::string, std::string, gen_params_base<double> *);
template int fourier_detector_response_equatorial<adouble>(adouble *, int , std::complex<adouble> *,std::string, std::string, gen_params_base<adouble> *);
//
template int fourier_detector_response<double>(double *, int, std::complex<double> *, std::string, std::string, gen_params_base<double>*);
template int fourier_detector_response<adouble>(adouble *, int, std::complex<adouble> *, std::string, std::string, gen_params_base<adouble>*);
//
template int fourier_detector_response_equatorial<double>(double *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *, double, double , double, double, std::string);
template int fourier_detector_response_equatorial<adouble>(adouble *, int, std::complex<adouble> *, std::complex<adouble> *, std::complex<adouble> *, adouble, adouble , adouble, double, std::string);
