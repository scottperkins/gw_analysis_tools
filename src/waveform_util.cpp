#include "waveform_util.h"
#include "util.h"
#include "waveform_generator.h"
#include "noise_util.h"
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



double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				std::complex<double> *data,
				std::string detector,
				std::string generation_method,
				gen_params param
				)
{
	
	/*produce noise curve*/
	double *noise = (double *)malloc(sizeof(double)*length);	
	populate_noise(frequencies,detector, noise, length);
	for(int i = 0; i<length; i++)
		noise[i] = noise[i]*noise[i];

	std::complex<double> q = Q(param.theta,param.phi,param.incl_angle);
	std::complex<double> *detector_response
			 = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	double *integrand
			 = (double *)malloc(sizeof(double)*length);
	/*produce the waveform*/
	fourier_waveform(frequencies, length, detector_response, generation_method, &param);

	/*Calculate the template snr integrand 4*Re(h* h /S(f)) - factor of q for the plus, cross modes 
 * 	effect on the detector*/
	for (int i = 0; i<length;i++)
		integrand[i] = 4.*real(conj(q*detector_response[i])*q*detector_response[i])/noise[i]; 
	double delta_f = frequencies[1]-frequencies[0];
	double snr_template;
	//snr_template = sqrt(trapezoidal_sum_uniform(delta_f,length, integrand));
	snr_template = sqrt(simpsons_sum(delta_f,length, integrand));

	fftw_complex *in, *out; 
	fftw_plan p;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
        p = fftw_plan_dft_1d(length, in, out,FFTW_FORWARD, FFTW_MEASURE);
	std::complex<double> g_tilde;
        for (int i=0;i<length; i++)
        {
                g_tilde = 4.*q*conj(data[i]) * detector_response[i] / noise[i];
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
	return max/snr_template;
}
double data_snr_maximized_extrinsic(double *frequencies,
				int length,
				double *data_real,
				double *data_imag,
				std::string detector,
				std::string generation_method,
				gen_params param
				)
{
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>) * length);
        for (int i =0; i<length; i++)
                data[i] = std::complex<double>(data_real[i],data_imag[i]);	
	double snr;
	snr = data_snr_maximized_extrinsic(frequencies,
				length,
				data,
				detector,
				generation_method,
				param);
	free(data);
	return snr;
}
/*! \brief Caclulates the snr given a detector and waveform (complex) and frequencies
 *      
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

/* \brief calculates the detector response for a given waveform and detector
 */
int fourier_detector_response(double *frequencies, /**<array of frequencies corresponding to waveform*/
			int length,/**< length of frequency/waveform arrays*/
			std::complex<double> *hplus, /*<precomputed plus polarization of the waveform*/ 
			std::complex<double> *hcross, /**<precomputed cross polarization of the waveform*/ 
			std::complex<double> *detector_response, /**< [out] detector response*/
			double theta, /**< polar angle theta in detector frame*/
			double phi, /**< azimuthal angle phi in detector frame*/ 
			double iota,/**< inclination angle of the binary*/ 
			std::string detector/**< detector - list of supported detectors in noise_util*/
			)
{
	int status=1;
	double ci = cos(iota);
	double fplus;
	double fcross;
	if(	detector == "LIGO" || 
		detector == "Livingston" || 
		detector == "Hanford" || 
		detector == "VIRGO")
	{
		fplus = right_interferometer_plus(theta,phi);
		fcross = right_interferometer_cross(theta,phi);	
	}
	for (int i =0; i <length; i++)
	{
		detector_response[i] = fplus * (1. + ci*ci)/2. * hplus[i] 
					+ (fcross * ci);
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
 * This is a wrapper that combines generation with response functions: if producing mulitple responses for one waveform (ie stacking Hanford, Livingston, and VIRGO), it will be considerably more efficient to calculate the waveform once, then combine each response manually - detector responses will be in the noise_util files
 */
int fourier_detector_response(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/ int length,/**<integer length of all the arrays*/
			std::complex<double> *response, /**< [out] complex array for the output plus polarization waveform*/
			std::string detector,
			std::string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	int status = 1;
	//generate waveform
	std::complex<double> *waveform_plus =
		(std::complex<double> *)malloc(sizeof(std::complex<double>) * length);
	std::complex<double> *waveform_cross=
		(std::complex<double> *)malloc(sizeof(std::complex<double>) * length);
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
			parameters->incl_angle, 
			detector 
			) ;
		
	//Deallocate memory
	free(waveform_plus);
	free(waveform_cross);

	return status;
}
