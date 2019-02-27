#include "waveform_util.h"
#include "util.h"
#include "waveform_generator.h"
#include "noise_util.h"
#include <fftw3.h>
#include <algorithm>
#include <complex>
#include <vector>



std::complex<double> Q(double theta, double phi, double iota)
{
	double ct = cos(theta);
	double cp2 = cos(2.*phi);
	double sp2 = sin(2.*phi);
	double ci = cos(iota);

	double Fplus = (1./2)*(1+ ct*ct)*cp2;
	double Fcross = ct * sp2;
	std::complex<double> Q = (1+ci*ci)/2. *Fplus + std::complex<double>(0,Fcross*ci);
	return Q;
}

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
