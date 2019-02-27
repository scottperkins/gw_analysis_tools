#include "mcmc_routines.h"
#include "waveform_generator.h"
#include "util.h"
#include "waveform_util.h"
#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <iostream>



//#####################################################################
//#####################################################################
//
//MEMORY LEAK IN MCMC ROUTINES -- MUST FIX
//
//#####################################################################
//#####################################################################

/*! \file 
 * Routines for implementation in MCMC algorithms
 *
 * */

/*!\brief Function to calculate the log Likelihood as defined by -1/2 (d-h|d-h) maximized over the extrinsic parameters phic and tc
 *
 * frequency array must be uniform spacing - this shouldn't be a problem when working with real data as DFT return uniform spacing
 */
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan
				)
{
	//Make template waveform
	gen_params params;
	params.mass1 = calculate_mass1(chirpmass, symmetric_mass_ratio);	
	params.mass2 = calculate_mass2(chirpmass, symmetric_mass_ratio);	
	params.Luminosity_Distance = 500;	
	std::complex<double> *template_strain = (std::complex<double>*)malloc(sizeof(std::complex<double>)*length);	
	params.spin1[0]= 0;
	params.spin1[1]= 0;
	params.spin1[2]= spin1;
	params.spin2[0]= 0;
	params.spin2[1]= 0;
	params.spin2[2]= spin2;
	params.phic=0;
	params.tc=0;
	params.NSflag=NSflag;
	//fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", m1,m2,Dl,spin1vec,
	//				spin2vec);
	fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", &params);
	
	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[1]-frequencies[0];
	double sum = 0;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(template_strain[i]*std::conj(template_strain[i]))/noise[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double normalizing_factor = SNR/sqrt(integral);
	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	//std::complex<double> g_tilde[length];
	std::complex<double> g_tilde;
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*normalizing_factor*conj(data[i]) * template_strain[i] / noise[i]; 
		plan->in[i][0] = real(g_tilde);
		plan->in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double) *length);
	
	fftw_execute(plan->p);
	
	for (int i=0;i<length; i++)
	{
		g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
	}

	double max = *std::max_element(g, g+length)*delta_f; 
	
	free(template_strain);
	free(g);
	free(integrand);
	return -(SNR*SNR - max);
}
				
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag)
{
	fftw_outline plan;
	initiate_likelihood_function(&plan, length);
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD(frequencies,
				length,
				data,
				noise,
				SNR,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				NSflag,
				&plan);
	deactivate_likelihood_function(&plan);
	free(data);
	return out;
}
double maximized_coal_log_likelihood_IMRPhenomD(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double SNR,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				bool NSflag,
				fftw_outline *plan)
{
	std::complex<double> *data = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD(frequencies,
				length,
				data,
				noise,
				SNR,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				NSflag,
				plan);
	free(data);
	return out;
}
				
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				int length,
				std::complex<double> *data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan)
{

	//Make template waveform
	gen_params params;
	params.mass1 = calculate_mass1(chirpmass, symmetric_mass_ratio);	
	params.mass2 = calculate_mass2(chirpmass, symmetric_mass_ratio);	
	params.Luminosity_Distance = Luminosity_Distance;	
	std::complex<double> *template_strain = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);	
	params.spin1[0]= 0;
	params.spin1[1]= 0;
	params.spin1[2]= spin1;
	params.spin2[0]= 0;
	params.spin2[1]= 0;
	params.spin2[2]= spin2;
	params.phic=0;
	params.tc=0;
	params.NSflag=NSflag;

	fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", &params);
	std::complex<double> q = Q(theta,phi,iota);
	for (int i = 0;i<length;i++)
		template_strain[i] = template_strain[i]*q;
	

	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[1]-frequencies[0];
	double sum = 0;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(template_strain[i]*std::conj(template_strain[i]))/noise[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	sum = 0;
	double integral2;
	for (int i =0;i< length;i++)
		integrand[i] = real(data[i]*std::conj(data[i]))/noise[i];
	//integral2 = 4*trapezoidal_sum_uniform(delta_f, length, integrand2);
	integral2 = 4*simpsons_sum(delta_f, length, integrand);
	double DD = integral2;

	//double normalizing_factor = SNR/sqrt(integral);
	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	std::complex<double> g_tilde;
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * template_strain[i] / noise[i]; 
		plan->in[i][0] = real(g_tilde);
		plan->in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double)*length);
	
	fftw_execute(plan->p);
	
	for (int i=0;i<length; i++)
	{
		g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
	}

	double max = *std::max_element(g, g+length)*delta_f; 

	free(template_strain);
	free(integrand);
	free(g);

	return -0.5*(DD+HH- 2*max);
}
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag)
{
	fftw_outline plan;
	initiate_likelihood_function(&plan,length);
	std::complex<double> *data = (std::complex<double> *) malloc(sizeof(std::complex<double> ) *length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD_Full_Param(frequencies,
				length,
				data,
				noise,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				Luminosity_Distance,
				theta,
				phi,
				iota,
				NSflag,
				&plan);
	deactivate_likelihood_function(&plan);

	free(data);

	return out;
}
double maximized_coal_log_likelihood_IMRPhenomD_Full_Param(double *frequencies,
				size_t length,
				double *real_data,
				double *imag_data,
				double *noise,
				double chirpmass,	/**< in solar masses*/
				double symmetric_mass_ratio, 
				double spin1,
				double spin2,
				double Luminosity_Distance,
				double theta,
				double phi,
				double iota,
				bool NSflag,
				fftw_outline *plan)
{
	std::complex<double> *data = (std::complex<double> *) malloc(sizeof(std::complex<double> ) *length);
	for (int i =0; i<length; i++)
		data[i] = std::complex<double>(real_data[i],imag_data[i]);
	double out =  maximized_coal_log_likelihood_IMRPhenomD_Full_Param(frequencies,
				length,
				data,
				noise,
				chirpmass,	
				symmetric_mass_ratio, 
				spin1,
				spin2,
				Luminosity_Distance,
				theta,
				phi,
				iota,
				NSflag,
				plan);

	free(data);

	return out;
}

void initiate_likelihood_function(fftw_outline *plan, int length)
{
	plan->in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	plan->p = fftw_plan_dft_1d(length, plan->in, plan->out,FFTW_FORWARD, FFTW_MEASURE);
}
void deactivate_likelihood_function(fftw_outline *plan)
{
	fftw_destroy_plan(plan->p);
	fftw_free(plan->in);
	fftw_free(plan->out);
	fftw_cleanup();
}
