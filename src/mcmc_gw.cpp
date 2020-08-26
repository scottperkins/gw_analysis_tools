#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "util.h"
#include "io_util.h"
#include "detector_util.h"
#include "ppE_utilities.h"
#include "waveform_util.h"
#include "ortho_basis.h"
#include "fisher.h"
#include "mcmc_sampler.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <mutex>
#include <thread>
#include <condition_variable>

/*! \file 
 * Routines for implementation in MCMC algorithms specific to GW CBC analysis
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
	params.phiRef=0;
	params.tc=0;
	params.NSflag1=NSflag;
	params.NSflag2=NSflag;
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
	allocate_FFTW_mem_forward(&plan, length);
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
	deallocate_FFTW_mem(&plan);
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
	params.phiRef=0.;
	params.tc=0.0;
	params.NSflag1=NSflag;
	params.NSflag2=NSflag;

	fourier_waveform(frequencies, length, template_strain, "IMRPhenomD", &params);
	std::complex<double> q = Q(theta,phi,iota);
	for (int i = 0;i<length;i++)
		template_strain[i] = template_strain[i]*q;
	

	////Calculate template snr and scale it to match the data snr
	////later, upgrade to non uniform spacing, cause why not
	//double delta_f = frequencies[1]-frequencies[0];
	//double sum = 0;
	//double *integrand = (double *)malloc(sizeof(double)*length);
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(template_strain[i]*std::conj(template_strain[i]))/noise[i];
	////double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	//double integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double HH = integral;

	//sum = 0;
	//double integral2;
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(data[i]*std::conj(data[i]))/noise[i];
	////integral2 = 4*trapezoidal_sum_uniform(delta_f, length, integrand2);
	//integral2 = 4*simpsons_sum(delta_f, length, integrand);
	//double DD = integral2;

	////###################################################################
	////testing
	////###################################################################
	//sum = 0;
	//double integral3;
	//for (int i =0;i< length;i++)
	//	integrand[i] = real(template_strain[i]*std::conj(data[i]))/noise[i];
	////integral2 = 4*trapezoidal_sum_uniform(delta_f, length, integrand2);
	//integral3 = 4*simpsons_sum(delta_f, length, integrand);
	//double HD = integral3;
	////###################################################################



	////double normalizing_factor = SNR/sqrt(integral);
	//
	////calculate the fourier transform that corresponds to maximizing over phic and tc
	////Use malloc at some point, not sure how long these arrays will be 
	//std::complex<double> g_tilde;
	//for (int i=0;i<length; i++)
	//{
	//	g_tilde = 4.*conj(data[i]) * template_strain[i] / noise[i]; 
	//	plan->in[i][0] = real(g_tilde);
	//	plan->in[i][1] = imag(g_tilde);
	//}

	//double *g = (double *)malloc(sizeof(double)*length);
	//
	//fftw_execute(plan->p);
	//
	//for (int i=0;i<length; i++)
	//{
	//	g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
	//}

	//double max = *std::max_element(g, g+length)*delta_f; 

	//free(integrand);
	//free(g);

	double out =  maximized_Log_Likelihood_aligned_spin_internal(data,
				noise,
				frequencies,
				template_strain,
				length,
				plan
				);
	free(template_strain);
	return out;
	//return -0.5*(HH- 2*HD);//TESTING
	//return -0.5*(HH- 2*max);
	//return -0.5*(DD+HH- 2*max);
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
	allocate_FFTW_mem_forward(&plan,length);
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
	deallocate_FFTW_mem(&plan);

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



/*! \brief routine to maximize over all extrinsic quantities and return the log likelihood
 *
 * IMRPhenomD -- maximizes over DL, phic, tc, \iota, \phi, \theta
 * IMRPhenomP -- maximizes over DL, phic,tc, \psi, \phi , \theta
 */
double maximized_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				//gen_params *params,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	double ll = 0;
	//if(	generation_method  == "IMRPhenomD" ||
	//	generation_method  == "ppE_IMRPhenomD_Inspiral" ||
	//	generation_method  == "ppE_IMRPhenomD_IMR" ){
	if(	generation_method.find("IMRPhenomD")!=std::string::npos){
		std::complex<double> *response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		fourier_detector_response_horizon(frequencies, length, response, detector, generation_method, params);
		ll = maximized_Log_Likelihood_aligned_spin_internal(data, psd, frequencies,response, length, plan);
		
		free(response);
		}
	//else if ( 	generation_method == "IMRPhenomPv2" ||
	//		generation_method == "ppE_IMRPhenomPv2_Inspiral" ||
	//		generation_method == "ppE_IMRPhenomPv2_IMR" ){
	if(	generation_method.find("IMRPhenomPv2")!=std::string::npos){
				
		std::complex<double> *hp =
				(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		std::complex<double> *hc =
				(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
		fourier_waveform(frequencies,length,hp,hc,generation_method,params);
		ll = maximized_Log_Likelihood_unaligned_spin_internal(data,psd,frequencies,hp,hc, length, plan);

		free(hp);
		free(hc);
	}
	return ll;
}

					
double maximized_coal_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				double *template_real,
				double *template_imag,
				fftw_outline *plan
				)
{
	
	std::complex<double> *data = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *h = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	for(int i =0; i<length; i ++){
		data[i] = std::complex<double>(data_real[i],data_imag[i]);
		h[i] = std::complex<double>(template_real[i],template_imag[i]);
	}
	double tc,phic,ll=0;
	ll = maximized_coal_Log_Likelihood_internal(data,
				psd,
				frequencies,
				h,
				length,
				plan,
				&tc,
				&phic);
	
	free(data);
	free(h);
	return ll;
}
double maximized_Log_Likelihood(double *data_real, 
				double *data_imag,
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	
	std::complex<double> *data = 
			(std::complex<double> *) malloc(sizeof(std::complex<double>)*length);
	for(int i =0; i<length; i ++){
		data[i] = std::complex<double>(data_real[i],data_imag[i]);
	}
	
	double ll = maximized_Log_Likelihood(data, 
				psd,
				frequencies,
				length,
				params,
				detector,
				generation_method,
				plan);
	free(data);
	return ll;
}

/*! \brief Function to maximize only over coalescence variables tc and phic, returns the maximum values used
 *
 */
double maximized_coal_Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan,
				double *tc,
				double *phiRef
				)
{
	std::complex<double> *response =
		(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
	fourier_detector_response(frequencies, length, response, detector, generation_method, params);
	double ll = maximized_coal_Log_Likelihood_internal(data, psd, frequencies,
				response, length, plan, tc, phiRef);
	free(response);
	return ll;
}

double maximized_coal_Log_Likelihood_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan,
				double *tc,
				double *phiRef
				)
{

	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[length/2]-frequencies[length/2-1];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	std::complex<double> g_tilde;

	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * detector_response[i] / psd[i]; 
		//g_tilde = conj(data[i]) * detector_response[i] / psd[i]; 
		//plan->in[i][0] = real(g_tilde);
		//plan->in[i][1] = imag(g_tilde);
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double)*length);
	std::complex<double> *gc = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	
	for (int i=0;i<length; i++)
	{
		g[i] = std::abs(std::complex<double>(out[i][0],out[i][1])) ;
		gc[i] = std::complex<double>(out[i][0],out[i][1]) ;
		//g[i] = plan->out[i][0]*plan->out[i][0]+plan->out[i][1]*plan->out[i][1] ;
		//g[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1] ;
	}
	//std::vector<int>::iterator max;
	double *max;
	max = std::max_element(g, g+length); 
	double max_val = *max*delta_f;
	//double max_val = *max;
	int max_index = std::distance(g,max);
	double Tau = 1./delta_f;
	*tc = (double)(max_index)/length * ( Tau );
	*phiRef = std::arg(gc[max_index]);
	//double max = *std::max_element(g, g+length)*delta_f; 

	free(integrand);
	free(g);
	free(gc);
	fftw_free(in);
	fftw_free(out);

	return -0.5*(HH- 2*max_val);
	//return .5*(max*max)/HH;
	//return .5*(max)/HH;
}

/*! \brief Unmarginalized log of the likelihood
 *
 */
double Log_Likelihood(std::complex<double> *data, 
				double *psd,
				double *frequencies,
				size_t length,
				gen_params_base<double> *params,
				std::string detector,
				std::string generation_method,
				fftw_outline *plan
				)
{
	double ll = 0;

	std::complex<double> *detect_response =
			(std::complex<double> *) malloc(sizeof(std::complex<double>) * length);
	fourier_detector_response(frequencies,length,detect_response,detector, generation_method,params);
	ll = Log_Likelihood_internal(data,psd,frequencies,detect_response, length, plan);

	//if(ll>0){
	//}
	free(detect_response);
	return ll;
}
/*! \brief Maximized match over coalescence variables - returns log likelihood NOT NORMALIZED for aligned spins
 *
 * Note: this function is not properly normalized for an absolute comparison. This is made for MCMC sampling, so to minimize time, constant terms like (Data|Data), which would cancel in the Metropolis-Hasting ratio, are left out for efficiency
 */
double maximized_Log_Likelihood_aligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *detector_response,
				size_t length,
				fftw_outline *plan
				)
{
	//Calculate template snr and scale it to match the data snr
	//later, upgrade to non uniform spacing, cause why not
	double delta_f = frequencies[1]-frequencies[0];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//double integral = 4.*trapezoidal_sum_uniform(delta_f, length, integrand);
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	
	//calculate the fourier transform that corresponds to maximizing over phic and tc
	//Use malloc at some point, not sure how long these arrays will be 
	std::complex<double> g_tilde;

	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * detector_response[i] / psd[i]; 
		//g_tilde = conj(data[i]) * detector_response[i] / psd[i]; 
		//plan->in[i][0] = real(g_tilde);
		//plan->in[i][1] = imag(g_tilde);
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}

	double *g = (double *)malloc(sizeof(double)*length);
	
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	
	for (int i=0;i<length; i++)
	{
		//g[i] = std::abs(std::complex<double>(plan->out[i][0],plan->out[i][1])) ;
		//g[i] = plan->out[i][0]*plan->out[i][0]+plan->out[i][1]*plan->out[i][1] ;
		g[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1] ;
	}

	double max = *std::max_element(g, g+length)*delta_f*delta_f; 
	//double max = *std::max_element(g, g+length)*length; 

	free(integrand);
	free(g);
	fftw_free(in);
	fftw_free(out);
	//std::cout<<.5*(max)/HH<<" "<<1./length<<" "<<delta_f*delta_f<<" "<<delta_f/length<<" "<<length<<" "<<1./(delta_f*length)<<std::endl;
	return .5*(max)/HH;
}

/*! \brief log likelihood function that maximizes over extrinsic parameters tc, phic, D, and phiRef, the reference frequency - for unaligned spins 
 *
 * Ref: arXiv 1603.02444v2
 */
double maximized_Log_Likelihood_unaligned_spin_internal(std::complex<double> *data,
				double *psd,
				double *frequencies,
				std::complex<double> *hplus,
				std::complex<double> *hcross,
				size_t length,
				fftw_outline *plan
				)
{
	double delta_f = frequencies[1]-frequencies[0];
	double *integrand = (double *)malloc(sizeof(double)*length);
	double integral;
	//NOTE: I'm setting detector to HANFORD TEMPORARILY for testing
	//
	//THIS NEEDS TO CHANGE
	//
	//
	//
	//calculate overall template snr squared HH = <H|H> NOT NECESSARY I think
	//std::complex<double> *detector_response = 
	//	(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	//fourier_detector_response(frequencies, length, hplus,hcross, detector_response, 0.0,0.0,"Hanford");

	//for (int i =0;i< length;i++)
	//	integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	//integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double HH = integral;

	//Calculate template snr for plus polarization sqrt(<H+|H+>)
	for (int i =0;i< length;i++){
		integrand[i] = real(hplus[i]*std::conj(hplus[i]))/psd[i];
	}
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HpHproot = sqrt(integral);

	//Calculate template snr for cross polarization sqrt(<Hx|Hx>)
	for (int i =0;i< length;i++)
		integrand[i] = real(hcross[i]*std::conj(hcross[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HcHcroot = sqrt(integral);
	
	//Rescale waveforms from hplus/cross to \hat{hplus/cross} 
	std::complex<double> *hpnorm = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *hcnorm = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	for (int i =0 ;i<length;i++)
	{
		hpnorm[i] = hplus[i]/HpHproot;
		hcnorm[i] = hcross[i]/HcHcroot;
	}

	//calculate \hat{rhoplus/cross} (just denoted rhoplus/cross) <d|hpnorm> 
	//To maximize of coalescence phase, this is an FFT (so its a vector of 
	//<d|h> at discrete tc
	double *rhoplus2 = 
		(double *)malloc(sizeof(double)*length);
	double *rhocross2 = 
		(double *)malloc(sizeof(double)*length);
	std::complex<double> *rhoplus = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *rhocross = 
		(std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	double *gammahat = 
		(double *)malloc(sizeof(double)*length);
	std::complex<double> g_tilde;
	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);	

	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * hpnorm[i] / psd[i]; 
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	for (int i=0;i<length; i++)
	{
		rhoplus[i] = delta_f*(std::complex<double>(out[i][0],out[i][1]));
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhoplus2[i] = delta_f*delta_f*(out[i][0]* out[i][0]+ out[i][1]* out[i][1]);
		
	}
	for (int i=0;i<length; i++)
	{
		g_tilde = 4.*conj(data[i]) * hcnorm[i] / psd[i]; 
		in[i][0] = real(g_tilde);
		in[i][1] = imag(g_tilde);
	}
	//fftw_execute(plan->p);
	fftw_execute_dft(plan->p, in, out);
	for (int i=0;i<length; i++)
	{
		rhocross[i] = delta_f*(std::complex<double>(out[i][0],out[i][1]));
		//Norm of the output, squared (Re{g}^2 + Im{g}^2)
		rhocross2[i] = delta_f*delta_f*(out[i][0]* out[i][0]+ out[i][1]* out[i][1]);
	}
	
	for (int i = 0; i <length;i++)
		gammahat[i] = real(rhoplus[i] * conj(rhocross[i]));
	
	
	for (int i =0;i< length;i++)
		integrand[i] = real(hpnorm[i]*std::conj(hcnorm[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double Ipc = integral;

	double *lambda = (double *)malloc(sizeof(double) * length);
	for(int i = 0; i < length; i++){
		lambda[i] =  (rhoplus2[i] + rhocross2[i] - 2*gammahat[i] * Ipc +
			sqrt( (rhoplus2[i] -rhocross2[i])*(rhoplus2[i] -rhocross2[i]) +
			4. * (Ipc*rhoplus2[i] - gammahat[i] ) * (Ipc*rhocross2[i] - gammahat[i])))
			/(1. - Ipc*Ipc);
	}
	double max = .25 * (*std::max_element(lambda, lambda+length)); 
	//double maxrhoplus2 =  (*std::max_element(rhoplus2, rhoplus2+length)); 
	//double maxrhocross2 =  (*std::max_element(rhocross2, rhocross2+length)); 
	//double maxgammahat =  (*std::max_element(gammahat, gammahat+length)); 

	free(integrand);
	free(hpnorm);
	free(hcnorm);
	free(rhoplus2);
	free(rhoplus);
	free(rhocross);
	free(rhocross2);
	free(gammahat);
	free(lambda);
	fftw_free(in);
	fftw_free(out);

	//std::cout<<max<<" "<<maxrhoplus2 <<" "<<maxrhocross2 <<" "<<maxgammahat << " "<<Ipc<<std::endl;
	//return -0.5*(HH- 2*max);
	return max;
}

/*! \brief Internal function for the unmarginalized log of the likelihood 
 *
 * .5 * ( ( h | h ) - 2 ( D | h ) )
 */
double Log_Likelihood_internal(std::complex<double> *data,
			double *psd,
			double *frequencies,
			std::complex<double> *detector_response,
			int length,
			fftw_outline *plan
			)
{
	
	double delta_f = frequencies[length/2]-frequencies[length/2-1];
	double sum = 0.;
	double *integrand = (double *)malloc(sizeof(double)*length);
	for (int i =0;i< length;i++)
		integrand[i] = real(detector_response[i]*std::conj(detector_response[i]))/psd[i];
	double integral = 4.*simpsons_sum(delta_f, length, integrand);
	double HH = integral;

	for (int i =0;i< length;i++)
		integrand[i] = real(data[i]*std::conj(detector_response[i]))/psd[i];
	integral = 4.*simpsons_sum(delta_f, length, integrand);
	double DH = integral;

	//for (int i =0;i< length;i++)
	//	integrand[i] = real(data[i]*std::conj(data[i]))/psd[i];
	//integral = 4.*simpsons_sum(delta_f, length, integrand);
	//double DD = integral;

	free(integrand);
	//std::cout<<HH<<" "<<DH<<" "<<-0.5*(HH- 2*DH)<<std::endl;
	return -0.5*(HH- 2*DH);
}


struct skysearch_params{
	std::complex<double> *hplus;
	std::complex<double> *hcross;
};

/*! \brief Takes in intrinsic parameters and conducts a sky search using uncorrelated MCMC
 *
 * Searches over RA, sin(DEC), PSI, INCLINATION, phiREF, tc,ln(DL)
 *
 * Only works with PhenomD
 *
 * hplus and hcross should be constructed with the value for inclination in initial_pos -- its scaled out internally
 *
 */
void SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	double *bppe,
	std::complex<double> *hplus,
	std::complex<double> *hcross,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//std::cout<<"MCMC GMST: "<<mcmc_gmst<<std::endl;
	mcmc_Nmod = Nmod;
	mcmc_bppe = bppe;
	mcmc_log_beta = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	//mcmc_save_waveform = true;
	//for(int i =1 ;i<mcmc_num_detectors; i++){
	//	if( mcmc_data_length[i] != mcmc_data_length[0] ||
	//		mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
	//		mcmc_frequencies[i][mcmc_data_length[i]-1] 
	//			!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
	//		mcmc_save_waveform= false;
	//	}
	//		
	//}

	bool local_seeding ;
	if(!seeding_var){
		local_seeding = true;
	}
	else
		local_seeding = false;

	PTMCMC_method_specific_prep("SkySearch", dimension, seeding_var, local_seeding);

	skysearch_params p;
	p.hplus = hplus;
	p.hcross = hcross;

	double cos_i_initial = (initial_pos[3]);
	double p_fac = 0.5*(1+cos_i_initial*cos_i_initial)/(exp(initial_pos[6])*MPC_SEC);
	double c_fac = (cos_i_initial)/(exp(initial_pos[6])*MPC_SEC);
	for(int i = 0 ; i<data_length[0]; i ++){
		p.hplus[i] /= p_fac;
		p.hcross[i] /= c_fac;
	}
	
	skysearch_params **P = new skysearch_params*[chain_N];
	for(int i = 0 ; i<chain_N; i ++){
		P[i]= &p;
	}
	

	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var, chain_temps, 
		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,max_chunk_size,chain_distribution_scheme,
		//log_prior,MCMC_likelihood_wrapper_SKYSEARCH, MCMC_fisher_wrapper_SKYSEARCH,(void **)P,numThreads, pool, 
		log_prior,MCMC_likelihood_wrapper_SKYSEARCH, NULL,(void **)P,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}

double MCMC_likelihood_wrapper_SKYSEARCH(double *param, mcmc_data_interface *interface,void *parameters)
{
	int dimension = interface->max_dim;
	double ll = 0 ; 
	skysearch_params *p = (skysearch_params *) parameters;
	std::complex<double> *hplus = new std::complex<double>[mcmc_data_length[0]];	
	std::complex<double> *hcross = new std::complex<double>[mcmc_data_length[0]];	
	double p_fac = 0.5*(1+param[3]*param[3])/(exp(param[6])*MPC_SEC);
	double c_fac = (param[3])/(exp(param[6])*MPC_SEC);
	for(int i = 0 ; i<mcmc_data_length[0]; i ++){
		hplus[i] = p->hplus[i]*p_fac;
		hcross[i] = p->hcross[i]*c_fac;
	}
	std::complex<double> *response = new std::complex<double>[mcmc_data_length[0]];	
	double delta_t,tc;
	//double tc_ref = param[5],tc;
	double T = 1./( mcmc_frequencies[0][1]-mcmc_frequencies[0][0]);
	double tc_ref = T-param[5];
	//double llvec[3];
	for(int i = 0 ; i<mcmc_num_detectors; i++){
		delta_t = DTOA_DETECTOR(param[0],asin(param[1]),mcmc_gmst, mcmc_detectors[0],mcmc_detectors[i]);
		//tc = tc_ref + delta_t;
		tc = tc_ref - delta_t;
		tc*=2.*M_PI;
		fourier_detector_response_equatorial(mcmc_frequencies[i],mcmc_data_length[i],hplus,hcross, response,param[0],asin(param[1]),param[2],mcmc_gmst, (double *)NULL,0.,0.,0.,0.,mcmc_detectors[i]);
		for(int j = 0 ; j<mcmc_data_length[i];j++){
			//response[j]*=exp(-std::complex<double>(
			//	0,tc*(mcmc_frequencies[i][j]-20.) - param[4] ));
			response[j]*=exp(std::complex<double>(
				0,tc*(mcmc_frequencies[i][j]) - param[4] ));
		}	
		ll += Log_Likelihood_internal(mcmc_data[i], 
			mcmc_noise[i],
			mcmc_frequencies[i],
			response,
			(size_t) mcmc_data_length[i],
			&mcmc_fftw_plans[i]
			);
		//llvec[i]=ll;
		
	}
	//debugger_print(__FILE__,__LINE__,ll);
	//std::cout<<ll<<" "<<llvec[0]<<" "<<llvec[1]-llvec[0]<<" "<<llvec[2]-llvec[1]<<std::endl;
	delete [] response;
	delete [] hplus;
	delete [] hcross;
	return ll;
}


/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void continue_RJPTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int ***status,
	int max_dim,
	int min_dim,
	int N_steps,
	int swp_freq,
	double(*log_prior)(double *param, int *status, mcmc_data_interface *interface,void *parameters),
	void (*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status,mcmc_data_interface *interface, void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	int Nmod,
	double *bppe,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename
	)
{
	std::fstream file_in;
	file_in.open(start_checkpoint_file,std::ios::in);
	std::string line, item;
	int chain_n;
	if(file_in){
		std::getline(file_in,line);
		std::stringstream lineStream(line);
		std::getline(lineStream, item, ',');
		std::getline(lineStream, item, ',');
		std::getline(lineStream, item, ',');
		chain_n = std::stod(item);
		
	}
	file_in.close();
	
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	mcmc_Nmod = Nmod;
	mcmc_bppe = bppe;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//Random numbers for RJstep:
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	mcmc_rvec = (gsl_rng **)malloc(sizeof(gsl_rng *)*chain_n);
	for(int i = 0 ; i<chain_n; i++){
		mcmc_rvec[i] = gsl_rng_alloc(T);
		gsl_rng_set(mcmc_rvec[i],i*11);
	}

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	//##################################################################
	//std::function<void(double*, double*, int*, int*, int, int, double)> RJstep;
	void (*RJstep)(double *, double *, int *, int *, mcmc_data_interface *,void *);
	if(!RJ_proposal){
		RJstep = RJPTMCMC_RJ_proposal;
	}
	else{
		RJstep = RJ_proposal;
	}
	//#################################################################
	bool local_seeding=false ;

	RJPTMCMC_method_specific_prep(generation_method, max_dim,min_dim, NULL, local_seeding);

	continue_RJPTMCMC_MH(start_checkpoint_file,output, status,N_steps,swp_freq,log_prior,
			RJPTMCMC_likelihood_wrapper, RJPTMCMC_fisher_wrapper,RJstep,user_parameters,numThreads, pool, 
			show_prog,statistics_filename,chain_filename,
			auto_corr_filename, likelihood_log_filename,final_checkpoint_filename);
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	for(int i = 0 ; i<chain_n; i++){
		gsl_rng_free(mcmc_rvec[i]);	
	}
	free(mcmc_rvec);
}
/*! \brief Wrapper for the RJPTMCMC_MH function, specifically for GW analysis
 *
 * Handles the details of setting up the MCMC sampler and wraps the fisher and log likelihood to conform to the format of the sampler
 *
 * *NOTE*  -- This sampler as a whole is NOT thread safe. There is global memory declared for each call to MCMC_MH_GW, so separate samplers should not be run in the same process space.
 *
 * Supported parameter combinations:
 *
 * IMRPhenomD - 8 dimensions --  
 * 	MIN DIMENSIONS	-- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, 
 * 	TRANSDIMENSIONAL DIMENSIONS	-- ppE parameters for the bppe array  specified
 *
 * If RJ_proposal is NULL, a default proposal is used.
 */
void RJPTMCMC_MH_GW(double ***output,
		int ***status,
		int max_dim,
		int min_dim,
		int N_steps,
		int chain_N,
		double *initial_pos,
		int *initial_status,
		double *seeding_var,	
		double *chain_temps,
		int swp_freq,
		double(*log_prior)(double *param, int *status, mcmc_data_interface *interface,void *parameters),
		void (*RJ_proposal)(double *current_param, double *proposed_param, int *current_status, int *proposed_status, mcmc_data_interface *interface,void *parameters),
		int numThreads,
		bool pool,
		bool show_prog,
		int num_detectors,
		std::complex<double> **data,
		double **noise_psd,
		double **frequencies,
		int *data_length,
		double gps_time,
		std::string *detectors,
		int Nmod_max,
		double *bppe,
		std::string generation_method,
		std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
		std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
		std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
		std::string likelihood_log_filename,
		std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
					)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	mcmc_Nmod_max = Nmod_max;
	mcmc_bppe = bppe;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;
	mcmc_max_dim = max_dim ;
	mcmc_min_dim = min_dim ;

	//Random numbers for RJstep:
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	mcmc_rvec = (gsl_rng **)malloc(sizeof(gsl_rng *)*chain_N);
	for(int i = 0 ; i<chain_N; i++){
		mcmc_rvec[i] = gsl_rng_alloc(T);
		gsl_rng_set(mcmc_rvec[i],i*11);
	}

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;
	
	//##################################################################
	//std::function<void(double*, double*, int*, int*, int, int, double)> RJstep;
	void (*RJstep)(double *, double *, int *, int *, mcmc_data_interface *,void *);
	if(!RJ_proposal){
		RJstep = RJPTMCMC_RJ_proposal;
	}
	else{
		RJstep = RJ_proposal;
	}
	//##################################################################
	RJPTMCMC_method_specific_prep(generation_method, max_dim, min_dim,seeding_var, local_seeding);
	RJPTMCMC_MH(output, status,max_dim,min_dim, N_steps, chain_N, initial_pos,initial_status,seeding_var, chain_temps, swp_freq,
		 log_prior,RJPTMCMC_likelihood_wrapper, RJPTMCMC_fisher_wrapper, RJstep,user_parameters,numThreads, pool, show_prog,statistics_filename,
		chain_filename,auto_corr_filename,likelihood_log_filename, checkpoint_file);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
	for(int i = 0 ; i<chain_N; i++){
		gsl_rng_free(mcmc_rvec[i]);	
	}
	free(mcmc_rvec);
}

/*! \brief Wrapper for the MCMC_MH function, specifically for GW analysis
 *
 * Handles the details of setting up the MCMC sampler and wraps the fisher and log likelihood to conform to the format of the sampler
 *
 * *NOTE*  -- This sampler is NOT thread safe. There is global memory declared for each call to MCMC_MH_GW, so separate samplers should not be run in the same process space
 *
 * Supported parameter combinations:
 *
 * IMRPhenomD - 4 dimensions -- ln chirpmass, eta, chi1, chi2
 *
 * IMRPhenomD - 7 dimensions -- ln D_L, tc, phic, ln chirpmass, eta, chi1, chi2
 *
 * IMRPhenomD - 9 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, psi
 * 
 * dCS_IMRPhenomD- 8 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, \alpha^2 (the coupling parameter)
 *
 * dCS_IMRPhenomD_root_alpha- 8 dimensions -- cos inclination, RA, DEC, ln D_L, ln chirpmass, eta, chi1, chi2, \sqrt \alpha (in km) (the coupling parameter)
 *
 * IMRPhenomPv2 - 9 dimensions -- cos J_N, ln chirpmass, eta, |chi1|, |chi1|, theta_1, theta_2, phi_1, phi_2
 */
void PTMCMC_MH_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	double *initial_pos,
	double *seeding_var,	
	double *chain_temps,
	int swp_freq,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string likelihood_log_filename,
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
				)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_gmst = 2.456551496456678;
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	//MCMC_modification_struct ms;
	//ms.ppE_Nmod = Nmod;
	//ms.bppe = bppe;
	//mcmc_mod_struct = &ms;
	mcmc_mod_struct=mod_struct;
	
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var, local_seeding);

	PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var, chain_temps, swp_freq,
		 log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, show_prog,statistics_filename,
		 //log_prior,MCMC_likelihood_wrapper, NULL,user_parameters,numThreads, pool, show_prog,statistics_filename,
		chain_filename, checkpoint_file);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}
/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(std::string checkpoint_file_start,
	mcmc_sampler_output *sampler_output,
	double **output,
	int N_steps,
	int max_chain_N_thermo_ensemble,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	int dimension;
	dimension_from_checkpoint_file(checkpoint_file_start,&dimension, &dimension);
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding=false;
	double *seeding_var=NULL;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var, local_seeding);

	continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated(checkpoint_file_start,sampler_output,output,  N_steps,  
		max_chain_N_thermo_ensemble, chain_temps, 
		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,max_chunk_size,chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,user_parameters,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
}
/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int corr_threshold,
	int corr_segments,
	double corr_converge_thresh,
	double corr_target_ac,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	bool GAUSS_QUAD = false;
	std::mutex fisher_mutex;

	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	mcmc_noise = noise_psd;	
	mcmc_init_pos = initial_pos;
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = mod_struct->ppE_Nmod;
	//mcmc_bppe = mod_struct->bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;

	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var, local_seeding);


	//######################################################
	int T = (int)(1./(mcmc_frequencies[0][1]-mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmc_num_detectors];
	double **burn_freqs = new double*[mcmc_num_detectors];
	double **burn_noise = new double*[mcmc_num_detectors];
	int *burn_lengths = new int[mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[num_detectors];
	for(int j = 0; j<mcmc_num_detectors; j++){
		burn_lengths[j] = mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmc_data[j][k];
				burn_freqs[j][ct] = mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	double *fish_freqs = NULL;
	double *fish_weights = NULL;
	double **fish_psd = NULL; //Needs to be interpolated from data given
	int fish_length = 100;
	double flow = mcmc_frequencies[0][0];
	double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	if(GAUSS_QUAD){
		fish_freqs = new double[fish_length];	
		fish_weights = new double[fish_length];	
		fish_psd = new double*[mcmc_num_detectors];	
		
		gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
		for(int i = 0 ; i<fish_length; i++){
			fish_freqs[i]=pow(10,fish_freqs[i]);
		}
		for(int i = 0 ; i<mcmc_num_detectors; i++){
			fish_psd[i] = new double[fish_length];
			gsl_interp_accel *accel = gsl_interp_accel_alloc();
			gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
			gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
			for(int j = 0 ; j<fish_length; j++){
				fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
			}
			gsl_spline_free(spline);
		}
	}
	//######################################################
	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= GAUSS_QUAD;
		user_parameters[i]->fish_freqs= fish_freqs;
		user_parameters[i]->fish_weights= fish_weights;
		user_parameters[i]->fish_psd= fish_psd;
		user_parameters[i]->fish_length= fish_length;

		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################


	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var, chain_temps, 
		swp_freq, t0, nu, corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,max_chunk_size,chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void**)user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,(void **)user_parameters,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	if(GAUSS_QUAD){
		delete [] fish_freqs;
		delete [] fish_weights;
		for(int i = 0 ;i<mcmc_num_detectors; i++){
			delete [] fish_psd[i];
		}
		delete [] fish_psd;
	}
	delete [] burn_data;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	//#################################################
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void PTMCMC_MH_dynamic_PT_alloc_GW(double ***output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var)
		local_seeding = true;
	else
		local_seeding = false;

	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var, local_seeding);

	PTMCMC_MH_dynamic_PT_alloc(output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var, chain_temps, 
		swp_freq, t0, nu, chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, checkpoint_filename);
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}



/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void continue_PTMCMC_MH_GW(std::string start_checkpoint_file,
	double ***output,
	int dimension,
	int N_steps,
	int swp_freq,
	double(*log_prior)(double *param, mcmc_data_interface *interface ,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string auto_corr_filename,
	std::string likelihood_log_filename,
	std::string final_checkpoint_filename,
	bool tune
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}
	void **user_parameters=NULL;
	mcmc_noise = noise_psd;	
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = Nmod;
	//mcmc_bppe = bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmc_num_detectors; i++){
		if( mcmc_data_length[i] != mcmc_data_length[0] ||
			mcmc_frequencies[i][0]!= mcmc_frequencies[0][0] ||
			mcmc_frequencies[i][mcmc_data_length[i]-1] 
				!= mcmc_frequencies[0][mcmc_data_length[0]-1]){
			mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding=false ;

	PTMCMC_method_specific_prep(generation_method, dimension, NULL, local_seeding);

	continue_PTMCMC_MH(start_checkpoint_file,output, N_steps,swp_freq,log_prior,
			MCMC_likelihood_wrapper, MCMC_fisher_wrapper,user_parameters,numThreads, pool, 
			//MCMC_likelihood_wrapper, NULL,numThreads, pool, 
			show_prog,statistics_filename,chain_filename,
			final_checkpoint_filename,tune);
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	free(plans);
}
/*! \brief Unpacks MCMC parameters for method specific initiation (RJ version)
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void RJPTMCMC_method_specific_prep(std::string generation_method, int max_dim, int min_dim,double *seeding_var, bool local_seeding)
{
	if(min_dim==9 && (generation_method =="ppE_IMRPhenomD_Inspiral"|| generation_method =="ppE_IMRPhenomD_IMR")){
		mcmc_intrinsic=false;
		std::cout<<"Sampling in parameters: cos inclination, RA, DEC, ln DL ,ln chirpmass, eta, chi1, chi2, psi";
		for(int i =0; i<mcmc_Nmod_max; i++){
			std::cout<<", beta"<<i<<" ("<<mcmc_bppe[i]<<")";
		}
		std::cout<<endl;
		mcmc_log_beta = false;
		if(local_seeding){
			seeding_var = new double[max_dim];
			seeding_var[0]=.1;
			seeding_var[1]=.5;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			seeding_var[4]=.1;
			seeding_var[5]=.1;
			seeding_var[6]=.1;
			seeding_var[7]=.1;
			seeding_var[8]=.1;
			for(int i =0; i< mcmc_Nmod_max;i++){
				seeding_var[9 + i]=2;
			}
		}
	}

}

/*! \brief Unpacks MCMC parameters for method specific initiation 
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void PTMCMC_method_specific_prep(std::string generation_method, int dimension,double *seeding_var, bool local_seeding)
{
	if(dimension==4 && generation_method =="IMRPhenomD"){
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=.5;
			seeding_var[1]=.1;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
		}
		mcmc_intrinsic=true;
	}
	else if(generation_method =="SkySearch"){
		std::cout<<"Sampling in parameters: "<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1;
			seeding_var[1]=.3;
			seeding_var[2]=1.;
			seeding_var[3]=1.;
			seeding_var[4]=1.;
			seeding_var[5]=1.;
			seeding_var[6]=1.;
		}
		mcmc_intrinsic=false;
	}
	else if(dimension==11 && generation_method =="IMRPhenomD"){
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
		mcmc_intrinsic=false;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=.1;
			seeding_var[1]=.1;
			seeding_var[2]=.1;
			seeding_var[3]=.1;
			seeding_var[4]=.5;
			seeding_var[5]=.1;
			seeding_var[6]=.1;
			seeding_var[7]=.1;
			seeding_var[8]=.1;
			seeding_var[9]=.1;
			seeding_var[10]=.5;
		}
	}
	else if(dimension==12 && generation_method =="dCS_IMRPhenomD"){
		mcmc_mod_struct->ppE_Nmod = 1;
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, sqrt(alpha) (km)"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.5;
			seeding_var[2]=1;
			seeding_var[3]=1;
			seeding_var[4]=1;
			seeding_var[5]=.1;
			seeding_var[6]=1;
			seeding_var[7]=1.;
			seeding_var[8]=.3;
			seeding_var[9]=.5;
			seeding_var[10]=.5;
			seeding_var[11]=1;
		}
	}
	//Sampling parameter is root alpha in km
	else if(dimension==12 && generation_method =="EdGB_IMRPhenomD"){
		mcmc_mod_struct->ppE_Nmod = 1;
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, sqrt(alpha) (km)"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.5;
			seeding_var[2]=1;
			seeding_var[3]=1;
			seeding_var[4]=1;
			seeding_var[5]=.1;
			seeding_var[6]=1;
			seeding_var[7]=1.;
			seeding_var[8]=.3;
			seeding_var[9]=.5;
			seeding_var[10]=.5;
			seeding_var[11]=1;
		}
	}
	else if(dimension==16 && generation_method =="dCS_IMRPhenomPv2"){
		mcmc_mod_struct->ppE_Nmod = 1;
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2, sqrt(alpha) (km)"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.5;
			seeding_var[2]=1;
			seeding_var[3]=1;
			seeding_var[4]=1;
			seeding_var[5]=.1;
			seeding_var[6]=1;
			seeding_var[7]=1.;
			seeding_var[8]=.3;
			seeding_var[9]=.5;
			seeding_var[10]=.5;
			seeding_var[11]=1.;
			seeding_var[12]=1.;
			seeding_var[13]=1.;
			seeding_var[14]=1.;
			seeding_var[15]=1.;
		}
	}
	//Sampling parameter is root alpha in km
	else if(dimension==16 && generation_method =="EdGB_IMRPhenomPv2"){
		mcmc_mod_struct->ppE_Nmod = 1;
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2, sqrt(alpha) (km)"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.5;
			seeding_var[2]=1;
			seeding_var[3]=1;
			seeding_var[4]=1;
			seeding_var[5]=.1;
			seeding_var[6]=1;
			seeding_var[7]=1.;
			seeding_var[8]=.3;
			seeding_var[9]=.5;
			seeding_var[10]=.5;
			seeding_var[11]=1.;
			seeding_var[12]=1.;
			seeding_var[13]=1.;
			seeding_var[14]=1.;
			seeding_var[15]=1.;
		}
	}
	else if(dimension==15 && generation_method =="IMRPhenomPv2"){
		mcmc_intrinsic=false;
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.5;
			seeding_var[2]=1;
			seeding_var[3]=1;
			seeding_var[4]=1;
			seeding_var[5]=.1;
			seeding_var[6]=1;
			seeding_var[7]=1.;
			seeding_var[8]=.3;
			seeding_var[9]=.5;
			seeding_var[10]=.5;
			seeding_var[11]=1.;
			seeding_var[12]=1.;
			seeding_var[13]=1.;
			seeding_var[14]=1.;
		}
	}
	else if(dimension==8 && generation_method =="IMRPhenomPv2"){
		mcmc_intrinsic=true;
		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		if(local_seeding){
			seeding_var = new double[dimension];
			seeding_var[0]=1.;
			seeding_var[1]=.3;
			seeding_var[2]=1.;
			seeding_var[3]=1.;
			seeding_var[4]=1.;
			seeding_var[5]=1.;
			seeding_var[6]=1.;
			seeding_var[7]=1.;
		}
	}
	else if( generation_method == "ppE_IMRPhenomPv2_Inspiral"
		|| generation_method == "ppE_IMRPhenomPv2_IMR"
		){
		if(dimension-mcmc_mod_struct->ppE_Nmod == 15){
			
		std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
			for(int i =0; i<mcmc_Nmod; i++){
				std::cout<<", beta"<<i;
			}
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.5;
				seeding_var[2]=1;
				seeding_var[3]=1;
				seeding_var[4]=1;
				seeding_var[5]=.1;
				seeding_var[6]=1;
				seeding_var[7]=1.;
				seeding_var[8]=.3;
				seeding_var[9]=.5;
				seeding_var[10]=.5;
				seeding_var[11]=1.;
				seeding_var[12]=1.;
				seeding_var[13]=1.;
				seeding_var[14]=1.;
				for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
					seeding_var[15 + i]=2;
				}
			}
		}
		if(dimension-mcmc_mod_struct->ppE_Nmod == 8){
			mcmc_intrinsic=true;
			std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
			for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
				std::cout<<", beta"<<i;
			}
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.3;
				seeding_var[2]=1.;
				seeding_var[3]=1.;
				seeding_var[4]=1.;
				seeding_var[5]=1.;
				seeding_var[6]=1.;
				seeding_var[7]=1.;
				for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
					seeding_var[8 + i]=2;
				}
			}
		}

	}
	else if(generation_method == "gIMRPhenomPv2" 
		){
		int totalmod = (mcmc_mod_struct->gIMR_Nmod_phi + mcmc_mod_struct->gIMR_Nmod_sigma + mcmc_mod_struct->gIMR_Nmod_beta + mcmc_mod_struct->gIMR_Nmod_alpha );
		if(dimension-totalmod == 15){
			
			mcmc_intrinsic=false;
			std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
			for(int i =0; i<totalmod; i++){
				std::cout<<", mod "<<i;
			}
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.5;
				seeding_var[2]=1;
				seeding_var[3]=1;
				seeding_var[4]=1;
				seeding_var[5]=.1;
				seeding_var[6]=1;
				seeding_var[7]=1.;
				seeding_var[8]=.3;
				seeding_var[9]=.5;
				seeding_var[10]=.5;
				seeding_var[11]=1.;
				seeding_var[12]=1.;
				seeding_var[13]=1.;
				seeding_var[14]=1.;
				for(int i =0; i< totalmod;i++){
					seeding_var[15 + i]=2;
				}
			}
		}
		else{
			mcmc_intrinsic=true;
			std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2";
			for(int i =0; i<totalmod; i++){
				std::cout<<", mod "<<i;
			}
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.3;
				seeding_var[2]=1.;
				seeding_var[3]=1.;
				seeding_var[4]=1.;
				seeding_var[5]=1.;
				seeding_var[6]=1.;
				seeding_var[7]=1.;
				for(int i =0; i< totalmod;i++){
					seeding_var[8 + i]=2;
				}
			}
		}

	}
	else if(generation_method == "gIMRPhenomD" 
		){
		int totalmod = (mcmc_mod_struct->gIMR_Nmod_phi + mcmc_mod_struct->gIMR_Nmod_sigma + mcmc_mod_struct->gIMR_Nmod_beta + mcmc_mod_struct->gIMR_Nmod_alpha );
		if(dimension-totalmod == 4){
			std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
			
			mcmc_intrinsic=true;
			for(int i =0; i<totalmod; i++){
				std::cout<<", beta"<<i;
			}
			
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.3;
				seeding_var[2]=1;
				seeding_var[3]=1;
				for(int i =0; i< totalmod;i++){
					seeding_var[4 + i]=2;
				}
			}
		}
		else{
			std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2";
			mcmc_intrinsic=false;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.5;
				seeding_var[2]=1;
				seeding_var[3]=1;
				seeding_var[4]=1;
				seeding_var[5]=.1;
				seeding_var[6]=1;
				seeding_var[7]=1.;
				seeding_var[8]=.3;
				seeding_var[9]=.5;
				seeding_var[10]=.5;
				for(int i =0; i< totalmod;i++){
					seeding_var[11 + i]=2;
				}
			}
		}

	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral"
		|| generation_method == "ppE_IMRPhenomD_IMR"
		){
		if(dimension-mcmc_mod_struct->ppE_Nmod == 11){
			
			std::cout<<"Sampling in parameters: RA, DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2";
			for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
				std::cout<<", beta"<<i;
			}
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=1.;
				seeding_var[1]=.5;
				seeding_var[2]=1;
				seeding_var[3]=1;
				seeding_var[4]=1;
				seeding_var[5]=.1;
				seeding_var[6]=1;
				seeding_var[7]=1.;
				seeding_var[8]=.3;
				seeding_var[9]=.5;
				seeding_var[10]=.5;
				for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
					seeding_var[11 + i]=2;
				}
			}
		}
		if(dimension-mcmc_mod_struct->ppE_Nmod == 4){
			mcmc_intrinsic=true;
			std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
			for(int i =0; i<mcmc_mod_struct->ppE_Nmod; i++){
				std::cout<<", beta"<<i;
			}
			std::cout<<endl;
			if(local_seeding){
				seeding_var = new double[dimension];
				seeding_var[0]=.1;
				seeding_var[1]=.1;
				seeding_var[2]=.1;
				seeding_var[3]=1;
				seeding_var[4]=.5;
				seeding_var[5]=.1;
				seeding_var[6]=.1;
				seeding_var[7]=.1;
				for(int i =0; i< mcmc_mod_struct->ppE_Nmod;i++){
					seeding_var[8 + i]=2;
				}
			}
		}

	}
	else{
		std::cout<<
			"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
		exit(1);
	}
}


void MCMC_fisher_wrapper(double *param,  double **output, mcmc_data_interface *interface,void *parameters)
{
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;
	int dimension = interface->max_dim;
	double *temp_params = new double[dimension];
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&params, dimension, mcmc_generation_method,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################
	repack_parameters(param, &params, 
		"MCMC_"+mcmc_generation_method, dimension, NULL);
	//repack_parameters(mcmc_init_pos, &params, 
	//	"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//std::cout<<"INCL angle fisher: "<<params.incl_angle<<std::endl;
	for(int j =0; j<dimension; j++){
		for(int k =0; k<dimension; k++)
		{
			output[j][k] =0;
		}
	} 
	double **temp_out = allocate_2D_array(dimension,dimension);
	for(int i =0 ; i <mcmc_num_detectors; i++){
		
		//Use AD 
		if(user_param->GAUSS_QUAD)
		{	
			std::unique_lock<std::mutex> lock{*(user_param->mFish)};
			//fisher_autodiff(mcmc_frequencies[i], mcmc_data_length[i],
			//	"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmc_noise[i]);
			fisher_autodiff(user_param->fish_freqs, user_param->fish_length,
				"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
				(gen_params *)(&params),  "GAUSSLEG",user_param->fish_weights,true,user_param->fish_psd[i]);
		}
		else{
			fisher_numerical(mcmc_frequencies[i], mcmc_data_length[i],
				"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
				&params, mcmc_deriv_order, NULL, NULL, mcmc_noise[i]);

		}
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++)
			{
				output[j][k] +=temp_out[j][k];
				//if(std::isnan(output[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmc_generation_method.find("Pv2") && !mcmc_intrinsic){
	if(!mcmc_intrinsic){
		output[0][0] += 1./(4*M_PI*M_PI);//RA
		output[1][1] += 1./4;//sin DEC
		output[2][2] += 1./(4*M_PI*M_PI);//psi
		output[3][3] += 1./(4);//cos iota
		output[4][4] += 1./(4*M_PI*M_PI);//phiref
		output[5][5] += 1./(.01);//tc
		output[8][8] += 1./.25;//eta
		output[9][9] += 1./4;//spin1
		output[10][10] += 1./4;//spin2
		if(mcmc_generation_method.find("PhenomPv2") != std::string::npos){
			output[11][11] += 1./4;//cos theta1
			output[12][12] += 1./4;//cos theta2
			output[13][13] += 1./(4*M_PI*M_PI);//phi1
			output[14][14] += 1./(4*M_PI*M_PI);//phi2
		}
	}
	else{
		if(mcmc_generation_method.find("PhenomPv2") != std::string::npos){
			output[1][1] =1./(.25) ;//eta
			output[2][2] =1./(4);//spin1
			output[3][3] =1./(4);//spin2
			output[4][4] =1./(4);//cos theta1
			output[5][5] =1./(4);//cos theta2
			output[6][6] =1./(4*M_PI*M_PI) ;//phi1
			output[7][7] =1./(4*M_PI*M_PI) ;//phi2
		}
		else if (mcmc_generation_method.find("PhenomD")!=std::string::npos){
			output[1][1] =1./(.25) ;//eta
			output[2][2] =1./(4) ;//spin1
			output[3][3] =1./(4) ;//spin2
	
		}
	}
	deallocate_2D_array(temp_out, dimension,dimension);
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		if(local_gen.find("ppE") != std::string::npos ||
			local_gen.find("dCS")!=std::string::npos||
			local_gen.find("EdGB")!=std::string::npos){
			delete [] params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}


//RA, DEC, and PSI were absorbed into gen_params structure -- remove from arguments
double MCMC_likelihood_extrinsic(bool save_waveform, gen_params_base<double> *parameters,std::string generation_method, int *data_length, double **frequencies, std::complex<double> **data, double **psd, std::string *detectors, fftw_outline *fftw_plans, int num_detectors, double RA, double DEC,double gps_time)
{
	double *phi = new double[num_detectors];
	double *theta = new double[num_detectors];
	//celestial_horizon_transform(RA,DEC, gps_time, detectors[0], &phi[0], &theta[0]);
	double tc_ref, phic_ref, ll=0, delta_t;
	double LISA_alpha0,LISA_phi0, LISA_thetal, LISA_phil;
	double *times=NULL;
	//Needs some work
	//#######################################################################3
	//if (save_waveform){
	if (false){
		std::complex<double> *hplus = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		std::complex<double> *hcross = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		std::complex<double> *response = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		parameters->tc=0;
		fourier_waveform(frequencies[0], data_length[0], 
			hplus, hcross,generation_method, parameters);
		//std::cout<<hplus[100]<<" "<<hcross[100]<<std::endl;	
		//fourier_detector_response(frequencies[0], data_length[0], 
		//	hplus, hcross, response, parameters->theta, parameters->phi, 
		//	detectors[0]);
		fourier_detector_response_equatorial(frequencies[0], data_length[0], 
			hplus, hcross, response, parameters->RA, parameters->DEC, parameters->psi,
			parameters->gmst,times, LISA_alpha0, LISA_phi0, LISA_thetal, LISA_phil,detectors[0]);
		ll += maximized_coal_Log_Likelihood_internal(data[0], 
				psd[0],
				frequencies[0],
				response,
				(size_t) data_length[0],
				&fftw_plans[0],
				&tc_ref,
				&phic_ref
				);
		
		//Use maximum LL phic if using PhenomD, but if using Pv2
		//PhiRef is already a parameter
		if(generation_method.find("IMRPhenomD")==std::string::npos){
			phic_ref=0;
		}
		//if(generation_method.find("IMRPhenomD")!=std::string::npos){
		//	parameters->phiRef=phic_ref;
		//}
		for(int i=1; i < num_detectors; i++){
			//celestial_horizon_transform(RA,DEC, gps_time, 
			//		mcmc_detectors[i], &phi[i], &theta[i]);
			//parameters->phi=phi[i];
			//parameters->theta=theta[i];
			//delta_t = DTOA(theta[0], theta[i], detectors[0], detectors[i]);
			delta_t = DTOA_DETECTOR(parameters->RA, parameters->DEC,mcmc_gmst, detectors[0], detectors[i]);
			//parameters->tc = tc_ref + delta_t;
			parameters->tc = tc_ref - delta_t;
			
			//fourier_detector_response(frequencies[i], 
			//	data_length[i], hplus, hcross, response, 
			//	parameters->theta, parameters->phi, parameters->psi,
			//	detectors[i]);
			fourier_detector_response_equatorial(frequencies[i], data_length[i], 
				hplus, hcross, response, parameters->RA, parameters->DEC, parameters->psi,
				parameters->gmst,times, LISA_alpha0, LISA_phi0, LISA_thetal, LISA_phil,detectors[i]);
			for(int j =0; j<data_length[i]; j++){
				//response[j] *=std::exp(std::complex<double>(0,-parameters->tc*2*M_PI*frequencies[i][j]+ phic_ref) );	
				//response[j] *=std::exp(std::complex<double>(0,-parameters->tc*2*M_PI*frequencies[i][j]) + phic_ref);	
				//response[j] *=std::exp(std::complex<double>(0,2*M_PI*parameters->tc*frequencies[i][j]) + phic_ref);	
				response[j] *=std::exp(std::complex<double>(0,parameters->tc*frequencies[i][j]) + phic_ref);	
			}
			ll += Log_Likelihood_internal(data[i], 
					psd[i],
					frequencies[i],
					response,
					(size_t) data_length[i],
					&fftw_plans[i]
					);
		}
		free(hplus); free(hcross); free(response);
	}
	//#######################################################################3
	//Generally, the data lengths don't have to be the same
	//if(generation_method.find("IMRPhenomPv2") !=std::string::npos){
	//double snr = 0;
	//if(false)
	{
		std::complex<double> *hplus = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		std::complex<double> *hcross = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		std::complex<double> *response = 
			(std::complex<double> *)malloc(sizeof(std::complex<double>)*
				data_length[0]);
		//tc_ref = parameters->tc;
		double T = 1./( frequencies[1]-frequencies[0]);
		tc_ref = T-parameters->tc;
		double tc = tc_ref;
		tc*=2.*M_PI;
		parameters->tc=0;
		fourier_waveform(frequencies[0], data_length[0], 
			hplus, hcross,generation_method, parameters);
		fourier_detector_response_equatorial(frequencies[0], data_length[0], 
			hplus, hcross, response, parameters->RA, parameters->DEC, parameters->psi,
			parameters->gmst,times, LISA_alpha0, LISA_phi0, LISA_thetal, LISA_phil,detectors[0]);
		for(int i = 0 ; i<data_length[0];i++){
			//response[i]*=exp(-std::complex<double>(0,tc*(frequencies[0][i]-parameters->f_ref)));
			//response[i]*=exp(std::complex<double>(0,tc*(frequencies[0][i]-parameters->f_ref)));
			response[i]*=exp(std::complex<double>(0,tc*(frequencies[0][i])));
		}	
		//Referecne detector first
		ll += Log_Likelihood_internal(data[0], 
				psd[0],
				frequencies[0],
				response,
				(size_t) data_length[0],
				&fftw_plans[0]
				);
		//snr+=pow_int(data_snr(frequencies[0],data_length[0],data[0],response,psd[0]),2);
		for(int i=1; i < num_detectors; i++){
			//celestial_horizon_transform(RA,DEC, gps_time, 
			//		detectors[i], &phi[i], &theta[i]);
			//parameters->phi=phi[i];
			//parameters->theta=theta[i];
			//delta_t = DTOA(theta[0], theta[i], detectors[0], detectors[i]);
			delta_t = DTOA_DETECTOR(parameters->RA, parameters->DEC,mcmc_gmst, detectors[0], detectors[i]);
			//tc = tc_ref + delta_t;
			tc = tc_ref - delta_t;
			fourier_detector_response_equatorial(frequencies[i], data_length[i], 
				hplus, hcross, response, parameters->RA, parameters->DEC, 
				parameters->psi,parameters->gmst,times, LISA_alpha0, 
				LISA_phi0, LISA_thetal, LISA_phil,detectors[i]);
			tc*=2.*M_PI;
			for(int j = 0 ; j<data_length[i];j++){
				//response[j]*=exp(-std::complex<double>(
				//	0,tc*(frequencies[i][j]-parameters->f_ref)
				//	));
				//response[j]*=exp(std::complex<double>(
				//	0,tc*(frequencies[i][j]-parameters->f_ref)
				//	));
				response[j]*=exp(std::complex<double>(
					0,tc*(frequencies[i][j])
					));
			}	
			//snr+=pow_int(data_snr(frequencies[i],data_length[i],data[i],response,psd[i]),2);
		
			ll += Log_Likelihood_internal(data[i], 
				psd[i],
				frequencies[i],
				response,
				(size_t) data_length[i],
				&fftw_plans[i]
				);
		}
		free( hplus);
		free( hcross);
		free( response);
	}
	delete [] phi;
	delete [] theta;
	return ll;
	//return 2;
}
/*! \brief utility to do MCMC specific transformations on the input param vector before passing to the repacking utillity
 *
 * Returns the local generation method to be used in the LL functions
 */
std::string MCMC_prep_params(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method, MCMC_modification_struct *mod_struct)
{
	if(mcmc_intrinsic) gen_params->sky_average = true;
	else gen_params->sky_average = false;
	gen_params->f_ref = 20;
	gen_params->shift_time = true;
	gen_params->shift_phase = true;
	//gen_params->shift_time = false;
	//gen_params->shift_phase = false;
	gen_params->gmst = mcmc_gmst;
	gen_params->equatorial_orientation=false;
	gen_params->horizon_coord=false;

	gen_params->NSflag1 = mod_struct->NSflag1;
	gen_params->NSflag2 = mod_struct->NSflag2;
	//gen_params->NSflag1 = false;
	//gen_params->NSflag2 = false;
	if(check_mod(generation_method)){
		if(generation_method.find("ppE") != std::string::npos ||
			generation_method.find("dCS") !=std::string::npos||
			generation_method.find("EdGB") != std::string::npos){
			gen_params->bppe=mcmc_mod_struct->bppe;
			gen_params->Nmod=mcmc_mod_struct->ppE_Nmod;
			gen_params->betappe=new double[gen_params->Nmod];
		}
		else if(generation_method.find("gIMR") != std::string::npos){
			gen_params->Nmod_phi=mcmc_mod_struct->gIMR_Nmod_phi;
			gen_params->phii=mcmc_mod_struct->gIMR_phii;
			if(gen_params->Nmod_phi !=0){
				gen_params->delta_phi=new double[gen_params->Nmod_phi];
			}
			gen_params->Nmod_sigma=mcmc_mod_struct->gIMR_Nmod_sigma;
			gen_params->sigmai=mcmc_mod_struct->gIMR_sigmai;
			if(gen_params->Nmod_sigma !=0){
				gen_params->delta_sigma=new double[gen_params->Nmod_sigma];
			}
			gen_params->Nmod_beta=mcmc_mod_struct->gIMR_Nmod_beta;
			gen_params->betai=mcmc_mod_struct->gIMR_betai;
			if(gen_params->Nmod_beta !=0){
				gen_params->delta_beta=new double[gen_params->Nmod_beta];
			}
			gen_params->Nmod_alpha=mcmc_mod_struct->gIMR_Nmod_alpha;
			gen_params->alphai=mcmc_mod_struct->gIMR_alphai;
			if(gen_params->Nmod_alpha !=0){
				gen_params->delta_alpha=new double[gen_params->Nmod_alpha];
			}
		}
	}
	for(int i = 0 ; i <dimension; i++){
		temp_params[i]=param[i];
	}
	return generation_method;
}
double MCMC_likelihood_wrapper(double *param, mcmc_data_interface *interface ,void *parameters)
{
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;

	int dimension = interface->max_dim;
	double ll = 0;
	double *temp_params = new double[dimension];
	//#########################################################################
	gen_params_base<double> gen_params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&gen_params, dimension, mcmc_generation_method,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################

	//repack_non_parameters(temp_params, &gen_params, 
		//"MCMC_"+mcmc_generation_method, dimension, NULL);
	repack_parameters(temp_params, &gen_params, 
		"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//return 1;
	std::complex<double> **local_data = mcmc_data;
	double **local_freqs = mcmc_frequencies;
	double **local_noise = mcmc_noise;
	int *local_lengths = mcmc_data_length;
	fftw_outline *local_plans = mcmc_fftw_plans;
	//if(interface->burn_phase && user_param->burn_data){
	if(false){
		local_data = user_param->burn_data;
		local_freqs = user_param->burn_freqs;
		local_noise = user_param->burn_noise;
		local_lengths = user_param->burn_lengths;
		local_plans = user_param->burn_plans;
	}
	if(mcmc_intrinsic){
		if(mcmc_generation_method.find("IMRPhenomD") != std::string::npos){
			if(!mcmc_save_waveform){
				for(int i=0; i < mcmc_num_detectors; i++){
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 10;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					std::complex<double> *response =
						(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
					fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmc_detectors[i], local_gen, &gen_params);
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
					//ll += maximized_Log_Likelihood(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		(size_t) mcmc_data_length[i],
					//		&gen_params,
					//		mcmc_detectors[i],
					//		local_gen,
					//		&mcmc_fftw_plans[i]
					//		);
					free(response);
				}
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 10;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				std::complex<double> *response =
					(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
				fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmc_detectors[0], local_gen, &gen_params);
				//std::complex<double> *hc =
				//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
				//std::complex<double> *hp =
				//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
				//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
					//ll += maximized_Log_Likelihood(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		(size_t) mcmc_data_length[i],
					//		&gen_params,
					//		mcmc_detectors[i],
					//		local_gen,
					//		&mcmc_fftw_plans[i]
					//		);
					//ll += maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		hp,
					//		hc,
					//		(size_t) mcmc_data_length[i],
					//		&mcmc_fftw_plans[i]
					//		);
				}
				free(response);
				//free(hp);
				//free(hc);

			}

		}
		else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
			//if(!mcmc_save_waveform){
			if(false){
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 20;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				std::complex<double> *hc =
					(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
				std::complex<double> *hp =
					(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
				fourier_waveform(local_freqs[0],local_lengths[0], hp,hc, local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							hp,
							hc,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
				}
				free(hp);
				free(hc);
				//std::cout<<ll<<std::endl;
			}

		}
	}
	else{
		double RA = gen_params.RA;
		double DEC = gen_params.DEC;
		double PSI = gen_params.psi;
		//if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){
		
		ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, 
			&gen_params,local_gen, local_lengths, 
			local_freqs, local_data, local_noise, mcmc_detectors, 
			local_plans, mcmc_num_detectors, RA, DEC,mcmc_gps_time);
		//ll=2;

		//ll = Log_Likelihood(mcmc_data[0], 
		//		mcmc_noise[0],
		//		mcmc_frequencies[0],
		//		mcmc_data_length[0],
		//		&gen_params,
		//		mcmc_detectors[0],
		//		local_gen,
		//		&mcmc_fftw_plans[0]
		//		);

		//}
		//else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){

		//}
	}
	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		if( local_gen.find("ppE") != std::string::npos ||
			local_gen.find("dCS") != std::string::npos ||
			local_gen.find("EdGB") != std::string::npos){
			delete [] gen_params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] gen_params.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] gen_params.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] gen_params.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] gen_params.delta_alpha;
			}

		}
	}
	return ll;

}
//#########################################################
double RJPTMCMC_likelihood_wrapper(double *param, 
	int *status,
	mcmc_data_interface *interface,
	void *parameters
	) 
{
	int dimension = interface->max_dim;
	double ll = 0 ;
	if(!mcmc_intrinsic && (mcmc_generation_method == "ppE_IMRPhenomD_Inspiral"||mcmc_generation_method == "ppE_IMRPhenomD_IMR")){
	//if(false){

		std::string local_method ;
		if(mcmc_generation_method == "ppE_IMRPhenomD_Inspiral"){
			local_method = "ppE_IMRPhenomD_Inspiral";
		}
		else if(mcmc_generation_method == "ppE_IMRPhenomD_IMR"){
			local_method = "ppE_IMRPhenomD_IMR";
		}
		//unpack parameter vector
		double incl = acos(param[0]);
		double RA = param[1];
		double DEC = param[2];
		double DL = std::exp(param[3]);
		double chirpmass = std::exp(param[4]);
		double eta = param[5];
		double chi1 = param[6];
		double chi2 = param[7];
		double psi = param[8];
		//double lnalpha2 = param[8];
		int mods_ct = 0;
		double *local_bppe = new double[mcmc_Nmod_max];
		double local_beta[mcmc_Nmod_max] ;
		for(int i = mcmc_min_dim ; i < mcmc_max_dim; i++){
			if(status[i] == 1){
				local_bppe[mods_ct] = mcmc_bppe[i-mcmc_min_dim];
				local_beta[mods_ct] = param[i];
				mods_ct++;
			}			
		}
		if(mods_ct ==0){
			local_method="IMRPhenomD";
		}

		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = DL;
		//The rest is maximized over for this option
		parameters.tc = 0;
		parameters.phiRef = 0;
		parameters.incl_angle = incl;
		//parameters.phi=phi[0];
		//parameters.theta=theta[0];
		parameters.gmst=mcmc_gmst;
		parameters.RA=RA;
		parameters.DEC=DEC;
		parameters.NSflag1 = false;
		parameters.NSflag2 = false;
		parameters.sky_average = false;
		parameters.psi = psi;
		if(mods_ct !=0){
			parameters.Nmod =mods_ct;
			parameters.betappe = new double[mods_ct];
			parameters.bppe =new double[mods_ct] ;
			//parameters.betappe[0] = lnalpha2;
			for (int j = 0 ; j<mods_ct; j++){
				parameters.betappe[j] = local_beta[j];
				parameters.bppe[j] = local_bppe[j];
			}
		}
		ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, &parameters,local_method, mcmc_data_length, mcmc_frequencies, mcmc_data, mcmc_noise, mcmc_detectors, mcmc_fftw_plans, mcmc_num_detectors, RA, DEC,mcmc_gps_time);
		if(mods_ct !=0){
			delete [] parameters.betappe;
			delete [] parameters.bppe;
		}
		delete [] local_bppe;
	}
	else if(!mcmc_intrinsic && (mcmc_generation_method == "ppE_IMRPhenomPv2_Inspiral"||mcmc_generation_method == "ppE_IMRPhenomPv2_IMR")){
		std::string local_method ;
		if(mcmc_generation_method == "ppE_IMRPhenomPv2_Inspiral"){
			local_method = "ppE_IMRPhenomPv2_Inspiral";
		}
		else if(mcmc_generation_method == "ppE_IMRPhenomPv2_IMR"){
			local_method = "ppE_IMRPhenomPv2_IMR";
		}
		//else if(false){	
		//unpack parameter vector
		//All parameters defined at 20Hz
		double incl = acos(param[0]);
		double RA = param[1];
		double DEC = param[2];
		double DL = std::exp(param[3]);
		double chirpmass = std::exp(param[4]);
		double eta = param[5];
		double chi1 = param[6];
		double chi2 = param[7];
		double theta1 = param[8];//Polar angles of spins relative to L
		double theta2 = param[9];
		double phi1 = param[10];//Azimuthal angles of spins relative to L
		double phi2 = param[11];
		double phiref = param[12];//Orbital phase at fref (20Hz) -- Someday, find a way to maximize this out
		double psi = param[13];//Polarization angle
		double delta_t = 0;
		double tc_ref =0;
		double fref = 20;
		//transform
		double spin1[3];
		double spin2[3];
		double spin1sphr[3] =  {chi1, theta1, phi1};
		double spin2sphr[3] =  {chi2, theta2, phi2};
		transform_sph_cart( &spin1sphr[0], &spin1[0]);
		transform_sph_cart( &spin2sphr[0], &spin2[0]);
	
		int mods_ct = 0;
		double *local_bppe = new double[mcmc_Nmod_max];
		double local_beta[mcmc_Nmod_max] ;
		for(int i = mcmc_min_dim ; i < mcmc_max_dim; i++){
			if(status[i] == 1){
				local_bppe[mods_ct] = mcmc_bppe[i-mcmc_min_dim];
				local_beta[mods_ct] = param[i];
				mods_ct++;
			}			
		}
		if(mods_ct ==0){
			local_method="IMRPhenomPv2";
		}

		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = spin1[0];
		parameters.spin1[1] = spin1[1];
		parameters.spin1[2] = spin1[2];
		parameters.spin2[0] = spin2[0];
		parameters.spin2[1] = spin2[1];
		parameters.spin2[2] = spin2[2];
		parameters.Luminosity_Distance = DL;
		//The rest is maximized over for this option
		parameters.tc = 0;
		//parameters.phiRef = 0;
		parameters.phiRef=phiref;
		parameters.f_ref=fref;
		parameters.incl_angle = incl;
		//parameters.phi=phi[0];
		//parameters.theta=theta[0];
		parameters.phi=0;
		parameters.theta=0;
		parameters.RA=RA;
		parameters.DEC=DEC;
		parameters.gmst=mcmc_gmst;
		parameters.NSflag1 = false;
		parameters.NSflag2 = false;
		parameters.sky_average = false;
		if(mods_ct !=0){
			parameters.Nmod =mods_ct;
			parameters.betappe = new double[mods_ct];
			parameters.bppe =new double[mods_ct] ;
			//parameters.betappe[0] = lnalpha2;
			for (int j = 0 ; j<mods_ct; j++){
				parameters.betappe[j] = local_beta[j];
				parameters.bppe[j] = local_bppe[j];
			}
		}
		ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, &parameters,local_method, mcmc_data_length, mcmc_frequencies, mcmc_data, mcmc_noise, mcmc_detectors, mcmc_fftw_plans, mcmc_num_detectors, RA, DEC,mcmc_gps_time);
		if(mods_ct !=0){
			delete [] parameters.betappe;
			delete [] parameters.bppe;
		}
		delete [] local_bppe;
	}
	//std::cout<<ll<<std::endl;
	return ll;

}
void RJPTMCMC_RJ_proposal(double *current_params, 
	double *proposed_params, 
	int *current_status, 
	int *proposed_status,
	mcmc_data_interface *interface,
	void *parameters
	) 
{
	int max_dim = interface->max_dim;
	int step_width = interface->RJ_step_width;
	int chain_id = interface->chain_number;
	//Sanity check -- only add or remove if you can
	if(mcmc_max_dim == mcmc_min_dim){
	//if(true){
		for(int i = 0 ; i < max_dim; i ++){
			proposed_params[i]=current_params[i];
			proposed_status[i]=current_status[i];
		}
	}
	else if(!mcmc_intrinsic && (mcmc_generation_method == "ppE_IMRPhenomD_Inspiral"||mcmc_generation_method == "ppE_IMRPhenomD_IMR" || mcmc_generation_method=="ppE_IMRPhenomPv2_Inspiral"|| mcmc_generation_method=="ppE_IMRPhenomPv2_IMR")){
		//copy over original
		for(int i = 0 ; i < max_dim; i ++){
			proposed_params[i]=current_params[i];
			proposed_status[i]=current_status[i];
		}
		//Find modifications currently in use
		int mods_in_use = 0 ;
		int mod_ids[mcmc_max_dim - mcmc_min_dim];
		int nonmod_ids[mcmc_max_dim - mcmc_min_dim];
		int mods_not_in_use = 0;
		for(int i = mcmc_min_dim ; i<mcmc_max_dim; i++){
			if(current_status[i] == 1){
				mod_ids[mods_in_use] = i;
				mods_in_use ++;
			}
			else{
				nonmod_ids[mods_not_in_use] = i;
				mods_not_in_use ++;

			}
		}
		//Only add Modification
		if(mods_in_use ==0){
			//pick random modification
			int alpha = (mcmc_max_dim - mcmc_min_dim)*gsl_rng_uniform(mcmc_rvec[chain_id]);
			//Pick random value centered on 0
			double beta = gsl_ran_gaussian(mcmc_rvec[chain_id], step_width);
			proposed_status[mcmc_min_dim+alpha] = 1;
			proposed_params[mcmc_min_dim +alpha] = beta;
			
			
		}
		//Only remove modification
		else if(mods_in_use == (mcmc_max_dim - mcmc_min_dim)){
			//pick random modification
			int alpha = (mcmc_max_dim - mcmc_min_dim)*gsl_rng_uniform(mcmc_rvec[chain_id]);
			//Pick random value centered on 0
			proposed_status[mcmc_min_dim+alpha] = 0;
			proposed_params[mcmc_min_dim +alpha] =0;

		}
		//50/50 to do either
		else{
			double delta = gsl_rng_uniform(mcmc_rvec[chain_id]);
			//add modification
			if(delta<.5){
				//pick random modification that is not in use
				//int alpha = ( (mcmc_max_dim - mcmc_min_dim) - mods_in_use)
				int alpha = ( mods_not_in_use)
					*gsl_rng_uniform(mcmc_rvec[chain_id]);
				//Pick random value centered on 0
				double beta = gsl_ran_gaussian(mcmc_rvec[chain_id], step_width);
				proposed_status[nonmod_ids[alpha]] = 1;
				proposed_params[nonmod_ids[alpha]] = beta;

			}
			//remove modification
			else{
				//pick random modification that is in use
				int alpha = (mods_in_use)*gsl_rng_uniform(mcmc_rvec[chain_id]);
				//Pick random value centered on 0
				proposed_status[mod_ids[alpha]] = 0;
				proposed_params[mod_ids[alpha]] =0;

			}
		}
	}
}

void RJPTMCMC_fisher_wrapper(double *param, int *status,  double **output, mcmc_data_interface *interface,void *parameters)
{
	int min_dim = interface->min_dim;
	int chain_id = interface->chain_number;
	if(!mcmc_intrinsic && (mcmc_generation_method =="ppE_IMRPhenomD_Inspiral"||mcmc_generation_method =="ppE_IMRPhenomD_IMR")){	
		std::string local_method = "IMRPhenomD";
		//unpack parameter vector
		double incl = acos(param[0]);
		double RA = param[1];
		double DEC = param[2];
		double DL = std::exp(param[3]);
		double chirpmass = std::exp(param[4]);
		double eta = param[5];
		double chi1 = param[6];
		double chi2 = param[7];
		double delta_t = 0;
		double tc_ref =0;
		double phic_ref =0;
		double psi = param[8];
	
		double *phi = new double[mcmc_num_detectors];
		double *theta = new double[mcmc_num_detectors];
		//celestial_horizon_transform(RA,DEC, mcmc_gps_time, "Hanford", &phi[0], &theta[0]);

	
		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = 0;
		parameters.spin1[1] = 0;
		parameters.spin1[2] = chi1;
		parameters.spin2[0] = 0;
		parameters.spin2[1] = 0;
		parameters.spin2[2] = chi2;
		parameters.Luminosity_Distance = DL;
		//The rest is maximized over for this option
		parameters.tc = 0;
		parameters.phiRef = 0;
		parameters.incl_angle = incl;
		parameters.phi=0;
		parameters.theta=0;
		parameters.NSflag1 = false;
		parameters.NSflag2 = false;
		parameters.sky_average = false;
		parameters.psi = psi;
		
		for(int j =0; j<min_dim; j++){
			for(int k =0; k<min_dim; k++)
			{
				output[j][k] =0;
			}
		} 
		double **temp_out = allocate_2D_array(min_dim,min_dim);
		for (int i =0; i<mcmc_num_detectors; i++){
			celestial_horizon_transform(RA,DEC, mcmc_gps_time, mcmc_detectors[i], &phi[i], &theta[i]);
			parameters.phi = phi[i];
			parameters.theta = theta[i];
			fisher_numerical(mcmc_frequencies[i], mcmc_data_length[i],
				"MCMC_"+local_method+"_Full", 
				mcmc_detectors[i], mcmc_detectors[0],temp_out, 9, &parameters, mcmc_deriv_order,
				NULL, NULL, mcmc_noise[i]);
			for(int j =0; j<min_dim; j++){
				for(int k =0; k<min_dim; k++)
				{
					output[j][k] +=temp_out[j][k];
				}
			} 
		}
		delete [] phi; delete [] theta;

		deallocate_2D_array(temp_out, min_dim,min_dim);
	}
	else if(!mcmc_intrinsic && (mcmc_generation_method =="ppE_IMRPhenomPv2_Inspiral"||mcmc_generation_method =="ppE_IMRPhenomPv2_IMR")){	
		std::string local_method = "IMRPhenomPv2";
		//unpack parameter vector
		double incl = acos(param[0]);
		double RA = param[1];
		double DEC = param[2];
		double DL = std::exp(param[3]);
		double chirpmass = std::exp(param[4]);
		double eta = param[5];
		double chi1 = param[6];
		double chi2 = param[7];
		double theta1 = param[8];//Polar angles of spins relative to L
		double theta2 = param[9];
		double phi1 = param[10];//Azimuthal angles of spins relative to L
		double phi2 = param[11];
		double phiref = param[12];//Orbital phase at fref (20Hz) -- Someday, find a way to maximize this out
		double psi = param[13];//Orbital phase at fref (20Hz) -- Someday, find a way to maximize this out
		double delta_t = 0;
		double tc_ref =0;
		double fref = 20;
		//transform
		double spin1[3];
		double spin2[3];
		double spin1sphr[3] =  {chi1, theta1, phi1};
		double spin2sphr[3] =  {chi2, theta2, phi2};
		transform_sph_cart( &spin1sphr[0], &spin1[0]);
		transform_sph_cart( &spin2sphr[0], &spin2[0]);
		double *phi = new double[mcmc_num_detectors];
		double *theta = new double[mcmc_num_detectors];
	
		//double *phi = new double[mcmc_num_detectors];
		//double *theta = new double[mcmc_num_detectors];
		//celestial_horizon_transform(RA,DEC, mcmc_gps_time, mcmc_detectors[0], &phi[0], &theta[0]);

		//create gen_param struct
		gen_params parameters; 
		parameters.mass1 = calculate_mass1(chirpmass, eta);
		parameters.mass2 = calculate_mass2(chirpmass, eta);
		parameters.spin1[0] = spin1[0];
		parameters.spin1[1] = spin1[1];
		parameters.spin1[2] = spin1[2];
		parameters.spin2[0] = spin2[0];
		parameters.spin2[1] = spin2[1];
		parameters.spin2[2] = spin2[2];
		parameters.Luminosity_Distance = DL;
		//The rest is maximized over for this option
		parameters.tc = 0;
		//parameters.phiRef = 0;
		parameters.phiRef=phiref;
		parameters.f_ref=fref;
		parameters.incl_angle = incl;
		parameters.psi = psi;
		//parameters.phi=phi[0];
		//parameters.theta=theta[0];
		//parameters.phi=0;
		//parameters.theta=0;
		parameters.RA = RA;
		parameters.DEC = DEC;
		parameters.gmst = mcmc_gmst;
		parameters.NSflag1 = false;
		parameters.NSflag2 = false;
		parameters.sky_average = false;
		
		for(int j =0; j<mcmc_min_dim; j++){
			for(int k =0; k<mcmc_min_dim; k++)
			{
				output[j][k] =0;
			}
		} 
		double **temp_out = allocate_2D_array(mcmc_min_dim,mcmc_min_dim);
		for (int i =0; i<mcmc_num_detectors; i++){
			celestial_horizon_transform(RA,DEC, mcmc_gps_time, mcmc_detectors[i], &phi[i], &theta[i]);
			parameters.phi = phi[i];
			parameters.theta = theta[i];
			fisher_numerical(mcmc_frequencies[i], mcmc_data_length[i],
				"MCMC_"+local_method+"_Full", 
				mcmc_detectors[i], mcmc_detectors[0],temp_out, 14, &parameters, mcmc_deriv_order,
				NULL, NULL, mcmc_noise[i]);
			//double dphi_dra, dtheta_dra,dphi_ddec, dtheta_ddec;
			//derivative_celestial_horizon_transform(RA,DEC,mcmc_gps_time,
			//	mcmc_detectors[i], &dphi_dra, &dtheta_dra, 
			//	&dphi_ddec,&dtheta_ddec);
			for(int j =0; j<mcmc_min_dim; j++){
				for(int k =0; k<mcmc_min_dim; k++)
				{
					output[j][k] +=temp_out[j][k];
					//std::cout<<j<<" "<<k<<" "<<output[j][k]<<std::endl;
				}
			} 
		}

		delete [] phi; delete [] theta;

		deallocate_2D_array(temp_out, mcmc_min_dim,mcmc_min_dim);
	}

}
			

