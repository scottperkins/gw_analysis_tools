#include <fisher.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <math.h>
#include <string>
#include "util.h"
#include "detector_util.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "waveform_generator.h"
#include "waveform_util.h"


using namespace std;





/*!\file 
 *
 * All subroutines associated with waveform differentiation and Fisher analysis
 */

/*!\brief Calculates the fisher matrix for the given arguments
 *
 * Utilizes numerical derivatives
 */
void fisher(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	string generation_method, 
	string detector, 
	double **output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params *parameters,
	//double *parameters,
	int *amp_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method -- if using numerical derivatives or speed isn't that important, just set to NULL*/
	int *phase_tapes,/**< if speed is required, precomputed tapes can be used - assumed the user knows what they're doing, no checks done here to make sure that the number of tapes matches the requirement by the generation_method*/
	double *noise
	)
{
	//populate noise and frequency
	double internal_noise[length];
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
	double **amplitude_deriv = (double **)malloc(dimension*sizeof(**amplitude_deriv));
	for (int i = 0; i<dimension; i++)
		amplitude_deriv[i] = (double *)malloc(length*sizeof(double));
	double **phase_deriv = (double **)malloc(dimension*sizeof(**phase_deriv));
	for (int i = 0; i<dimension; i++)
		phase_deriv[i] = (double *)malloc(length*sizeof(double));
	double *amplitude = (double*)malloc(length*sizeof(double));
	double *integrand = (double*)malloc(length*sizeof(double));
	

	calculate_derivatives(amplitude_deriv, 
			phase_deriv, 
			amplitude,
			frequency,
			length, 
			detector, 
			generation_method,
			parameters);

	//calulate fisher elements
	for (int j=0;j<dimension; j++)
	{
		for (int k = 0; k<j; k++)
		{
			for (int i =0;i<length;i++)
			{
				integrand[i] = 
					real( (amplitude_deriv[j][i]*amplitude_deriv[k][i]
					+amplitude[i]*amplitude[i]*
					phase_deriv[j][i]*phase_deriv[k][i])/internal_noise[i]);
			}
			output[j][k] = 4*simpsons_sum(
						frequency[1]-frequency[0], length, integrand);	
			output[k][j] = output[j][k];
		}

	}

	for (int j = 0; j<dimension; j ++)
	{

		for (int i =0;i<length;i++)
			integrand[i] = 
				real( (amplitude_deriv[j][i]*amplitude_deriv[j][i]
					+amplitude[i]*amplitude[i]*phase_deriv[j][i]*
					phase_deriv[j][i])/internal_noise[i]);
		output[j][j] = 4*simpsons_sum(
					frequency[1]-frequency[0], length, integrand);	
	}
		


	for (int i =0;i<dimension;i++)
	{
		free( amplitude_deriv[i]);
		free( phase_deriv[i]);
	}
	free(amplitude_deriv);
	free(phase_deriv);
	free(amplitude);
	free(integrand);
}


/*! \brief Abstraction layer for handling the case separation for the different waveforms
 *
 */
void calculate_derivatives(double  **amplitude_deriv, 
       	double **phase_deriv,
       	double *amplitude,
       	double *frequencies,
       	int length, 
       	string detector, 
       	string  gen_method,
       	gen_params *parameters)
{
	//Finite difference spacing
	double epsilon = 1e-7;
	double epsilonnaught = 1e-7;
	double *amplitude_plus_plus = (double *) malloc(sizeof(double)*length);
	double *amplitude_plus_minus = (double *) malloc(sizeof(double)*length);
	double *amplitude_cross_plus = (double *) malloc(sizeof(double)*length);
	double *amplitude_cross_minus = (double *) malloc(sizeof(double)*length);
	double *phase_plus_plus = (double *) malloc(sizeof(double)*length);
	double *phase_plus_minus = (double *) malloc(sizeof(double)*length);
	double *phase_cross_plus = (double *) malloc(sizeof(double)*length);
	double *phase_cross_minus = (double *) malloc(sizeof(double)*length);
	
	if (parameters->sky_average && gen_method == "IMRPhenomD"){
		IMRPhenomD<double> model;
		int dimension = 7;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		param_in[0] = parameters_in.A0;//seconds
		param_in[1] = parameters_in.tc;
		param_in[2] = parameters_in.phic;
		param_in[3] = parameters_in.chirpmass;//seconds
		param_in[4] = parameters_in.eta;
		param_in[5] = parameters_in.chi_s;
		param_in[6] = parameters_in.chi_a;

		waveform_params.sky_average=parameters->sky_average;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			model.change_parameter_basis(param_p, param_out, parameters_in.sky_average);
			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;


			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				//amplitude_cross_plus,
				gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				//amplitude_cross_plus,
				gen_method,
				&waveform_params);	

			model.change_parameter_basis(param_m, param_out,parameters_in.sky_average);
			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;
			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				//amplitude_cross_plus,
				gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				//amplitude_cross_plus,
				gen_method,
				&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[0][l] = amplitude_deriv[0][l]*param_in[0] ;
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
	}
	else if (!parameters->sky_average  && gen_method == "IMRPhenomD"){
		IMRPhenomD<double> model;
		int dimension = 11;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		param_in[0] = parameters->incl_angle;
		param_in[1] = parameters->RA;
		param_in[2] = parameters->DEC;
		param_in[3] = parameters_in.DL/MPC_SEC;//MPC
		param_in[4] = parameters_in.chirpmass;//seconds
		param_in[5] = parameters_in.eta;
		param_in[6] = parameters_in.spin1z;
		param_in[7] = parameters_in.spin2z;
		param_in[8] = parameters->phiRef;
		param_in[9] = parameters_in.tc;
		param_in[10] = parameters->psi;

		waveform_params.sky_average=parameters->sky_average;
		waveform_params.f_ref=parameters->f_ref;
		waveform_params.NSflag = parameters->NSflag;
		waveform_params.gmst = parameters->gmst;
		waveform_params.shift_time =false;
		std::complex<double> *response = new std::complex<double>[length];
		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			waveform_params.mass1 = calculate_mass1(param_p[4],param_p[5])/MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_p[4],param_p[5])/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_p[3];
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_p[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_p[7];
			waveform_params.phiRef = param_p[8];
			waveform_params.tc= param_p[9];
			waveform_params.incl_angle = param_p[0];
			waveform_params.RA = param_p[1];
			waveform_params.DEC = param_p[2];
			waveform_params.psi = param_p[10];

			fourier_detector_response_equatorial(frequencies, length, 
				response,detector, gen_method, &waveform_params);

			for (int k =0; k<length; k++){
				amplitude_plus_plus[k] =  std::abs(response[k]);
				phase_plus_plus[k] =  std::arg(response[k]);
			}
			//##############################################
			waveform_params.mass1 = calculate_mass1(param_m[4],param_m[5])/MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_m[4],param_m[5])/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_m[3];
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_m[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_m[7];
			waveform_params.phiRef = param_m[8];
			waveform_params.tc= param_m[9];
			waveform_params.incl_angle = param_m[0];
			waveform_params.RA = param_m[1];
			waveform_params.DEC = param_m[2];
			waveform_params.psi = param_m[10];

			fourier_detector_response_equatorial(frequencies, length, 
				response,detector, gen_method, &waveform_params);

			for (int k =0; k<length; k++){
				amplitude_plus_minus[k] =  std::abs(response[k]);
				phase_plus_minus[k] =  std::arg(response[k]);
			}
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
		delete [] response;
	}
	if (gen_method == "ppE_IMRPhenomD_Inspiral"|| gen_method=="ppE_IMRPhenomD_IMR"){
		IMRPhenomD<double> model;
		int dimension = 7+parameters->Nmod;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		param_in[0] = parameters_in.A0;//seconds
		param_in[1] = parameters_in.tc;
		param_in[2] = parameters_in.phic;
		param_in[3] = parameters_in.chirpmass;//seconds
		param_in[4] = parameters_in.eta;
		param_in[5] = parameters_in.chi_s;
		param_in[6] = parameters_in.chi_a;
		for(int i = 0 ; i<parameters->Nmod; i++){
			param_in[7+i] = parameters->betappe[i];
		}

		waveform_params.sky_average=parameters->sky_average;
		waveform_params.bppe = parameters->bppe;
		waveform_params.Nmod = parameters->Nmod;
		waveform_params.NSflag = parameters->NSflag;
		waveform_params.betappe = new double[waveform_params.Nmod];

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			model.change_parameter_basis(param_p, param_out, parameters_in.sky_average);
			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			for(int j = 0 ; j<waveform_params.Nmod ; j++){
				waveform_params.betappe[j] = param_p[7+j];
			}


			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				gen_method,
				&waveform_params);	

			model.change_parameter_basis(param_m, param_out,parameters_in.sky_average);
			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			for(int j = 0 ; j<waveform_params.Nmod ; j++){
				waveform_params.betappe[j] = param_m[7+j];
			}
			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				gen_method,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				gen_method,
				&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
				
		}
		delete [] waveform_params.betappe;
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[0][l] = amplitude_deriv[0][l]*param_in[0] ;
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
	}
	//*NOTE* this is not a good, rigorous fisher matrix, as some extrinsic parameters
	//are chosen randomly. This is only to inform MCMC steps, and is just an estimate
	if (gen_method == "MCMC_IMRPhenomPv2"){
		std::string local_method = "IMRPhenomPv2";
		double *temp = (double*)malloc(sizeof(double)*length);
		fourier_detector_amplitude_phase(frequencies, 
			length,
			amplitude,
			temp,
			//amplitude_cross_plus,
			detector,
			local_method,
			parameters);	
		IMRPhenomPv2<double> model;
		int dimension = 9;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		double spin1sph[3];
		double spin2sph[3];
		transform_cart_sph(parameters->spin1, spin1sph);
		transform_cart_sph(parameters->spin2, spin2sph);
		double chirpmass = calculate_chirpmass(parameters->mass1, parameters->mass2);
		double eta = calculate_eta(parameters->mass1, parameters->mass2);
		param_in[0] = cos(parameters->incl_angle);
		param_in[1] = chirpmass; //Sol mass
		param_in[2] = eta;
		param_in[3] = spin1sph[0];
		param_in[4] = spin2sph[0];
		param_in[5] = spin1sph[1];
		param_in[6] = spin2sph[1];
		param_in[7] = spin1sph[2];
		param_in[8] = spin2sph[2];

		waveform_params.sky_average=parameters->sky_average;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;
			if(param_p[2]>.25) param_p[2] = .25;

				
			waveform_params.mass1 =calculate_mass1(param_p[1],param_p[2]);
			waveform_params.mass2 =calculate_mass2(param_p[1],param_p[2]);
			waveform_params.Luminosity_Distance=parameters->Luminosity_Distance;
			double param_in_spin1[3] = {param_p[3],param_p[5],param_p[7]};
			transform_sph_cart(param_in_spin1, waveform_params.spin1);
			double param_in_spin2[3] = {param_p[4],param_p[6],param_p[8]};
			transform_sph_cart(param_in_spin2, waveform_params.spin2);
			waveform_params.phic = parameters->phic;
			waveform_params.tc=parameters->tc;
			waveform_params.incl_angle = std::acos(param_p[0]);
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;


			fourier_detector_amplitude_phase(frequencies, 
				length,
				amplitude_plus_plus,
				phase_plus_plus,
				//amplitude_cross_plus,
				detector,
				gen_method,
				&waveform_params);	
			//fourier_phase(frequencies, 
			//	length,
			//	phase_plus_plus,
			//	//amplitude_cross_plus,
			//	gen_method,
			//	&waveform_params);	

			waveform_params.mass1 =calculate_mass1(param_m[1],param_m[2]);
			waveform_params.mass2 =calculate_mass2(param_m[1],param_m[2]);
			waveform_params.Luminosity_Distance=parameters->Luminosity_Distance;

			param_in_spin1[0] = param_m[3];
			param_in_spin1[1] = param_m[5];
			param_in_spin1[2]=param_m[7];
			transform_sph_cart(param_in_spin1, waveform_params.spin1);
			param_in_spin1[0] = param_m[4];
			param_in_spin1[1] = param_m[6];
			param_in_spin1[2]=param_m[8];
			transform_sph_cart(param_in_spin2, waveform_params.spin2);

			waveform_params.phic = parameters->phic;
			waveform_params.tc=parameters->tc;
			waveform_params.incl_angle = std::acos(param_m[0]);
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;
			fourier_detector_amplitude_phase(frequencies, 
				length,
				amplitude_plus_minus,
				phase_plus_minus,
				//amplitude_cross_plus,
				detector,
				gen_method,
				&waveform_params);	
			//fourier_detector_phase(frequencies, 
			//	length,
			//	phase_plus_minus,
			//	//amplitude_cross_plus,
			//	gen_method,
			//	&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
			//ofstream ampfile
			//for (int l = 0;l<length;l++)
			//{
			//	
			//}
			
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[0][l] = amplitude_deriv[0][l]*param_in[0] ;
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
		free(temp);
	}
	if (gen_method == "MCMC_IMRPhenomD_Full"){
		fourier_amplitude(frequencies, 
			length,
			amplitude,
			//amplitude_cross_plus,
			"IMRPhenomD",
			parameters);	
		std::complex<double> Qtemp = 
			Q(parameters->theta, parameters->phi, parameters->incl_angle, parameters->psi);
		double a_corr = std::abs(Qtemp);
		for (int k =0; k<length; k++)
			amplitude[k] = a_corr * amplitude[k];

		IMRPhenomD<double> model;
		int dimension = 9;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		double chirpmass = calculate_chirpmass(parameters->mass1, parameters->mass2);
		double eta = calculate_eta(parameters->mass1, parameters->mass2);
		param_in[0] = cos(parameters->incl_angle);
		param_in[1] = parameters->theta;
		param_in[2] = parameters->phi;
		param_in[3] = parameters->Luminosity_Distance;//MPC
		param_in[4] = chirpmass; //Sol mass
		param_in[5] = eta;
		param_in[6] = parameters->spin1[2];
		param_in[7] = parameters->spin2[2];
		param_in[8] = parameters->psi;
		waveform_params.sky_average=parameters->sky_average;
		double m1, m2,Fpp,Fcc, a_corr_p, a_corr_m, p_corr_p, p_corr_m;
		std::complex<double> Qp, Qm;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			//cos \iota must lie within certain range
			if(param_p[0]>1.)param_p[0]=1.;
			else if( param_p[1]<-1.) param_p[0]=-1.;
			if(param_m[0]>1.)param_m[0]=1.;
			else if( param_m[1]<-1.) param_m[0]=-1.;

			
			waveform_params.mass1 = calculate_mass1(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_p[3];//MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_p[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_p[7];
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_p[0]);
			waveform_params.theta = param_p[1];
			waveform_params.phi = param_p[2];
			waveform_params.NSflag = parameters->NSflag;
			waveform_params.psi = param_p[8];


			Qp = Q(waveform_params.theta, waveform_params.phi, waveform_params.incl_angle, waveform_params.psi);
			a_corr_p = std::abs(Qp);
			p_corr_p = std::arg(Qp);
			

			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	


			waveform_params.mass1 = calculate_mass1(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_m[3];//MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_m[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_m[7];
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_m[0]);
			waveform_params.theta = param_m[1];
			waveform_params.phi = param_m[2];
			waveform_params.NSflag = parameters->NSflag;
			waveform_params.psi = param_m[8];

			Qm = Q(waveform_params.theta, waveform_params.phi, waveform_params.incl_angle, waveform_params.psi);
			a_corr_m = std::abs(Qm);
			p_corr_m = std::arg(Qm);

			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (a_corr_p*amplitude_plus_plus[l] -a_corr_m*amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (p_corr_p+phase_plus_plus[l] -p_corr_m -phase_plus_minus[l])/(2*epsilon);
			}
			
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
	}
	if (gen_method == "MCMC_IMRPhenomPv2_Full"){
		std::string local_method = "IMRPhenomPv2";
		std::complex<double> *response = new std::complex<double> [length];
		fourier_detector_response(frequencies, 
			length,
			response,
			detector,
			local_method,
			parameters);	
		for (int k =0; k<length; k++){
			amplitude[k] =  std::abs(response[k]);
		}

		IMRPhenomD<double> model;
		int dimension = 14;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		double chirpmass = calculate_chirpmass(parameters->mass1, parameters->mass2);
		double eta = calculate_eta(parameters->mass1, parameters->mass2);
		double spin1sph[3];
		double spin2sph[3];
		
		transform_cart_sph(parameters->spin1, spin1sph);
		transform_cart_sph(parameters->spin2, spin2sph);
		
		param_in[0] = cos(parameters->incl_angle);
		param_in[1] = parameters->RA;
		param_in[2] = parameters->DEC;
		param_in[3] = parameters->Luminosity_Distance;//MPC
		param_in[4] = chirpmass; //Sol mass
		param_in[5] = eta;
		param_in[6] = spin1sph[0];
		param_in[7] = spin2sph[0];
		param_in[8] = spin1sph[1];
		param_in[9] = spin2sph[1];
		param_in[10] = spin1sph[2];
		param_in[11] = spin2sph[2];
		param_in[12] = parameters->phiRef;
		param_in[13] = parameters->psi;
		waveform_params.sky_average=parameters->sky_average;
		waveform_params.NSflag=parameters->NSflag;
		waveform_params.gmst=parameters->gmst;
		waveform_params.f_ref = parameters->f_ref;
		double m1, m2,Fpp,Fcc, a_corr_p, a_corr_m, p_corr_p, p_corr_m;
		std::complex<double> Qp, Qm;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			//cos \iota must lie within certain range
			if(param_p[0]>1.)param_p[0]=1.;
			else if( param_p[1]<-1.) param_p[0]=-1.;
			if(param_m[0]>1.)param_m[0]=1.;
			else if( param_m[1]<-1.) param_m[0]=-1.;
			if(param_p[5]>.25) param_p[5] = .25;
			if(param_p[6]>.95) param_p[6] = .95;
			if(param_p[6]<-.95) param_p[6] = -.95;
			if(param_p[7]>.95) param_p[7] = .95;
			if(param_p[7]<-.95) param_p[7] = -.95;

			waveform_params.mass1 = calculate_mass1(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_p[3];//MPC_SEC;
			double param_in_spin1p[3] = {param_p[6],param_p[8],param_p[10]};
			transform_sph_cart(param_in_spin1p, waveform_params.spin1);
			double param_in_spin2p[3] = {param_p[7],param_p[9],param_p[11]};
			transform_sph_cart(param_in_spin2p, waveform_params.spin2);
			//waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_p[0]);
			//waveform_params.theta = param_p[1];
			//waveform_params.phi = param_p[2];
			waveform_params.RA = param_p[1];
			waveform_params.DEC = param_p[2];
			
			waveform_params.phiRef = param_p[12];
			waveform_params.psi = param_p[13];

			fourier_detector_response_equatorial(frequencies, 
				length,
				response,
				detector,
				local_method,
				&waveform_params);	
			for (int k =0; k<length; k++){
				amplitude_plus_plus[k] =  std::abs(response[k]);
				phase_plus_plus[k] =  std::arg(response[k]);
			}


			waveform_params.mass1 = calculate_mass1(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_m[3];//MPC_SEC;
			double param_in_spin1m[3] = {param_m[6],param_m[8],param_m[10]};
			transform_sph_cart(param_in_spin1m, waveform_params.spin1);
			double param_in_spin2m[3] = {param_m[7],param_m[9],param_m[11]};
			transform_sph_cart(param_in_spin2m, waveform_params.spin2);
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_m[0]);
			//waveform_params.theta = param_m[1];
			//waveform_params.phi = param_m[2];
			waveform_params.RA = param_m[1];
			waveform_params.DEC = param_m[2];
			waveform_params.phiRef = param_m[12];
			waveform_params.psi = param_m[13];

			fourier_detector_response_equatorial(frequencies, 
				length,
				response,
				detector,
				local_method,
				&waveform_params);	
			for (int k =0; k<length; k++){
				amplitude_plus_minus[k] =  std::abs(response[k]);
				phase_plus_minus[k] =  std::arg(response[k]);
			}

			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
			
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
		delete [] response;
	}
	if (gen_method == "MCMC_ppE_IMRPhenomPv2_Inspiral_Full"||
		gen_method =="MCMC_ppE_IMRPhenomPv2_IMR_Full" ||
		gen_method =="MCMC_dCS_IMRPhenomPv2_root_alpha_Full" ||
		gen_method =="MCMC_EdGB_IMRPhenomPv2_root_alpha_Full" 
		){
		std::string local_gen;
		bool log_scaling = false;
		if(gen_method == "MCMC_dCS_IMRPhenomPv2_log_Full"){
			local_gen = "dCS_IMRPhenomPv2";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_dCS_IMRPhenomPv2_Full"){
			local_gen = "dCS_IMRPhenomPv2";
		}
		else if(gen_method == "MCMC_dCS_IMRPhenomPv2_root_alpha_Full"){
			local_gen = "dCS_IMRPhenomPv2";
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomPv2_log_Full"){
			local_gen = "EdGB_IMRPhenomPv2";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomPv2_Full"){
			local_gen = "EdGB_IMRPhenomPv2";
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomPv2_root_alpha_Full"){
			local_gen = "EdGB_IMRPhenomPv2";
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomPv2_IMR_log_Full"){
			local_gen = "ppE_IMRPhenomPv2_IMR";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomPv2_Inspiral_log_Full"){
			local_gen = "ppE_IMRPhenomPv2_Inspiral";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomPv2_IMR_Full"){
			local_gen = "ppE_IMRPhenomPv2_IMR";
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomPv2_Inspiral_Full"){
			local_gen = "ppE_IMRPhenomPv2_Inspiral";
		}
		std::complex<double> *response = new std::complex<double> [length];
		fourier_detector_response(frequencies, 
			length,
			response,
			detector,
			local_gen,
			parameters);	
		for (int k =0; k<length; k++){
			amplitude[k] =  std::abs(response[k]);
		}

		int dimension = 14+parameters->Nmod;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		double chirpmass = calculate_chirpmass(parameters->mass1, parameters->mass2);
		double eta = calculate_eta(parameters->mass1, parameters->mass2);
		double spin1sph[3];
		double spin2sph[3];
		
		transform_cart_sph(parameters->spin1, &spin1sph[0]);
		transform_cart_sph(parameters->spin2, &spin2sph[0]);
		
		param_in[0] = cos(parameters->incl_angle);
		param_in[1] = parameters->RA;
		param_in[2] = parameters->DEC;
		param_in[3] = parameters->Luminosity_Distance;//MPC
		param_in[4] = chirpmass; //Sol mass
		param_in[5] = eta;
		param_in[6] = spin1sph[0];
		param_in[7] = spin2sph[0];
		param_in[8] = spin1sph[1];
		param_in[9] = spin2sph[1];
		param_in[10] = spin1sph[2];
		param_in[11] = spin2sph[2];
		param_in[12] = parameters->phiRef;
		param_in[13] = parameters->psi;
		for(int i = 0; i<parameters->Nmod; i++){
			param_in[14+i] = parameters->betappe[i];
		}
		waveform_params.sky_average=parameters->sky_average;
		waveform_params.bppe = parameters->bppe;//new int[parameters->Nmod];
		//for(int i =0 ;i<parameters->Nmod; i++){
		//	waveform_params.bppe[i] = parameters->bppe[i];
		//}
		double m1, m2,Fpp,Fcc, a_corr_p, a_corr_m, p_corr_p, p_corr_m;
		std::complex<double> Qp, Qm;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			//cos \iota must lie within certain range
			if(param_p[0]>1.)param_p[0]=1.;
			else if( param_p[1]<-1.) param_p[0]=-1.;
			if(param_m[0]>1.)param_m[0]=1.;
			else if( param_m[1]<-1.) param_m[0]=-1.;
			if(param_p[5]>.25) param_p[5] = .25;
			if(param_p[6]>.95) param_p[6] = .95;
			if(param_p[6]<-.95) param_p[6] = -.95;
			if(param_p[7]>.95) param_p[7] = .95;
			if(param_p[7]<-.95) param_p[7] = -.95;

			
			waveform_params.mass1 = calculate_mass1(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_p[3];//MPC_SEC;
			double param_in_spin1p[3] = {param_p[6],param_p[8],param_p[10]};
			transform_sph_cart(param_in_spin1p, waveform_params.spin1);
			double param_in_spin2p[3] = {param_p[7],param_p[9],param_p[11]};
			transform_sph_cart(param_in_spin2p, waveform_params.spin2);
			//waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_p[0]);
			//waveform_params.theta = param_p[1];
			//waveform_params.phi = param_p[2];
			waveform_params.RA = param_p[1];
			waveform_params.DEC = param_p[2];
			waveform_params.gmst = parameters->gmst;
			waveform_params.phiRef = param_p[12];
			waveform_params.f_ref = parameters->f_ref;
			waveform_params.NSflag = parameters->NSflag;
			waveform_params.sky_average = parameters->sky_average;
			waveform_params.psi = param_p[13];
			waveform_params.betappe = new double[parameters->Nmod];
			if(gen_method == "MCMC_dCS_IMRPhenomPv2_root_alpha_Full"||gen_method == "MCMC_EdGB_IMRPhenomPv2_root_alpha_Full"){
				waveform_params.betappe[0]=pow(param_p[14]/(3.e5), 4.);
			}


			fourier_detector_response_equatorial(frequencies, 
				length,
				response,
				detector,
				local_gen,
				&waveform_params);	
			for (int k =0; k<length; k++){
				amplitude_plus_plus[k] =  std::abs(response[k]);
				phase_plus_plus[k] =  std::arg(response[k]);
			}


			waveform_params.mass1 = calculate_mass1(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_m[3];//MPC_SEC;
			double param_in_spin1m[3] = {param_m[6],param_m[8],param_m[10]};
			transform_sph_cart(param_in_spin1m, waveform_params.spin1);
			double param_in_spin2m[3] = {param_m[7],param_m[9],param_m[11]};
			transform_sph_cart(param_in_spin2m, waveform_params.spin2);
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_m[0]);
			//waveform_params.theta = param_m[1];
			//waveform_params.phi = param_m[2];
			waveform_params.RA = param_m[1];
			waveform_params.DEC = param_m[2];
			waveform_params.gmst = parameters->gmst;
			waveform_params.phiRef = param_m[12];
			waveform_params.f_ref = parameters->f_ref;
			waveform_params.NSflag = parameters->NSflag;
			waveform_params.psi = param_m[13];
			for (int j =0; j<parameters->Nmod; j++){
				waveform_params.betappe[j] = param_m[14+j];
			}
			if(gen_method == "MCMC_dCS_IMRPhenomPv2_root_alpha_Full"||gen_method == "MCMC_EdGB_IMRPhenomPv2_root_alpha_Full"){
				waveform_params.betappe[0]=pow(param_m[14]/(3.e5), 4.);
			}

			fourier_detector_response_equatorial(frequencies, 
				length,
				response,
				detector,
				local_gen,
				&waveform_params);	
			for (int k =0; k<length; k++){
				amplitude_plus_minus[k] =  std::abs(response[k]);
				phase_plus_minus[k] =  std::arg(response[k]);
			}

			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
			
			delete [] waveform_params.betappe;
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
		delete [] response;
		//delete [] waveform_params.bppe;
	}
	else if (gen_method == "MCMC_dCS_IMRPhenomD_log_Full" 
		|| gen_method == "MCMC_dCS_IMRPhenomD_Full"
		|| gen_method == "MCMC_dCS_IMRPhenomD_root_alpha_Full"
		|| gen_method == "MCMC_EdGB_IMRPhenomD_log_Full"
		|| gen_method == "MCMC_EdGB_IMRPhenomD_Full"
		|| gen_method == "MCMC_EdGB_IMRPhenomD_root_alpha_Full"
		|| gen_method == "MCMC_ppE_IMRPhenomD_Inspiral_log_Full"
		|| gen_method == "MCMC_ppE_IMRPhenomD_IMR_log_Full"
		|| gen_method == "MCMC_ppE_IMRPhenomD_Inspiral_Full"
		|| gen_method == "MCMC_ppE_IMRPhenomD_IMR_Full"
		){
		std::string local_gen;
		bool log_scaling = false;
		if(gen_method == "MCMC_dCS_IMRPhenomD_log_Full"){
			local_gen = "dCS_IMRPhenomD";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_dCS_IMRPhenomD_Full"){
			local_gen = "dCS_IMRPhenomD";
		}
		else if(gen_method == "MCMC_dCS_IMRPhenomD_root_alpha_Full"){
			local_gen = "dCS_IMRPhenomD";
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomD_log_Full"){
			local_gen = "EdGB_IMRPhenomD";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomD_Full"){
			local_gen = "EdGB_IMRPhenomD";
		}
		else if(gen_method == "MCMC_EdGB_IMRPhenomD_root_alpha_Full"){
			local_gen = "EdGB_IMRPhenomD";
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomD_IMR_log_Full"){
			local_gen = "ppE_IMRPhenomD_IMR";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomD_Inspiral_log_Full"){
			local_gen = "ppE_IMRPhenomD_Inspiral";
			log_scaling = true;
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomD_IMR_Full"){
			local_gen = "ppE_IMRPhenomD_IMR";
		}
		else if(gen_method == "MCMC_ppE_IMRPhenomD_Inspiral_Full"){
			local_gen = "ppE_IMRPhenomD_Inspiral";
		}
		fourier_amplitude(frequencies, 
			length,
			amplitude,
			//amplitude_cross_plus,
			local_gen,
			parameters);	
		std::complex<double> Qtemp = 
			Q(parameters->theta, parameters->phi, parameters->incl_angle);
		double a_corr = std::abs(Qtemp);
		for (int k =0; k<length; k++)
			amplitude[k] = a_corr * amplitude[k];

		IMRPhenomD<double> model;
		int dimension = 8+parameters->Nmod;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		double chirpmass = calculate_chirpmass(parameters->mass1, parameters->mass2);
		double eta = calculate_eta(parameters->mass1, parameters->mass2);
		param_in[0] = cos(parameters->incl_angle);
		param_in[1] = parameters->theta;
		param_in[2] = parameters->phi;
		param_in[3] = parameters->Luminosity_Distance;//MPC
		param_in[4] = chirpmass; //Sol mass
		param_in[5] = eta;
		param_in[6] = parameters->spin1[2];
		param_in[7] = parameters->spin2[2];
		for(int i = 0; i<parameters->Nmod; i++){
			param_in[8+i] = parameters->betappe[i];
		}
		waveform_params.sky_average=parameters->sky_average;
		waveform_params.bppe = new int[parameters->Nmod];
		for(int i =0 ;i<parameters->Nmod; i++){
			waveform_params.bppe[i] = parameters->bppe[i];
		}
		double m1, m2,Fpp,Fcc, a_corr_p, a_corr_m, p_corr_p, p_corr_m;
		std::complex<double> Qp, Qm;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;

			//cos \iota must lie within certain range
			if(param_p[0]>1.)param_p[0]=1.;
			else if( param_p[1]<-1.) param_p[0]=-1.;
			if(param_m[0]>1.)param_m[0]=1.;
			else if( param_m[1]<-1.) param_m[0]=-1.;
			if(param_p[5]>.25) param_p[5] = .25;

			
			waveform_params.mass1 = calculate_mass1(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_p[4],param_p[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_p[3];//MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_p[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_p[7];
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_p[0]);
			waveform_params.theta = param_p[1];
			waveform_params.phi = param_p[2];
			waveform_params.NSflag = parameters->NSflag;
			waveform_params.betappe = new double[parameters->Nmod];
			for (int j =0; j<parameters->Nmod; j++){
				waveform_params.betappe[j] = param_p[8+j];
			}
			if(gen_method == "MCMC_dCS_IMRPhenomD_root_alpha_Full"||gen_method == "MCMC_EdGB_IMRPhenomD_root_alpha_Full"){
				waveform_params.betappe[0]=pow(param_p[8]/(3.e5), 4.);
			}


			Qp = Q(waveform_params.theta, waveform_params.phi, waveform_params.incl_angle);
			a_corr_p = std::abs(Qp);
			p_corr_p = std::arg(Qp);
			

			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				//amplitude_cross_plus,
				local_gen,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				//amplitude_cross_plus,
				local_gen,
				&waveform_params);	


			waveform_params.mass1 = calculate_mass1(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.mass2 = calculate_mass2(param_m[4],param_m[5]);//MSOL_SEC;
			waveform_params.Luminosity_Distance=param_m[3];//MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_m[6];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_m[7];
			waveform_params.phic = 0;
			waveform_params.tc= 0;
			waveform_params.incl_angle = acos(param_m[0]);
			waveform_params.theta = param_m[1];
			waveform_params.phi = param_m[2];
			waveform_params.NSflag = parameters->NSflag;
			for (int j =0; j<parameters->Nmod; j++){
				waveform_params.betappe[j] = param_m[8+j];
			}
			if(gen_method == "MCMC_dCS_IMRPhenomD_root_alpha_Full"||gen_method == "MCMC_EdGB_IMRPhenomD_root_alpha_Full"){
				waveform_params.betappe[0]=pow(param_m[8]/(3.e5), 4.);
			}

			Qm = Q(waveform_params.theta, waveform_params.phi, waveform_params.incl_angle);
			a_corr_m = std::abs(Qm);
			p_corr_m = std::arg(Qm);

			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				local_gen,
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				local_gen,
				&waveform_params);	
			delete [] waveform_params.betappe;
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (a_corr_p*amplitude_plus_plus[l] -a_corr_m*amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (p_corr_p+phase_plus_plus[l] -p_corr_m -phase_plus_minus[l])/(2*epsilon);
			}
			
		
				
		}
		delete [] waveform_params.bppe;
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
			if(log_scaling){
				for (int j =0 ; j<parameters->Nmod; j++){
					amplitude_deriv[8+j][l] = 
						amplitude_deriv[8+j][l]*param_in[8+j];	
					phase_deriv[8+j][l] = 
						phase_deriv[8+j][l]*param_in[8+j];	
				}
			}
		}
	}
	else if (gen_method == "MCMC_IMRPhenomD"){
		IMRPhenomD<double> model;
		int dimension = 7;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		//parameters_in = parameters_in.populate_source_parameters(parameters->mass1, 
		//		parameters->mass2,
		//		parameters->Luminosity_Distance,
		//		parameters->spin1,
		//		parameters->spin2,
		//		parameters->phic,
		//		parameters->tc);
		
		//#################################
		//TESTING 
		//parameters->sky_average=true;
		//#################################
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		param_in[0] = parameters_in.DL/MPC_SEC;//MPC -- too big a number for numdiff?
		param_in[1] = parameters_in.tc;
		param_in[2] = parameters_in.phic;
		param_in[3] = parameters_in.chirpmass;//seconds
		param_in[4] = parameters_in.eta;
		param_in[5] = parameters_in.chi_s+parameters_in.chi_a;
		param_in[6] = parameters_in.chi_s-parameters_in.chi_a;

		waveform_params.sky_average=parameters->sky_average;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			//if(i==0) epsilon = 1e-25;
			//else epsilon = epsilonnaught;
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;
			//param_p[i] = param_in[i]*(1. + epsilon);
			//param_m[i] = param_in[i] *(1- epsilon);

			if(param_p[4]>.25) param_p[4] = .25;

			model.change_parameter_basis(param_p, param_out, parameters_in.sky_average);
			//Spins already individual, reverse, but change_param_basis
			//assumes you input chi_s and chi_a
			param_out[2] = param_p[0];
			param_out[3] = param_p[5];	
			param_out[4] = param_p[6];	

			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2];
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;


			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	

			model.change_parameter_basis(param_m, param_out,parameters_in.sky_average);
			//Spins already individual, reverse, but change_param_basis
			//assumes you input chi_s and chi_a
			param_out[3] = param_m[5];	
			param_out[4] = param_m[6];	
			param_out[2] = param_p[0];

			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[2];
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = param_out[5];
			waveform_params.tc= param_out[6];
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;
			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
			//ofstream ampfile
			//for (int l = 0;l<length;l++)
			//{
			//	
			//}
			
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[0][l] = amplitude_deriv[0][l]*param_in[0] ;
			amplitude_deriv[3][l] = amplitude_deriv[3][l]*param_in[3] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
		}
			
	fourier_amplitude(frequencies, 
		length,
		amplitude,
		//amplitude_cross_plus,
		"IMRPhenomD",
		parameters);	
			
	}	
	else if (gen_method == "MCMC_IMRPhenomD_single_detect"){
		IMRPhenomD<double> model;
		int dimension = 4;
		source_parameters<double> parameters_in;
		gen_params waveform_params;
		//parameters_in = parameters_in.populate_source_parameters(parameters->mass1, 
		//		parameters->mass2,
		//		parameters->Luminosity_Distance,
		//		parameters->spin1,
		//		parameters->spin2,
		//		parameters->phic,
		//		parameters->tc);
		
		//#################################
		//TESTING 
		//parameters->sky_average=true;
		//#################################
		//parameters->sky_average=false;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		//param_in[0] = parameters_in.DL/MPC_SEC;//MPC -- too big a number for numdiff?
		param_in[0] = parameters_in.chirpmass;//seconds
		param_in[1] = parameters_in.eta;
		param_in[2] = parameters_in.chi_s+parameters_in.chi_a;
		param_in[3] = parameters_in.chi_s-parameters_in.chi_a;

		waveform_params.sky_average=parameters->sky_average;

		for (int i =0; i<dimension; i++){
			for( int j =0;j<dimension;j++){
				param_p[j] = param_in[j] ;
				param_m[j] = param_in[j] ;
			}
			if(param_p[1]>.25) param_p[1] = .25;
			//if(i==0) epsilon = 1e-25;
			//else epsilon = epsilonnaught;
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;
			//param_p[i] = param_in[i]*(1. + epsilon);
			//param_m[i] = param_in[i] *(1- epsilon);

			//model.change_parameter_basis(param_p, param_out, parameters_in.sky_average);
			//param_out[0] = param_p[0];
			param_out[0] = calculate_mass1(param_p[0],param_p[1]);
			param_out[1] = calculate_mass2(param_p[0],param_p[1]);
			param_out[2] = param_p[2];
			param_out[3] = param_p[3];

			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=parameters->Luminosity_Distance;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[2];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[3];

			waveform_params.phic = parameters->phic;
			waveform_params.tc= parameters->tc;
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;


			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_plus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	

			//model.change_parameter_basis(param_m, param_out,parameters_in.sky_average);
			//param_out[0] = param_m[0];
			//param_out[1] = calculate_mass1(param_m[1],param_m[2]);
			//param_out[2] = calculate_mass2(param_m[1],param_m[2]);
			//param_out[3] = param_m[3];
			//param_out[4] = param_m[4];

			param_out[0] = calculate_mass1(param_m[0],param_m[1]);
			param_out[1] = calculate_mass2(param_m[0],param_m[1]);
			param_out[2] = param_m[2];
			param_out[3] = param_m[3];

			waveform_params.mass1 = param_out[0]/MSOL_SEC;
			waveform_params.mass2 = param_out[1]/MSOL_SEC;
			waveform_params.Luminosity_Distance=parameters->Luminosity_Distance;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[2];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[3];

			waveform_params.phic = parameters->phic;
			waveform_params.tc= parameters->tc;
			waveform_params.incl_angle = parameters->incl_angle;
			waveform_params.theta = parameters->theta;
			waveform_params.phi = parameters->phi;
			waveform_params.NSflag = parameters->NSflag;


			fourier_amplitude(frequencies, 
				length,
				amplitude_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			fourier_phase(frequencies, 
				length,
				phase_plus_minus,
				//amplitude_cross_plus,
				"IMRPhenomD",
				&waveform_params);	
			for (int l =0;l<length;l++)
			{
				amplitude_deriv[i][l] = (amplitude_plus_plus[l] -amplitude_plus_minus[l])/(2*epsilon);
				phase_deriv[i][l] = (phase_plus_plus[l] -phase_plus_minus[l])/(2*epsilon);
			}
			//ofstream ampfile
			//for (int l = 0;l<length;l++)
			//{
			//	
			//}
			
		
				
		}
		//Normalize for log factors
		for (int l =0;l<length;l++)
		{
			amplitude_deriv[0][l] = amplitude_deriv[0][l]*param_in[0] ;
			//amplitude_deriv[1][l] = amplitude_deriv[1][l]*param_in[1] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			//phase_deriv[1][l] = phase_deriv[1][l]*param_in[1] ;
		}
		fourier_amplitude(frequencies, 
			length,
			amplitude,
			//amplitude_cross_plus,
			"IMRPhenomD",
			parameters);	
			
	}	
	fourier_amplitude(frequencies, 
		length,
		amplitude,
		//amplitude_cross_plus,
		gen_method,
		parameters);	

	free(amplitude_plus_plus);
	free(amplitude_plus_minus);
	free(amplitude_cross_plus);
	free(amplitude_cross_minus);
	free(phase_plus_plus);
	free(phase_plus_minus);
	free(phase_cross_plus);
	free(phase_cross_minus);
}


/*!\brief Calculates the fisher matrix for the given arguments to within numerical error using automatic differention - slower than the numerical version
 *
 * Build  around  ADOL-C
 */
void fisher_autodiff(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	std::string generation_method, 
	std::string detector, 
	double **output,/**< double [dimension][dimension]*/
	int dimension, 
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
		//for(int i = 0 ; i < length;i++)
		//{
		//	internal_noise[i] = noise[i];
		//}
	}
	else{
		internal_noise = new double[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
		local_noise=true;
	}
		
	//populate derivatives
	double *integrand = (double*)malloc(length*sizeof(double));

	std::complex<double> **response_deriv = new std::complex<double>*[dimension];
	for(int i =0 ;i<dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}
	
	if (parameters->sky_average && generation_method == "IMRPhenomD")
	{
		double **amplitude_deriv = (double **)malloc(dimension*sizeof(**amplitude_deriv));
		for (int i = 0; i<dimension; i++)
			amplitude_deriv[i] = (double *)malloc(length*sizeof(double));
		double **phase_deriv = (double **)malloc(dimension*sizeof(**phase_deriv));
		for (int i = 0; i<dimension; i++)
			phase_deriv[i] = (double *)malloc(length*sizeof(double));
		double *amplitude = (double*)malloc(length*sizeof(double));
		IMRPhenomD<double> model;
		model.fisher_calculation_sky_averaged(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
		for(int j = 0  ;j<dimension; j++){
			for(int i = 0 ; i<length ; i++){
				response_deriv[j][i] = amplitude_deriv[j][i] - std::complex<double>(0,1)*amplitude[i]*phase_deriv[j][i];
			}
		}
		for (int i =0;i<dimension;i++)
		{
			free( amplitude_deriv[i]);
			free( phase_deriv[i]);
		}
		free(amplitude_deriv);
		free(phase_deriv);
		free(amplitude);
	}
	else if (parameters->sky_average &&generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		double **amplitude_deriv = (double **)malloc(dimension*sizeof(**amplitude_deriv));
		for (int i = 0; i<dimension; i++)
			amplitude_deriv[i] = (double *)malloc(length*sizeof(double));
		double **phase_deriv = (double **)malloc(dimension*sizeof(**phase_deriv));
		for (int i = 0; i<dimension; i++)
			phase_deriv[i] = (double *)malloc(length*sizeof(double));
		double *amplitude = (double*)malloc(length*sizeof(double));
		ppE_IMRPhenomD_Inspiral<double> ppemodel;
		ppemodel.fisher_calculation_sky_averaged(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
		for(int j = 0  ;j<dimension; j++){
			for(int i = 0 ; i<length ; i++){
				response_deriv[j][i] = amplitude_deriv[j][i] - std::complex<double>(0,1)*amplitude[i]*phase_deriv[j][i];
			}
		}
		for (int i =0;i<dimension;i++)
		{
			free( amplitude_deriv[i]);
			free( phase_deriv[i]);
		}
		free(amplitude_deriv);
		free(phase_deriv);
		free(amplitude);
	}
	else if (parameters->sky_average &&generation_method == "ppE_IMRPhenomD_IMR")
	{
		double **amplitude_deriv = (double **)malloc(dimension*sizeof(**amplitude_deriv));
		for (int i = 0; i<dimension; i++)
			amplitude_deriv[i] = (double *)malloc(length*sizeof(double));
		double **phase_deriv = (double **)malloc(dimension*sizeof(**phase_deriv));
		for (int i = 0; i<dimension; i++)
			phase_deriv[i] = (double *)malloc(length*sizeof(double));
		double *amplitude = (double*)malloc(length*sizeof(double));
		ppE_IMRPhenomD_IMR<double> ppemodel;
		ppemodel.fisher_calculation_sky_averaged(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
		for(int j = 0  ;j<dimension; j++){
			for(int i = 0 ; i<length ; i++){
				response_deriv[j][i] = amplitude_deriv[j][i] - std::complex<double>(0,1)*amplitude[i]*phase_deriv[j][i];
			}
		}
		for (int i =0;i<dimension;i++)
		{
			free( amplitude_deriv[i]);
			free( phase_deriv[i]);
		}
		free(amplitude_deriv);
		free(phase_deriv);
		free(amplitude);
	}
	else{
		calculate_derivatives_autodiff(frequency,length, dimension,generation_method, parameters, response_deriv, NULL, detector);
	}
	
	//calulate fisher elements
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
					/internal_noise[i]);
			}
			
			output[j][k] = 4*simpsons_sum(
						frequency[1]-frequency[0], length, integrand);	
			output[k][j] = output[j][k];
		}

	}

	for (int j = 0; j<dimension; j ++)
	{

		for (int i =0;i<length;i++){
				integrand[i] = 
					real( (response_deriv[j][i]*std::conj(response_deriv[j][i]))
					/internal_noise[i]);
		}
		output[j][j] = 4*simpsons_sum(
					frequency[1]-frequency[0], length, integrand);	
	}
		


	free(integrand);
	if(local_noise){delete [] internal_noise;}
	for(int i =0 ;i<dimension; i++){
		delete [] response_deriv[i];
	}
	delete [] response_deriv;
}

/*! \brief Calculate the full fisher matrix for a precessing waveform
 *
 * Calculates the fisher without maximizing the extrinsic parameters
 *
 * 15 dimensions: 
 *
 * cos \thetaJ
 *
 * RA
 *
 * DEC
 *
 * DL
 *
 * chirpmass
 *
 * eta
 *
 * spin1
 *
 * spin2
 *
 * theta1
 *
 * theta2
 *
 * phi1
 *
 * phi2
 *
 * phiRef
 *
 * psi
 *
 * all at reference frequency of 20hz
 */
void calculate_derivatives_autodiff(double *frequency,
	int length,
	int dimension,
	std::string generation_method,
	gen_params *parameters,
	std::complex<double> **waveform_deriv,
	int *waveform_tapes,/*<< Waveform tapes -- length=6*/
	std::string detector
	)
{
	//Transform gen_params to double vectors
	double vec_parameters[dimension+1];
	bool log_factors[dimension];
	int boundary_num= boundary_number(generation_method);
	double *freq_boundaries=new double[boundary_num];
	double *grad_freqs=new double[boundary_num];
	unpack_parameters(vec_parameters,log_factors, freq_boundaries,grad_freqs,boundary_num,parameters, generation_method, dimension);
	
	//calculate_derivative tapes
	int tapes[boundary_num];
	for(int i =0; i<boundary_num; i++){
		tapes[i]=i;
		trace_on(tapes[i]);
		adouble avec_parameters[dimension+1];
		avec_parameters[0] <<=grad_freqs[i];
		for(int j = 1; j <= dimension; j++){
			avec_parameters[j]<<=vec_parameters[j];	
		}
		//Repack parameters
		gen_params_base<adouble> a_parameters;
		//############################################
		a_parameters.sky_average = parameters->sky_average;
		a_parameters.f_ref = parameters->f_ref;
		a_parameters.gmst = parameters->gmst;
		a_parameters.NSflag = parameters->NSflag;
		a_parameters.shift_time = false;
		//############################################
		adouble afreq;
		repack_parameters(avec_parameters,&a_parameters, &afreq, generation_method, dimension);
		std::complex<adouble> a_response;
		//HERE
		int status  = fourier_detector_response_equatorial(&afreq, 1, &a_response, detector, generation_method, &a_parameters);

		double response[2];
		real(a_response) >>=  response[0];	
		imag(a_response) >>=  response[1];	

		trace_off();
		
		
	}

	//Evaluate derivative tapes
	int dep = 2;//Output is complex
	int indep = dimension+1;//First element is for frequency
	bool eval = false;//Keep track of when a boundary is hit
	double **jacob = allocate_2D_array(dep,indep);
	for(int k = 0 ;k <length; k++){
		vec_parameters[0]=frequency[k];
		for(int n = 0 ; n<boundary_num; n++){
			if(vec_parameters[0]<freq_boundaries[n]){
				jacobian(tapes[n], dep, indep, vec_parameters, jacob);
				for(int i =0; i<dimension; i++){
					waveform_deriv[i][k] = jacob[0][i+1] 
						+ std::complex<double>(0,1)*jacob[1][i+1];
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
				waveform_deriv[j][i] = (vec_parameters[j+1])*waveform_deriv[j][i];
			}
		}
	}
	deallocate_2D_array(jacob,dep,indep);
	if(!freq_boundaries){
		delete [] freq_boundaries;
	}
	if(!grad_freqs){
		delete [] grad_freqs;
	}

}
int boundary_number(std::string method)
{
	if(method == "IMRPhenomPv2"||
		method=="IMRPhenomD"||
		method=="ppE_IMRPhenomD_IMR"||
		method=="ppE_IMRPhenomD_Inspiral"||
		method=="ppE_IMRPhenomPv2_Inspiral"||
		method=="ppE_IMRPhenomPv2_IMR"){
		return 5;
	}
}
/*! \brief Transforms input gen_params into several base class arrays for use with adolc 
 *
 * This is one of the places where the generation-method/dimension/sky_average specific modifications should go
 *
 * DOES allocate memory for freq_boundaries and grad_freqs that must be deallocated by the use
 */
void unpack_parameters(double *parameters, 
	bool *log_factors, 
	double *freq_boundaries, 
	double *grad_freqs, 
	int boundary_num, 
	gen_params_base<double> *input_params, 
	std::string generation_method, 
	int dimension)
{
	source_parameters<double> s_param;
	s_param = source_parameters<double>::populate_source_parameters(input_params);
	lambda_parameters<double> lambda;
	//incl, RA, DEC, DL, chirpmass, eta, spin1, spin2, theta1, 
	//theta2, phi1, phi2, phiRef, tc, psi
	if(generation_method =="IMRPhenomPv2" && !input_params->sky_average){
		for(int i = 0 ; i<dimension; i++){
			log_factors[i] = false;
		}
		log_factors[3] = true;//Distance
		log_factors[4] = true;//chirpmass
		//###########################################
		//###########################################
		IMRPhenomPv2<double> modelp;
		modelp.assign_lambda_param(&s_param, &lambda);
		modelp.post_merger_variables(&s_param);
		double M = s_param.M;
		double fRD = s_param.fRD;
		double fpeak = modelp.fpeak(&s_param, &lambda);;
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
		grad_freqs[0] = freq_boundaries[0]*.9;
		for(int i = 1 ; i<boundary_num; i++){
			grad_freqs[i] = freq_boundaries[i-1]+(double)(freq_boundaries[i]-freq_boundaries[i-1])/2.;
		}
		//###########################################
		double spin1spher[3];
		double spin2spher[3];
		transform_cart_sph(input_params->spin1, spin1spher);
		transform_cart_sph(input_params->spin2, spin2spher);
		parameters[0]=grad_freqs[0];
		parameters[1]=input_params->incl_angle;
		parameters[2]=input_params->RA;
		parameters[3]=input_params->DEC;
		parameters[4]=s_param.DL;
		parameters[5]=s_param.chirpmass;
		parameters[6]=s_param.eta;
		parameters[7]=spin1spher[0];
		parameters[8]=spin2spher[0];
		parameters[9]=spin1spher[1];
		parameters[10]=spin2spher[1];
		parameters[11]=spin1spher[2];
		parameters[12]=spin2spher[2];
		parameters[13]=input_params->phiRef;
		parameters[14]=input_params->tc;
		parameters[15]=input_params->psi;
		
	}
	//incl, RA,DEC,DL,chirpmass,eta, spin1,spin2,phiRef,tc,psi
	else if(generation_method =="IMRPhenomD" && !input_params->sky_average){
		for(int i = 0 ; i<dimension; i++){
			log_factors[i] = false;
		}
		log_factors[3] = true;//Distance
		log_factors[4] = true;//chirpmass
		//###########################################
		//###########################################
		IMRPhenomD<double> modelp;
		modelp.assign_lambda_param(&s_param, &lambda);
		modelp.post_merger_variables(&s_param);
		double M = s_param.M;
		double fRD = s_param.fRD;
		double fpeak = modelp.fpeak(&s_param, &lambda);;
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
		grad_freqs[0] = freq_boundaries[0]*.9;
		for(int i = 1 ; i<boundary_num; i++){
			grad_freqs[i] = freq_boundaries[i-1]+(double)(freq_boundaries[i]-freq_boundaries[i-1])/2.;
		}
		//###########################################
		double spin1spher[3];
		double spin2spher[3];
		parameters[0]=grad_freqs[0];
		parameters[1]=input_params->incl_angle;
		parameters[2]=input_params->RA;
		parameters[3]=input_params->DEC;
		parameters[4]=s_param.DL;
		parameters[5]=s_param.chirpmass;
		parameters[6]=s_param.eta;
		parameters[7]=s_param.spin1z;
		parameters[8]=s_param.spin2z;
		parameters[9]=input_params->phiRef;
		parameters[10]=input_params->tc;
		parameters[11]=input_params->psi;
		
	}
	
	
}
/*! \brief Repack the parameters from an adouble vector to a gen_params_base<adouble> object and freqeuncy 
 *
 * This is one of the places where the generation-method/dimension/sky_average specific modifications should go
 */
void repack_parameters(adouble *avec_parameters, gen_params_base<adouble> *a_params, adouble *freq, std::string generation_method, int dim)
{
	if(generation_method =="IMRPhenomPv2" && !a_params->sky_average){
		*freq = avec_parameters[0];
		a_params->mass1 = calculate_mass1(avec_parameters[5],avec_parameters[6])/MSOL_SEC;
		a_params->mass2 = calculate_mass2(avec_parameters[5],avec_parameters[6])/MSOL_SEC;
		a_params->Luminosity_Distance = avec_parameters[4]/MPC_SEC;
		a_params->RA = avec_parameters[2];
		a_params->DEC = avec_parameters[3];
		a_params->psi = avec_parameters[15];
		a_params->phiRef = avec_parameters[13];
		a_params->tc = avec_parameters[14];
		adouble spin1sph[3] = {avec_parameters[7],avec_parameters[9],avec_parameters[11]};
		adouble spin2sph[3] = {avec_parameters[8],avec_parameters[10],avec_parameters[12]};
		transform_sph_cart(spin1sph,a_params->spin1);
		transform_sph_cart(spin2sph,a_params->spin2);
		a_params->incl_angle=avec_parameters[1];
	}
	else if(generation_method =="IMRPhenomD" && !a_params->sky_average){
		*freq = avec_parameters[0];
		a_params->mass1 = calculate_mass1(avec_parameters[5],avec_parameters[6])/MSOL_SEC;
		a_params->mass2 = calculate_mass2(avec_parameters[5],avec_parameters[6])/MSOL_SEC;
		a_params->Luminosity_Distance = avec_parameters[4]/MPC_SEC;
		a_params->RA = avec_parameters[2];
		a_params->DEC = avec_parameters[3];
		a_params->psi = avec_parameters[11];
		a_params->phiRef = avec_parameters[9];
		a_params->tc = avec_parameters[10];
		adouble spin1sph[3] = {avec_parameters[7],0,0};
		adouble spin2sph[3] = {avec_parameters[8],0,0};
		transform_sph_cart(spin1sph,a_params->spin1);
		transform_sph_cart(spin2sph,a_params->spin2);
		a_params->incl_angle=avec_parameters[1];
	}

}
