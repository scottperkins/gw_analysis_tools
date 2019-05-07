#include <fisher.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <math.h>
#include <string>
#include "util.h"
#include "noise_util.h"
#include "IMRPhenomD.h"
#include "ppE_IMRPhenomD.h"
#include "waveform_generator.h"


using namespace std;





/*!\file 
 *
 * All subroutines associated with waveform differentiation and Fisher analysis
 */

/*!\brief Calculates the fisher matrix for the given arguments
 */
void fisher(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	string generation_method, 
	string detector, 
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
	double internal_noise[length];
	if (noise==NULL)
	{
		//double noise[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
	}
	else
		for(int i = 0 ; i < length;i++)
		{
			internal_noise[i] = noise[i];
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

	//PROBLEM 
	//calulate fisher elements
	for (int j=0;j<dimension; j++)
	{
		for (int k = 0; k<j; k++)
		{
			for (int i =0;i<length;i++)
			{
				integrand[i] = 
					real( (amplitude_deriv[j][i]*amplitude_deriv[k][i]
					+amplitude[i]*amplitude[i]*phase_deriv[j][i]*phase_deriv[k][i])
					/internal_noise[i]);
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
					+amplitude[i]*amplitude[i]*phase_deriv[j][i]*phase_deriv[j][i])
					/internal_noise[i]);
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


/* \brief Abstraction layer for handling the case separation for the different waveforms
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
	double epsilon = 1e-5;
	double epsilonnaught = 1e-7;
	double *amplitude_plus_plus = (double *) malloc(sizeof(double)*length);
	double *amplitude_plus_minus = (double *) malloc(sizeof(double)*length);
	double *amplitude_cross_plus = (double *) malloc(sizeof(double)*length);
	double *amplitude_cross_minus = (double *) malloc(sizeof(double)*length);
	double *phase_plus_plus = (double *) malloc(sizeof(double)*length);
	double *phase_plus_minus = (double *) malloc(sizeof(double)*length);
	double *phase_cross_plus = (double *) malloc(sizeof(double)*length);
	double *phase_cross_minus = (double *) malloc(sizeof(double)*length);
	
	if (gen_method == "IMRPhenomD"){
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
			//if(i==0) epsilon = 1e-25;
			//else epsilon = epsilonnaught;
			param_p[i] = param_in[i] + epsilon;
			param_m[i] = param_in[i] - epsilon;
			//param_p[i] = param_in[i]*(1. + epsilon);
			//param_m[i] = param_in[i] *(1- epsilon);

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
			amplitude_deriv[4][l] = amplitude_deriv[4][l]*param_in[4] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[3][l] = phase_deriv[3][l]*param_in[3] ;
			phase_deriv[4][l] = phase_deriv[4][l]*param_in[4] ;
		}
			
	}	
	else if (gen_method == "MCMC_IMRPhenomD_single_detect"){
		IMRPhenomD<double> model;
		int dimension = 5;
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
		parameters->sky_average=false;
		parameters_in = parameters_in.populate_source_parameters(parameters); 
		double param_p[dimension] ;
		double param_m[dimension] ;
		double param_in[dimension] ;
		double param_out[dimension] ;
		param_in[0] = parameters_in.DL;//seconds
		param_in[1] = parameters_in.chirpmass;//seconds
		param_in[2] = parameters_in.eta;
		param_in[3] = parameters_in.chi_s+parameters_in.chi_a;
		param_in[4] = parameters_in.chi_s-parameters_in.chi_a;

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

			//model.change_parameter_basis(param_p, param_out, parameters_in.sky_average);
			param_out[0] = param_p[0];
			param_out[1] = calculate_mass1(param_p[1],param_p[2]);
			param_out[2] = calculate_mass2(param_p[1],param_p[2]);
			param_out[3] = param_p[3];
			param_out[4] = param_p[4];
			waveform_params.mass1 = param_out[1]/MSOL_SEC;
			waveform_params.mass2 = param_out[2]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[0]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
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
			param_out[0] = param_m[0];
			param_out[1] = calculate_mass1(param_m[1],param_m[2]);
			param_out[2] = calculate_mass2(param_m[1],param_m[2]);
			param_out[3] = param_m[3];
			param_out[4] = param_m[4];
			waveform_params.mass1 = param_out[1]/MSOL_SEC;
			waveform_params.mass2 = param_out[2]/MSOL_SEC;
			waveform_params.Luminosity_Distance=param_out[0]/MPC_SEC;
			waveform_params.spin1[0]=0;
			waveform_params.spin1[1]=0;
			waveform_params.spin1[2]=param_out[3];
			waveform_params.spin2[0]=0;
			waveform_params.spin2[1]=0;
			waveform_params.spin2[2]=param_out[4];
			waveform_params.phic = 0;
			waveform_params.tc= 0;
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
			amplitude_deriv[1][l] = amplitude_deriv[1][l]*param_in[1] ;
			//amplitude_deriv[2][l] = amplitude_deriv[2][l]*param_in[2] ;
			phase_deriv[0][l] = phase_deriv[0][l]*param_in[0] ;
			phase_deriv[1][l] = phase_deriv[1][l]*param_in[1] ;
			//phase_deriv[2][l] = phase_deriv[2][l]*param_in[2] ;
		}
			
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
 */
void fisher_autodiff(double *frequency, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	string generation_method, 
	string detector, 
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
	double internal_noise[length];
	if (noise==NULL)
	{
		//double noise[length];
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
	}
	else
		for(int i = 0 ; i < length;i++)
		{
			internal_noise[i] = noise[i];
		}
		
	//populate derivatives
	double **amplitude_deriv = (double **)malloc(dimension*sizeof(**amplitude_deriv));
	for (int i = 0; i<dimension; i++)
		amplitude_deriv[i] = (double *)malloc(length*sizeof(double));
	double **phase_deriv = (double **)malloc(dimension*sizeof(**phase_deriv));
	for (int i = 0; i<dimension; i++)
		phase_deriv[i] = (double *)malloc(length*sizeof(double));
	double *amplitude = (double*)malloc(length*sizeof(double));
	double *integrand = (double*)malloc(length*sizeof(double));
	
	if (generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> model;
		model.fisher_calculation(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
	}
	else if (generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodel;
		ppemodel.fisher_calculation(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
	}
	else if (generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodel;
		ppemodel.fisher_calculation(frequency, 
			length, 
			parameters,
			amplitude_deriv, 
			phase_deriv, 
			amplitude, 
			amp_tapes, 
			phase_tapes
			);
	}

	//PROBLEM 
	//calulate fisher elements
	for (int j=0;j<dimension; j++)
	{
		for (int k = 0; k<j; k++)
		{
			for (int i =0;i<length;i++)
			{
				integrand[i] = 
					real( (amplitude_deriv[j][i]*amplitude_deriv[k][i]
					+amplitude[i]*amplitude[i]*phase_deriv[j][i]*phase_deriv[k][i])
					/internal_noise[i]);
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
					+amplitude[i]*amplitude[i]*phase_deriv[j][i]*phase_deriv[j][i])
					/internal_noise[i]);
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



//############################################################################
//outdated

//void intialize_tape()
//{
//	int dimension = 7;
//	double A0 = 3.13e-22; 
//	double tc = 0;
//	double phic = 0;
//	double chirpmass = 6e-05;
//	double symm = .2222222;
//	double chi_s = .55;
//	double chi_a = -.35;
//	double params[dimension] = {A0,tc,phic,chirpmass,symm,chi_s,chi_a};
//	trace_on(1);
//	
//	deriv_params[0] <<= frequencies[2];
//	for ( int i =0; i<dimension; i ++)
//		deriv_params[i+1]<<= params[i];
//	
//	f = deriv_params[0];
//	for ( int i =0; i<dimension; i ++)
//		intermediate_params[i]= deriv_params[i+1];
//	y = amp_fun(f, intermediate_params);
//	y>>= out;
//	delete[] deriv_params;
//	delete[] intermediate_params;
//	trace_off();
//}
//void amplitude_derivative(double *frequencies, 	
//			int length, 
//			double **amp_derivative,
//			int dimension,
//			double *params,
//			adouble (*amp_fun)(adouble f, adouble *parameters)
//			)
//{
//	adouble *intermediate_params = new adouble[dimension];
//	adouble *deriv_params = new adouble[dimension+1]; 
//	adouble f;
//	adouble y ;
//	double out ;
//
//	trace_on(1);
//
//	deriv_params[0] <<= frequencies[2];
//	for ( int i =0; i<dimension; i ++)
//		deriv_params[i+1]<<= params[i];
//
//	f = deriv_params[0];
//	for ( int i =0; i<dimension; i ++)
//		intermediate_params[i]= deriv_params[i+1];
//	y = amp_fun(f, intermediate_params);
//	y>>= out;
//	delete[] deriv_params;
//	delete[] intermediate_params;
//	trace_off(1);
//
//	double evaluate_params[dimension +1] ;
//	//double amp_deriv[dimension+1];
//	double **amp_deriv = (double**)malloc(sizeof(**amp_deriv));
//	amp_deriv[0] = (double *)malloc(dimension* sizeof(double));
//	for ( int k =0; k<dimension; k ++)
//		evaluate_params[k+1]= params[k];
//	for (int i=0;i<length; i++)
//	{
//		evaluate_params[0] = frequencies[i];
//		//gradient(1,dimension+1,evaluate_params , amp_deriv);
//		jacobian(1,1,dimension+1,evaluate_params , amp_deriv);
//		for(int j=0;j<dimension;j++)
//		{
//			//amp_derivative[j][i] = amp_deriv[j+1];
//			amp_derivative[j][i] = amp_deriv[0][j+1];
//		} 
//	}
//	
//	
//}
//
//void waveform_derivative(double *frequencies, 	
//			int length, 
//			std::complex<double> **waveform_derivative,
//			int dimension,
//			double *params,
//			std::complex<adouble> (*waveform)(adouble f, adouble *parameters)
//			)
//{
//	adouble *intermediate_params = new adouble[dimension];
//	//double waveform_deriv[2][dimension];
//	double **jac = (double **) malloc(dimension*sizeof(**jac));
//	for (int i =0; i<dimension; i++)
//		jac[i] = (double *)malloc(2*sizeof(double));
//	adouble *deriv_params = new adouble[dimension+1]; 
//	adouble f;
//	std::complex<adouble> y ;
//	double reout ;
//	double imout ;
//
//	trace_on(1);
//
//	deriv_params[0] <<= frequencies[2];
//	for ( int i =0; i<dimension; i ++)
//		deriv_params[i+1]<<= params[i];
//
//	f = deriv_params[0];
//	for ( int i =0; i<dimension; i ++)
//		intermediate_params[i]= deriv_params[i+1];
//	y = waveform(f, intermediate_params);
//	std::real(y)>>= reout;
//	std::imag(y)>>= imout;
//	delete[] deriv_params;
//	delete[] intermediate_params;
//	trace_off();
//
//	double evaluate_params[dimension +1] ;
//	for ( int k =0; k<dimension; k ++)
//		evaluate_params[k+1]= params[k];
//	for (int i=0;i<length; i++)
//	{
//		evaluate_params[0] = frequencies[i];
//		jacobian(1,2,dimension+1,evaluate_params , jac);
//		for(int j=0;j<dimension;j++)
//		{
//			waveform_derivative[j][i] = std::complex<double>(jac[0][j+1],jac[1][j+1]);
//		} 
//	}
//	free(jac);
//	
//	
//}
//adouble IMRPhenomD_amp_parameter_tranform(adouble f, /**< frequency in Hz*/
//					adouble *parameters/**< adouble array of source parameters - [A0,tc,phic, chirpmass, symmetric mass ratio, chi_s, chi_a]*/
//				)
//{
//	//Parameter transformation
//	adouble m1, m2, DL, spin1,spin2;
//	m1 = calculate_mass1(parameters[3],parameters[4]);
//	m2 = calculate_mass2(parameters[3],parameters[4]);
//	adouble chirpmass_sec = parameters[3]*MSOL_SEC;
//	DL =1/(parameters[0] * sqrt(30/M_PI) /(chirpmass_sec*chirpmass_sec)/
//			pow(M_PI*chirpmass_sec,-7./6))/MPC_SEC  ;
//	spin1 = parameters[5]+parameters[6];
//	spin2 = parameters[5]-parameters[6];
//	adouble spinvec1[3] ={0.,0.,spin1};
//	adouble spinvec2[3] ={0.,0.,spin2};
//
//	source_parameters<adouble> s_params;
//	
//	s_params = s_params.populate_source_parameters(m1, m2, DL, spinvec1, spinvec2, parameters[2],parameters[1]);
//
//	IMRPhenomD<adouble> model;	
//	return model.construct_amplitude(f,&s_params);
//	
//}
//
//std::complex<adouble> IMRPhenomD_waveform_parameter_tranform(adouble f, /**< frequency in Hz*/
//					adouble *parameters/**< adouble array of source parameters - [A0,tc,phic, chirpmass, symmetric mass ratio, chi_s, chi_a]*/
//				)
//{
//	//Parameter transformation
//	adouble m1, m2, DL, spin1,spin2;
//	m1 = calculate_mass1(parameters[3],parameters[4]);
//	m2 = calculate_mass2(parameters[3],parameters[4]);
//	adouble chirpmass_sec = parameters[3]*MSOL_SEC;
//	DL =1/(parameters[0] * sqrt(30/M_PI) /(chirpmass_sec*chirpmass_sec)/
//			pow(M_PI*chirpmass_sec,-7./6))/MPC_SEC  ;
//	spin1 = parameters[5]+parameters[6];
//	spin2 = parameters[5]-parameters[6];
//	adouble spinvec1[3] ={0.,0.,spin1};
//	adouble spinvec2[3] ={0.,0.,spin2};
//
//	source_parameters<adouble> s_params;
//	
//	s_params = s_params.populate_source_parameters(m1, m2, DL, spinvec1, spinvec2, parameters[2],parameters[1]);
//
//	IMRPhenomD<adouble> model;	
//	return model.construct_waveform(f,&s_params);
//	
//}
