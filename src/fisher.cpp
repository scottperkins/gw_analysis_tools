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
