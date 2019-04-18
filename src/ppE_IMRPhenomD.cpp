#include "ppE_IMRPhenomD.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>


/*! \file 
 * File for the implementation of the ppE formalism for testing GR
 *
 * Extends the IMRPhenomD template to include non-GR phase terms
 */

/*!\brief Overloaded method for the inspiral portion of the phase
 */
template<class T>
T ppE_IMRPhenomD_Inspiral<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
	IMRPhenomD<T> model;
	//T piMFcube = powers->MFthird * powers->PIcube;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins(f, param, pn_coeff, lambda,powers);
	T phaseout= gr_ins;
	for(int i = 0; i<param->Nmod; i++)
		//phaseout += pow((param->chirpmass *M_PI * f),param->bppe[i]/3.) * param->betappe[i];
		phaseout =phaseout +  pow_int((PIMFcube),param->bppe[i]) * param->betappe[i];
	//return (gr_ins + pow((param->chirpmass *M_PI * f),param->bppe/3.) * param->betappe);
	return phaseout;

}
template<class T>
T ppE_IMRPhenomD_Inspiral<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda)
{
	IMRPhenomD<T> model;
	T gr_ins = model.Dphase_ins(f, param, pn_coeff,lambda);
	T phaseout= gr_ins;
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  ((T)param->bppe[i]) / 3. * pow( f ,(((T)param->bppe[i])/3.-1.))* pow((param->chirpmass *M_PI ),((T)param->bppe[i])/3.) * param->betappe[i];
	//return gr_ins + param->bppe / 3. * pow( f ,(param->bppe/3.-1))* pow((param->chirpmass *M_PI ),param->bppe/3.) * param->betappe;
	std::cout<<"Dphase test"<<std::endl;	
	return phaseout;

}
//#####################################################################
//Added additional ppE parameter to derivative calculations

template <class T>
void ppE_IMRPhenomD_Inspiral<T>::fisher_calculation(double *frequency, 
			int length, 
			gen_params *parameters,
			double **amplitude_deriv, 
			double **phase_deriv, 
			double *amplitude, 
			int *amp_tapes, 
			int *phase_tapes
			)
{
		ppE_IMRPhenomD_Inspiral<double> modeld;
		int dimension = 7 + parameters->Nmod;
		//populate model
		source_parameters<double> input_params;
		//double spin1vec[3] = {0,0,parameters[3]};
		//double spin2vec[3] = {0,0,parameters[4]};
		//###########################################################################
		//DOUBLE CHECK ORDER OF PARAMETERS
		//input_params = source_parameters<double>::populate_source_parameters(parameters[0],
		//	parameters[1],parameters[2],spin1vec,spin2vec,parameters[5],parameters[6]);
		input_params = source_parameters<double>::populate_source_parameters(parameters->mass1,
			parameters->mass2,parameters->Luminosity_Distance,parameters->spin1,
			parameters->spin2,parameters->phic,parameters->tc);
		//Need the splitting frequency	
		lambda_parameters<double> lambda, *lambda_ptr;
		modeld.assign_lambda_param(&input_params, &lambda);
		modeld.post_merger_variables(&input_params);
		input_params.f1 = 0.014/(input_params.M);
		input_params.f3 = modeld.fpeak(&input_params, &lambda);
		input_params.f1_phase = 0.018/(input_params.M);
		input_params.f2_phase = input_params.fRD/2.;
		input_params.betappe = parameters->betappe;
		input_params.bppe = parameters->bppe;
		input_params.Nmod = parameters->Nmod;
		//###########################################################################
		//populate derivative
		this->construct_amplitude_derivative(frequency, 
				length,
				dimension, 
				amplitude_deriv,
				&input_params,
				amp_tapes
				);
		
		//PROBLEM SECTION
		this->construct_phase_derivative(frequency, 
				length,
				dimension, 
				phase_deriv,
				&input_params,
				phase_tapes
				);
		modeld.construct_amplitude(frequency, length, amplitude, &input_params);
		//LOG FACTORS for A0, chirpmass, and eta)
		
		for (int i = 0;i <length; i++)
		{
			amplitude_deriv[0][i] = (input_params.A0)*amplitude_deriv[0][i];
			amplitude_deriv[3][i] = (input_params.chirpmass)*amplitude_deriv[3][i];
			amplitude_deriv[4][i] = (input_params.eta)*amplitude_deriv[4][i];
			phase_deriv[0][i] = (input_params.A0)*phase_deriv[0][i];
			phase_deriv[3][i] = (input_params.chirpmass)*phase_deriv[3][i];
			phase_deriv[4][i] = (input_params.eta)*phase_deriv[4][i];
		}
		

}
/*! \brief Creates the tapes for derivatives of the amplitude
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void ppE_IMRPhenomD_Inspiral<T>::amplitude_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	ppE_IMRPhenomD_Inspiral<double> modeld;
	//Find splitting frequencies between tapes
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1 = 0.014/(input_params->M);
	input_params->f3 = modeld.fpeak(input_params, &lambda);
	double f2 = (input_params->f1 + input_params->f3)/2.;
	double freqs[3] = {input_params->f1 * .9, f2, input_params->f3 * 1.1};
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double amp_out;

		adouble freq  ;
		adouble parameters[7+input_params->Nmod];
		//adouble amp ;
		ppE_IMRPhenomD_Inspiral<adouble> model;
		adouble new_params[7+input_params->Nmod];

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phic;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;

		for(int j = 0; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params-> betappe[j];
		//parameters[3] = parameters[3]/MSOL_SEC;
		model.change_parameter_basis(parameters, new_params);

		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		source_parameters<adouble> intermediate_params;
		//intermediate_params.betappe = new adouble[input_params->Nmod];	
		//intermediate_params.bppe = new int[input_params->Nmod];	
		//intermediate_params.betappe = (adouble *)malloc(sizeof(adouble)*input_params->Nmod);
		//intermediate_params.betappe[0] = parameters[7];
		//intermediate_params.bppe = (int *)malloc(sizeof(int)*input_params->Nmod);
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		intermediate_params.Nmod = input_params->Nmod;
		intermediate_params = intermediate_params.populate_source_parameters(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6]);
		for(int j = 0; j<input_params->Nmod; j++){
			//intermediate_params.betappe[j] = parameters[7+j];
			//intermediate_params.bppe[j] = input_params->bppe[j];
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		//intermediate_params.betappe = parameters[7];
		//intermediate_params.bppe = input_params->bppe;

		adouble freqs[1] = {freq};
		adouble amp_temp[1];
		model.construct_amplitude(freqs, 1,  amp_temp, &intermediate_params);
		amp_temp[0]>>=amp_out;
		trace_off();
		//delete [] intermediate_params.betappe;
		//delete [] intermediate_params.bppe;
		//free(intermediate_params.betappe);
		//free(intermediate_params.bppe);
	}
	
}

/*! \brief Creates the tapes for derivatives of phase
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void ppE_IMRPhenomD_Inspiral<T>::phase_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	ppE_IMRPhenomD_Inspiral<double> modeld;
	//Find splitting frequencies between tapes
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1_phase = 0.018/(input_params->M);
	input_params->f2_phase = input_params->fRD/2.;
	double f2 = (input_params->f1_phase + input_params->f2_phase)/2.;
	double freqs[3] = {input_params->f1_phase * .9, f2, input_params->f2_phase * 1.1};
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double phase_out;

		adouble freq  ;
		adouble parameters[7+input_params->Nmod];
		adouble phase ;
		ppE_IMRPhenomD_Inspiral<adouble> model;

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phic;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		for (int j =0 ; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params->betappe[j];
		//parameters[7] <<= input_params-> betappe;
		//parameters[3] = parameters[3]/MSOL_SEC;
		adouble new_params[7+input_params->Nmod];
		model.change_parameter_basis(parameters, new_params);
		source_parameters<adouble> intermediate_params;
		//intermediate_params.betappe = new adouble[input_params->Nmod];	
		//intermediate_params.bppe = new int[input_params->Nmod];	
		intermediate_params.Nmod = input_params->Nmod;
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		intermediate_params = intermediate_params.populate_source_parameters(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6]);
		for(int j = 0; j<input_params->Nmod; j++){
			//intermediate_params.betappe[j] = parameters[7+j];
			//intermediate_params.bppe[j] = input_params->bppe[j];
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		//intermediate_params.betappe = parameters[7];
		//intermediate_params.bppe = input_params->bppe;
		adouble freqs[1] = {freq};
		adouble phase_temp[1];
		//model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		phase_temp[0]>>=phase_out;
		trace_off();
		//delete [] intermediate_params.betappe;
		//delete [] intermediate_params.bppe;
	}
	
}
/*! \brief Construct the derivative of the amplitude for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void ppE_IMRPhenomD_Inspiral<T>::construct_amplitude_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **amplitude_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	if(tapes==NULL)
	{
		int tape_temp[3];
		tape_temp[0]=20;
		tape_temp[1]=21;
		tape_temp[2]=22;
		this->amplitude_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phic;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for (int i = 0 ; i<input_params->Nmod; i++){
		evaluate_params[8+i] = input_params-> betappe[i];
	}
	//evaluate_params[8] = input_params-> betappe;
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = 0;
			} 

		}
		else{	
		if(evaluate_params[0]<input_params->f1)
		{
			gradient(tapes[0],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		else if(evaluate_params[0] <input_params->f3)
		{
			gradient(tapes[1],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		else
		{
			gradient(tapes[2],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		}
	}
	
}
/*! \brief Construct the derivative of the phase for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void ppE_IMRPhenomD_Inspiral<T>::construct_phase_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **phase_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	if(tapes==NULL)
	{
		int tape_temp[3];
		tape_temp[0]=23;
		tape_temp[1]=24;
		tape_temp[2]=25;
		this->phase_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phic;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for (int i = 0 ; i<input_params->Nmod; i++){
		evaluate_params[8+i] = input_params-> betappe[i];
	}
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = 0;
			} 

		}
		else{	
		if(evaluate_params[0]<input_params->f1_phase)
		{
			gradient(tapes[0],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else if(evaluate_params[0] <input_params->f2_phase)
		{
			gradient(tapes[1],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else
		{
			gradient(tapes[2],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		}
	}
	
}
//#############################################################################################################################
template<class T>
T ppE_IMRPhenomD_IMR<T>::phase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T gr_mr = model.phase_mr(f, param, lambda);
	T phaseout = gr_mr;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		//phaseout+= pow((param->chirpmass *M_PI * f),param->bppe[i]/3.) * param->betappe[i];
		phaseout+= pow_int((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;
	//return (gr_mr + pow((param->chirpmass *M_PI * f),param->bppe/3.) * param->betappe);

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::Dphase_mr(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T Dgr_mr = model.Dphase_mr(f, param, lambda);
	T phaseout = Dgr_mr;
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,(param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),param->bppe[i]/3.) * param->betappe[i];
	return phaseout;
	//return Dgr_mr + param->bppe / 3. * pow( f ,(param->bppe/3.-1))* pow((param->chirpmass *M_PI ),param->bppe/3.) * param->betappe;

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::Dphase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T Dgr_int = model.Dphase_int(f, param, lambda);
	T phaseout = Dgr_int;
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= param->bppe[i] / 3. * pow( f ,(param->bppe[i]/3.-1))* pow((param->chirpmass *M_PI ),param->bppe[i]/3.) * param->betappe[i];
	return phaseout;
	//return Dgr_int + param->bppe / 3. * pow( f ,(param->bppe/3.-1))* pow((param->chirpmass *M_PI ),param->bppe/3.) * param->betappe;

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T gr_int = model.phase_int(f, param, lambda);
	T phaseout = gr_int;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		//phaseout+= pow((param->chirpmass *M_PI * f),param->bppe[i]/3.) * param->betappe[i];
		phaseout+= pow_int((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;
	//return (gr_int + pow((param->chirpmass *M_PI * f),param->bppe/3.) * param->betappe);

}
//#####################################################################
//Added additional ppE parameter to derivative calculations

template <class T>
void ppE_IMRPhenomD_IMR<T>::fisher_calculation(double *frequency, 
			int length, 
			gen_params *parameters,
			double **amplitude_deriv, 
			double **phase_deriv, 
			double *amplitude, 
			int *amp_tapes, 
			int *phase_tapes
			)
{
		ppE_IMRPhenomD_IMR<double> modeld;
		int dimension = 7+parameters->Nmod;
		//populate model
		source_parameters<double> input_params;
		//double spin1vec[3] = {0,0,parameters[3]};
		//double spin2vec[3] = {0,0,parameters[4]};
		//###########################################################################
		//DOUBLE CHECK ORDER OF PARAMETERS
		//input_params = source_parameters<double>::populate_source_parameters(parameters[0],
		//	parameters[1],parameters[2],spin1vec,spin2vec,parameters[5],parameters[6]);
		input_params = source_parameters<double>::populate_source_parameters(parameters->mass1,
			parameters->mass2,parameters->Luminosity_Distance,parameters->spin1,
			parameters->spin2,parameters->phic,parameters->tc);
		//Need the splitting frequency	
		lambda_parameters<double> lambda, *lambda_ptr;
		modeld.assign_lambda_param(&input_params, &lambda);
		modeld.post_merger_variables(&input_params);
		input_params.f1 = 0.014/(input_params.M);
		input_params.f3 = modeld.fpeak(&input_params, &lambda);
		input_params.f1_phase = 0.018/(input_params.M);
		input_params.f2_phase = input_params.fRD/2.;
		input_params.betappe = parameters->betappe;
		input_params.bppe = parameters->bppe;
		input_params.Nmod = parameters->Nmod;
		//###########################################################################
		//populate derivative
		this->construct_amplitude_derivative(frequency, 
				length,
				dimension, 
				amplitude_deriv,
				&input_params,
				amp_tapes
				);
		
		//PROBLEM SECTION
		this->construct_phase_derivative(frequency, 
				length,
				dimension, 
				phase_deriv,
				&input_params,
				phase_tapes
				);
		modeld.construct_amplitude(frequency, length, amplitude, &input_params);
		//LOG FACTORS for A0, chirpmass, and eta)
		
		for (int i = 0;i <length; i++)
		{
			amplitude_deriv[0][i] = (input_params.A0)*amplitude_deriv[0][i];
			amplitude_deriv[3][i] = (input_params.chirpmass)*amplitude_deriv[3][i];
			amplitude_deriv[4][i] = (input_params.eta)*amplitude_deriv[4][i];
			phase_deriv[0][i] = (input_params.A0)*phase_deriv[0][i];
			phase_deriv[3][i] = (input_params.chirpmass)*phase_deriv[3][i];
			phase_deriv[4][i] = (input_params.eta)*phase_deriv[4][i];
		}
		

}
/*! \brief Creates the tapes for derivatives of the amplitude
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void ppE_IMRPhenomD_IMR<T>::amplitude_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	ppE_IMRPhenomD_IMR<double> modeld;
	//Find splitting frequencies between tapes
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1 = 0.014/(input_params->M);
	input_params->f3 = modeld.fpeak(input_params, &lambda);
	double f2 = (input_params->f1 + input_params->f3)/2.;
	double freqs[3] = {input_params->f1 * .9, f2, input_params->f3 * 1.1};
	int Nmod = input_params->Nmod;
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double amp_out;

		adouble freq  ;
		adouble parameters[7+Nmod];
		//adouble amp ;
		ppE_IMRPhenomD_IMR<adouble> model;
		adouble new_params[7 + Nmod];

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phic;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		for(int j =0 ; i<Nmod ; i++){
			parameters[7+j] <<= input_params-> betappe[j];
		}
		//parameters[3] = parameters[3]/MSOL_SEC;
		model.change_parameter_basis(parameters, new_params);

		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		source_parameters<adouble> intermediate_params;
		intermediate_params.betappe = new adouble[Nmod];	
		intermediate_params.bppe = new int[Nmod];	
		intermediate_params.Nmod = input_params->Nmod;
		intermediate_params = intermediate_params.populate_source_parameters(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6]);
		for(int j = 0; i<Nmod; i++) {
			intermediate_params.betappe[j] = parameters[6+j];
			intermediate_params.bppe[j] = input_params->bppe[j];
		}
			//intermediate_params.betappe = parameters[7];
			//intermediate_params.bppe = input_params->bppe;

		adouble freqs[1] = {freq};
		adouble amp_temp[1];
		model.construct_amplitude(freqs, 1,  amp_temp, &intermediate_params);
		amp_temp[0]>>=amp_out;
		trace_off();
		delete [] intermediate_params.betappe;
		delete [] intermediate_params.bppe;
	}
	
}

/*! \brief Creates the tapes for derivatives of phase
 *
 * For efficiency in long runs of large sets of fishers, the tapes can be precomputed and reused
 */
template <class T>
void ppE_IMRPhenomD_IMR<T>::phase_tape(source_parameters<double> *input_params, /**< source parameters structure of the desired source*/
				int *tape /**<tape ids*/
				)
{
	ppE_IMRPhenomD_IMR<double> modeld;
	//Find splitting frequencies between tapes
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(input_params, &lambda);
	modeld.post_merger_variables(input_params);
	input_params->f1_phase = 0.018/(input_params->M);
	input_params->f2_phase = input_params->fRD/2.;
	double f2 = (input_params->f1_phase + input_params->f2_phase)/2.;
	double freqs[3] = {input_params->f1_phase * .9, f2, input_params->f2_phase * 1.1};
	int Nmod = input_params->Nmod;
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double phase_out;

		adouble freq  ;
		adouble parameters[8+Nmod];
		adouble phase ;
		ppE_IMRPhenomD_IMR<adouble> model;

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phic;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		//parameters[7] <<= input_params-> betappe;
		for(int i =0 ; i<Nmod ; i++){
			parameters[7+i] <<= input_params-> betappe[i];
		}
		//parameters[3] = parameters[3]/MSOL_SEC;
		adouble new_params[8+Nmod];
		model.change_parameter_basis(parameters, new_params);
		source_parameters<adouble> intermediate_params;
		intermediate_params.betappe = new adouble[Nmod];	
		intermediate_params.bppe = new int[Nmod];	
		intermediate_params.Nmod = input_params->Nmod;
		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		intermediate_params = intermediate_params.populate_source_parameters(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6]);
		for (int j = 0 ; j<Nmod; j++){
			intermediate_params.betappe[j] = parameters[7+j];
			intermediate_params.bppe[j] = input_params->bppe[j];
		}
		//intermediate_params.betappe = parameters[7];
		//intermediate_params.bppe = input_params->bppe;
		adouble freqs[1] = {freq};
		adouble phase_temp[1];
		//model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		phase_temp[0]>>=phase_out;
		trace_off();
		delete [] intermediate_params.betappe;
		delete [] intermediate_params.bppe;
	}
	
}
/*! \brief Construct the derivative of the amplitude for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void ppE_IMRPhenomD_IMR<T>::construct_amplitude_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **amplitude_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	if(tapes==NULL)
	{
		int tape_temp[3];
		tape_temp[0]=20;
		tape_temp[1]=21;
		tape_temp[2]=22;
		this->amplitude_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phic;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for(int i = 0 ; i<input_params->Nmod; i++)
		evaluate_params[8+i] = input_params-> betappe[i];
	//evaluate_params[8] = input_params-> betappe;
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = 0;
			} 

		}
		else{	
		if(evaluate_params[0]<input_params->f1)
		{
			gradient(tapes[0],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		else if(evaluate_params[0] <input_params->f3)
		{
			gradient(tapes[1],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		else
		{
			gradient(tapes[2],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				amplitude_derivative[j][i] = grad[j+1];
			} 
		}
		}
	}
	
}
/*! \brief Construct the derivative of the phase for a given source evaluated by the given frequency
 *
 * Order of output: dh/d \theta : \theta \el {A0,tc, phic, chirp mass, eta, symmetric spin, antisymmetric spin}
 * 
 */
template <class T>
void ppE_IMRPhenomD_IMR<T>::construct_phase_derivative(double *frequencies, /**< input array of frequency*/
				int length,/**< length of the frequency array*/
				int dimension,/** < dimension of the fisher*/
				double **phase_derivative,/**< output array for all the derivatives double[dimension][length]*/
				source_parameters<double> *input_params,/**< Source parameters structure for the source*/
				int *tapes /**<int array of tape ids, if NULL, these will be calculated*/
				)
{
	if(tapes==NULL)
	{
		int tape_temp[3];
		tape_temp[0]=23;
		tape_temp[1]=24;
		tape_temp[2]=25;
		this->phase_tape(input_params,tape_temp);
		tapes = tape_temp;
	}
	
	double evaluate_params[dimension +1] ;
	double grad[dimension+1];
	evaluate_params[1] = input_params-> A0;
	evaluate_params[2] = input_params-> phic;
	evaluate_params[3] = input_params-> tc;
	evaluate_params[4] = input_params-> chirpmass;
	evaluate_params[5] = input_params-> eta;
	evaluate_params[6] = input_params-> chi_s;
	evaluate_params[7] = input_params-> chi_a;
	for(int i = 0 ; i<input_params->Nmod;i++)
		evaluate_params[8+i] = input_params-> betappe[i];
	//evaluate_params[8] = input_params-> betappe;
	for (int i=0;i<length; i++)
	{
		evaluate_params[0] = frequencies[i];
		if(evaluate_params[0]>.2/input_params->M){
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = 0;
			} 

		}
		else{	
		if(evaluate_params[0]<input_params->f1_phase)
		{
			gradient(tapes[0],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else if(evaluate_params[0] <input_params->f2_phase)
		{
			gradient(tapes[1],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		else
		{
			gradient(tapes[2],dimension+1, evaluate_params, grad);
			for(int j=0;j<dimension;j++)
			{
				phase_derivative[j][i] = grad[j+1];
			} 
		}
		}
	}
	
}
template class ppE_IMRPhenomD_Inspiral<double>;
template class ppE_IMRPhenomD_Inspiral<adouble>;
template class ppE_IMRPhenomD_IMR<double>;
template class ppE_IMRPhenomD_IMR<adouble>;
