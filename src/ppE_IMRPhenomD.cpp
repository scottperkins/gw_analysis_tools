#include "ppE_IMRPhenomD.h"
#include <math.h>
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "util.h"

/*! \file 
 * File for the implementation of the ppE formalism for testing GR
 *
 * Extends the IMRPhenomD template to include non-GR phase terms
 *
 * Supported waveforms: ppE Inspiral, ppE IMR, dCS, EdGB
 */

template<class T>
T dCS_IMRPhenomD<T>::dCS_phase_mod( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T redshiftedM = M/(1.+Z);
	T phase_mod = (param->betappe[0]);
	T out =  16.*M_PI*phase_mod/(pow_int(redshiftedM,4)) *this->dCS_phase_factor(param);
	return out;
} 

template<class T>
T dCS_IMRPhenomD<T>::dCS_phase_factor(source_parameters<T> *param)
{
	T g=0;
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T coeff1 = -5./8192.;
	T coeff2 = 15075./114688.;
	T m1 = calculate_mass1(chirpmass,eta);
	T m2 = calculate_mass2(chirpmass,eta);
	T m = m1+m2;
	T chi1 = param->chi_s+param->chi_a;
	T chi2 = param->chi_s-param->chi_a;
	if(param->NSflag1 ){
		chi1 = 0;
	}
	if(param->NSflag2 ){
		chi2 = 0;
	}
	T s1temp = 2.+2.*pow_int(chi1,4) - 2.*sqrt((1.-chi1*chi1)) - chi1*chi1 * ((3. - 2.*sqrt(1.-chi1*chi1)));
	T s2temp = 2.+2.*pow_int(chi2,4) - 2.*sqrt((1.-chi2*chi2)) - chi2*chi2 * ((3. - 2.*sqrt(1.-chi2*chi2)));
	chi1 +=1e-10;
	chi2 +=1e-10;
	T s1  = s1temp/(2.*chi1*chi1*chi1);
	T s2  = s2temp/(2.*chi2*chi2*chi2);
	//Neutron stars don't source scalar charge
	if(param->NSflag1){s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	g+=coeff1/(pow(eta,14./5.)) * pow((m1*s2 - m2 * s1),2.)/(m*m);
	g+=coeff2/(pow(eta,14./5.)) * (m2*m2* chi1*chi1 - 350./201. * m1*m2*chi1*chi2 + m1*m1 * chi2*chi2)/(m*m);
	return g;
}

template<class T>
int dCS_IMRPhenomD<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
	params->betappe[0] = this->dCS_phase_mod(params);
	return ppE_IMRPhenomD_Inspiral<T>::construct_waveform(frequencies, length, waveform, params);
}

template<class T>
int dCS_IMRPhenomD<T>::construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params)
{
	return ppE_IMRPhenomD_Inspiral<T>::construct_amplitude(frequencies, length, amplitude, params);
}

template<class T>
int dCS_IMRPhenomD<T>::construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params)
{
	params->betappe[0] = this->dCS_phase_mod(params);
	return ppE_IMRPhenomD_Inspiral<T>::construct_phase(frequencies, length, phase, params);
}


//########################################################################
//####################################################################
template<class T>
T EdGB_IMRPhenomD<T>::EdGB_phase_mod( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T redshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	return 16.*M_PI*phase_mod/(pow_int(redshiftedM,4)) * this->EdGB_phase_factor(param);
} 

template<class T>
T EdGB_IMRPhenomD<T>::EdGB_phase_factor( source_parameters<T> *param)
{
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T m1 = calculate_mass1(chirpmass, eta);
	T m2 = calculate_mass2(chirpmass, eta);
	T chi1 = param->chi_s + param->chi_a;
	T chi2 = param->chi_s - param->chi_a;
	if(param->NSflag1 ){
		chi1 = 0;
	}
	if(param->NSflag2 ){
		chi2 = 0;
	}
	T temp1 = 2.*(sqrt(1.-chi1*chi1) - 1. + chi1*chi1);
	T temp2 = 2.*(sqrt(1.-chi2*chi2) - 1. + chi2*chi2);
	chi1 += 1.e-10;
	chi2 += 1.e-10;
	T s1 = temp1/(chi1*chi1);
	T s2 = temp2/(chi2*chi2);
	if(param->NSflag1){s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	return (-5./7168.)* pow_int((m1*m1 * s2 - m2*m2 * s1),2) / (pow_int(M,4) * pow(eta,(18./5)));
} 
template<class T>
int EdGB_IMRPhenomD<T>::construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params)
{
	params->betappe[0] = this->EdGB_phase_mod(params);
	return ppE_IMRPhenomD_Inspiral<T>::construct_waveform(frequencies, length, waveform, params);
}

template<class T>
int EdGB_IMRPhenomD<T>::construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params)
{
	return ppE_IMRPhenomD_Inspiral<T>::construct_amplitude(frequencies, length, amplitude, params);
}

template<class T>
int EdGB_IMRPhenomD<T>::construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params)
{
	params->betappe[0] = this->EdGB_phase_mod(params);
	return ppE_IMRPhenomD_Inspiral<T>::construct_phase(frequencies, length, phase, params);
}
//####################################################################


/*!\brief Overloaded method for the inspiral portion of the phase
 */
template<class T>
T ppE_IMRPhenomD_Inspiral<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, 
		lambda_parameters<T> *lambda, useful_powers<T> *powers)
{
	IMRPhenomD<T> model;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	T gr_ins = model.phase_ins(f, param, pn_coeff, lambda,powers);
	T phaseout= gr_ins;
	
	for(int i = 0; i<param->Nmod; i++)
		phaseout =phaseout +  pow_int((PIMFcube),(int)param->bppe[i]) * param->betappe[i];
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
	return phaseout;

}
//#####################################################################
//Added additional ppE parameter to derivative calculations

template <class T>
void ppE_IMRPhenomD_Inspiral<T>::fisher_calculation_sky_averaged(double *frequency, 
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
	//###########################################################################
	input_params = source_parameters<double>::populate_source_parameters(parameters);
	input_params.shift_time = parameters->shift_time;
	//Need the splitting frequency	
	lambda_parameters<double> lambda, *lambda_ptr;
	modeld.assign_lambda_param(&input_params, &lambda);
	modeld.post_merger_variables(&input_params);
	input_params.f1 = 0.014/(input_params.M);
	input_params.f3 = modeld.fpeak(&input_params, &lambda);
	input_params.f1_phase = 0.018/(input_params.M);
	input_params.f2_phase = input_params.fRD/2.;
	input_params.Nmod = parameters->Nmod;
	input_params.bppe = parameters->bppe;
	input_params.betappe = parameters->betappe;
	//###########################################################################
	//populate derivative
	this->construct_amplitude_derivative(frequency, 
			length,
			dimension, 
			amplitude_deriv,
			&input_params,
			amp_tapes
			);
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
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;

		for(int j = 0; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params-> betappe[j];
		model.change_parameter_basis(parameters, new_params,true);

		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		source_parameters<adouble> intermediate_params;
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6], input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;
		for(int j = 0; j<input_params->Nmod; j++){
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		intermediate_params.Nmod = input_params->Nmod;

		adouble freqs[1] = {freq};
		adouble amp_temp[1];
		model.construct_amplitude(freqs, 1,  amp_temp, &intermediate_params);
		amp_temp[0]>>=amp_out;
		trace_off();
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
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		for (int j =0 ; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params->betappe[j];
		adouble new_params[7+input_params->Nmod];
		model.change_parameter_basis(parameters, new_params,true);
		source_parameters<adouble> intermediate_params;
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6],input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;
		for(int j = 0; j<input_params->Nmod; j++){
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		intermediate_params.Nmod = input_params->Nmod;
		adouble freqs[1] = {freq};
		adouble phase_temp[1];
		model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		phase_temp[0]>>=phase_out;
		trace_off();
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
	if(!tapes)
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
	evaluate_params[2] = input_params-> phiRef;
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
	if(!tapes)
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
	evaluate_params[2] = input_params-> phiRef;
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
		phaseout+= pow_int((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;

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

}
template<class T>
T ppE_IMRPhenomD_IMR<T>::phase_int(T f, source_parameters<T> *param, lambda_parameters<T> *lambda)
{
	ppE_IMRPhenomD_Inspiral<T> model;
	T gr_int = model.phase_int(f, param, lambda);
	T phaseout = gr_int;
	T PIMFcube = pow(M_PI * param->chirpmass * f, 1./3.);
	for (int i = 0 ; i<param->Nmod; i++)
		phaseout+= pow_int((PIMFcube),param->bppe[i]) * param->betappe[i];
	return phaseout;

}
//#####################################################################
//Added additional ppE parameter to derivative calculations

template <class T>
void ppE_IMRPhenomD_IMR<T>::fisher_calculation_sky_averaged(double *frequency, 
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
		input_params = source_parameters<double>::populate_source_parameters_old(parameters->mass1,
			parameters->mass2,parameters->Luminosity_Distance,parameters->spin1,
			parameters->spin2,parameters->phiRef,parameters->tc,parameters->sky_average);
		input_params.shift_time = parameters->shift_time;
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
	
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double amp_out;

		adouble freq  ;
		adouble parameters[7+input_params->Nmod];
		//adouble amp ;
		ppE_IMRPhenomD_IMR<adouble> model;
		adouble new_params[7+input_params->Nmod];

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;

		for(int j = 0; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params-> betappe[j];
		model.change_parameter_basis(parameters, new_params,true);

		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		source_parameters<adouble> intermediate_params;
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6], input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;
		for(int j = 0; j<input_params->Nmod; j++){
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		intermediate_params.Nmod = input_params->Nmod;

		adouble freqs[1] = {freq};
		adouble amp_temp[1];
		model.construct_amplitude(freqs, 1,  amp_temp, &intermediate_params);
		amp_temp[0]>>=amp_out;
		trace_off();
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
	for (int i =0; i<3;i ++)	
	{
		trace_on(tape[i]);
		double phase_out;

		adouble freq  ;
		adouble parameters[7+input_params->Nmod];
		adouble phase ;
		ppE_IMRPhenomD_IMR<adouble> model;

		freq <<= freqs[i];
		parameters[0] <<= input_params-> A0;
		parameters[1] <<= input_params-> phiRef;
		parameters[2] <<= input_params-> tc;
		parameters[3] <<= input_params-> chirpmass;
		parameters[4] <<= input_params-> eta;
		parameters[5] <<= input_params-> chi_s;
		parameters[6] <<= input_params-> chi_a;
		for (int j =0 ; j<input_params->Nmod; j++)
			parameters[7+j] <<= input_params->betappe[j];
		adouble new_params[7+input_params->Nmod];
		model.change_parameter_basis(parameters, new_params,true);
		source_parameters<adouble> intermediate_params;
		adouble betappe[input_params->Nmod];
		int bppe[input_params->Nmod];
		adouble spin1vec[3] = {0,0,new_params[3]};
		adouble spin2vec[3] = {0,0,new_params[4]};
		intermediate_params = intermediate_params.populate_source_parameters_old(new_params[0]/MSOL_SEC,
				new_params[1]/MSOL_SEC,new_params[2]/MPC_SEC,spin1vec,spin2vec,new_params[5],
				new_params[6],input_params->sky_average);
		intermediate_params.shift_time = input_params->shift_time;
		for(int j = 0; j<input_params->Nmod; j++){
			betappe[j] = parameters[7+j];
			bppe[j] = input_params->bppe[j];
		}
		intermediate_params.betappe = betappe;
		intermediate_params.bppe = bppe;
		intermediate_params.Nmod = input_params->Nmod;
		adouble freqs[1] = {freq};
		adouble phase_temp[1];
		model.construct_phase(freqs, 1,  phase_temp, &intermediate_params);
		phase_temp[0]>>=phase_out;
		trace_off();
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
	if(!tapes)
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
	evaluate_params[2] = input_params-> phiRef;
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
	if(!tapes)
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
	evaluate_params[2] = input_params-> phiRef;
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
template class dCS_IMRPhenomD<double>;
template class dCS_IMRPhenomD<adouble>;
template class EdGB_IMRPhenomD<double>;
template class EdGB_IMRPhenomD<adouble>;
