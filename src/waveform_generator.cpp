#include <iostream>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomP.h"
#include "gIMRPhenomD.h"
#include "util.h"
#include <complex>
#include <time.h>
#include <adolc/adouble.h>
using namespace std;

/*!\file
 *File that handles the construction of the (2,2) waveform as described by IMRPhenomD by Khan et. al. 
 * 
 * Builds a waveform for given DETECTOR FRAME parameters
 */

/*!\brief Function to produce the plus/cross polarizations of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 * This puts the responsibility on the user to pass the necessary parameters
 *
 *
 * *NEED TO OUTLINE OPTIONS FOR EACH METHOD IN DEPTH*
 *
 *
 * NEW PHASE OPTIONS for 
 *
 * PHENOMD ONLY:
 *
 * If phic is assigned, the reference frequency and reference phase are IGNORED.
 *
 * If Phic is unassigned, a reference phase AND a reference frequency are looked for.If no options are found, both are set to 0.
 *
 * If tc is assigned, it is used. 
 *
 * If tc is unassigned, the waveform is shifted so the merger happens at 0. 
 *
 * PhenomPv2:
 *
 * PhiRef and f_ref are required, phic is not an option. 
 *
 * tc, if specified, is used with the use of interpolation. If not, tc is set such that coalescence happens at t=0
 */
template<class T>
int fourier_waveform(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			std::complex<T> *waveform_plus, /**< [out] complex array for the output plus polarization waveform*/
			std::complex<T> *waveform_cross, /**< [out] complex array for the output cross polarization waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters/**<structure containing all the source parameters*/
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag1 || NSflag2)
	{
		cout<<"NS waveforms still under develpment - BH only"<<endl;
		return 0;
	}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	//double mass1 = parameters->mass1;
	//double mass2 = parameters->mass2;
	//double Luminosity_Distance = parameters->Luminosity_Distance;
	//double *spin1 = parameters->spin1;
	//double *spin2 = parameters->spin2;
	//double phi_c = parameters->phic;
	//double t_c = parameters->tc;
	source_parameters<T> params;
	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	params = params.populate_source_parameters(parameters);
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.incl_angle = parameters->incl_angle;
	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.sky_average = parameters->sky_average;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.dep_postmerger = parameters->dep_postmerger;
	if(generation_method == "IMRPhenomD")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		IMRPhenomD<T> modeld;
		status = modeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		ppE_IMRPhenomD_Inspiral<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		bool local_spline = false;
		dCS_IMRPhenomD_log<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];

		
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		bool local_spline = false;
		dCS_IMRPhenomD<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];

		
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		EdGB_IMRPhenomD_log<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		EdGB_IMRPhenomD<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		ppE_IMRPhenomD_IMR<T> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i] * std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
	}
	else if(generation_method == "gIMRPhenomD")
	{
		std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		gIMRPhenomD<T> gmodeld;
		params.delta_phi = parameters->delta_phi;
		params.delta_sigma = parameters->delta_sigma;
		params.delta_beta = parameters->delta_beta;
		params.delta_alpha = parameters->delta_alpha;
		params.phii = parameters->phii;
		params.sigmai = parameters->sigmai;
		params.betai = parameters->betai;
		params.alphai = parameters->alphai;
		params.Nmod_phi = parameters->Nmod_phi;
		params.Nmod_sigma = parameters->Nmod_sigma;
		params.Nmod_beta = parameters->Nmod_beta;
		params.Nmod_alpha = parameters->Nmod_alpha;
		status = gmodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<T>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* std::complex<T>(.5,0) *(std::complex<T>(1,0)+ci*ci);
		}
	}
	else if(generation_method == "IMRPhenomPv2")
	{
		//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);

		IMRPhenomPv2<T> modeld;
		//Calculate Waveform
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<T> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempPlus
					+std::complex<T>(sin(2.*params.zeta_polariz),0)*tempCross;
			waveform_cross[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempCross
					-std::complex<T>(sin(2.*params.zeta_polariz),0)*tempPlus;
		}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_Inspiral")
	{
		ppE_IMRPhenomPv2_Inspiral<T> modeld;
		//Initialize Pv2 specific params	

		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		//########################################
		//Check to see if thetaJN was used
		//If not, its calculated
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<T> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempPlus
					+std::complex<T>(sin(2.*params.zeta_polariz),0)*tempCross;
			waveform_cross[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempCross
					-std::complex<T>(sin(2.*params.zeta_polariz),0)*tempPlus;
		}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_IMR")
	{
		ppE_IMRPhenomPv2_IMR<T> modeld;
		//Initialize Pv2 specific params	

		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		//########################################
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		//Calculate Waveform
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<T> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempPlus
					+std::complex<T>(sin(2.*params.zeta_polariz),0)*tempCross;
			waveform_cross[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempCross
					-std::complex<T>(sin(2.*params.zeta_polariz),0)*tempPlus;
		}
	}
	else if(generation_method == "dCS_IMRPhenomPv2")
	{
		//########################################
		//convert betappe for dCS (alpha**2) to the full betappe
		//dCS only supports one modification
		dCS_IMRPhenomD<T> dcs_phenomd;
		params.Nmod = 1;
		params.bppe = new int[1];
		params.bppe[0] = -1;
		params.betappe = new T[1];
		params.betappe[0] = parameters->betappe[0];
		params.betappe[0] = dcs_phenomd.dCS_phase_mod(&params);
		//########################################

		ppE_IMRPhenomPv2_Inspiral<T> model;
		//Initialize Pv2 specific params	

		//########################################
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			model.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			model.PhenomPv2_Param_Transform(&params);
		}
		//Calculate Waveform
		status = model.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<T> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempPlus
					+std::complex<T>(sin(2.*params.zeta_polariz),0)*tempCross;
			waveform_cross[i] = std::complex<T>(cos(2.*params.zeta_polariz),0)*tempCross
					-std::complex<T>(sin(2.*params.zeta_polariz),0)*tempPlus;
		}
		delete [] params.bppe;
		delete [] params.betappe;
	}

	return status ;
}







int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *waveform_plus_real, /**< complex array for the output waveform*/
			double *waveform_plus_imag, /**< complex array for the output waveform*/
			double *waveform_cross_real, /**< complex array for the output waveform*/
			double *waveform_cross_imag, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	//std::complex<double> *waveform_plus = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	//std::complex<double> *waveform_cross = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveform_plus = new std::complex<double>[length];
	std::complex<double> *waveform_cross = new std::complex<double>[length];
	int status = fourier_waveform(frequencies, length, waveform_plus,waveform_cross, generation_method, parameters);
	for(int l = 0 ; l<length;l++){
		waveform_plus_real[l] = std::real(waveform_plus[l]);
		waveform_plus_imag[l] = std::imag(waveform_plus[l]);
		waveform_cross_real[l] = std::real(waveform_cross[l]);
		waveform_cross_imag[l] = std::imag(waveform_cross[l]);
	}
	delete [] waveform_plus;
	delete [] waveform_cross;
	return status;
}





/*!\brief Function to produce the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 *
 */
int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			std::complex<double> *waveform, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag1 || NSflag2)
	{
		cout<<"NS waveforms still under develpment - BH only"<<endl;
		return 0;
	}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	double mass1 = parameters->mass1;
	double mass2 = parameters->mass2;
	double Luminosity_Distance = parameters->Luminosity_Distance;
	double *spin1 = parameters->spin1;
	double *spin2 = parameters->spin2;
	double phi_c = parameters->phic;
	double t_c = parameters->tc;
	source_parameters<double> params;
	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	params = params.populate_source_parameters(parameters);

	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.dep_postmerger = parameters->dep_postmerger;
	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_waveform(frequencies, length, waveform, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		bool local_spline = false;
		dCS_IMRPhenomD_log<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		bool local_spline = false;
		dCS_IMRPhenomD<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		EdGB_IMRPhenomD_log<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		EdGB_IMRPhenomD<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
	}
	//else if(generation_method == "IMRPhenomPv2")
	//{
	//	IMRPhenomPv2<double> modeld;
	//	//Initialize Pv2 specific params	
	//	status = modeld.construct_waveform(frequencies, length, waveform, &params);	
	//}

	return status ;
}

int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *waveform_real, /**< complex array for the output waveform*/
			double *waveform_imag, /**< complex array for the output waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	std::complex<double> *waveform = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	int status = fourier_waveform(frequencies, length, waveform, generation_method, parameters);
	for(int i = 0 ;i<length;i++)
	{
		waveform_real[i] = real(waveform[i]);
		waveform_imag[i] = imag(waveform[i]);
	}
	free(waveform);
	return status;
}

/*!\brief Function to produce the amplitude of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_amplitude(T *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *amplitude, /**< output array for the amplitude*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag1 || NSflag2)
	{
		cout<<"NS waveforms still under develpment - BH only"<<endl;
		return 0;
	}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;
	params = params.populate_source_parameters(parameters);
	//params.f_ref = parameters->f_ref;
	//params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.dep_postmerger = parameters->dep_postmerger;
	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<T> modeld;
		status = modeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		params.betappe = parameters->betappe;
		ppE_IMRPhenomD_Inspiral<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		dCS_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		dCS_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		EdGB_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		EdGB_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		params.betappe = parameters->betappe;
		ppE_IMRPhenomD_IMR<T> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}

	return status ;
}
template int fourier_amplitude<double>(double *, int , double * ,std::string, gen_params_base<double> *);
template int fourier_amplitude<adouble>(adouble *, int , adouble * ,std::string, gen_params_base<adouble> *);
/*!\brief Function to produce the phase of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_phase(T *frequencies, /**<double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *phase, /**<output array for the phase*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag1 || NSflag2)
	{
		cout<<"NS waveforms still under develpment - BH only"<<endl;
		return 0;
	}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	//double mass1 = parameters->mass1;
	//double mass2 = parameters->mass2;
	//double Luminosity_Distance = parameters->Luminosity_Distance;
	//double *spin1 = parameters->spin1;
	//double *spin2 = parameters->spin2;
	//double phi_c = parameters->phic;
	//double t_c = parameters->tc;
	source_parameters<T> params;
	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	params = params.populate_source_parameters(parameters);
	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<T> modeld;
		status = modeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_Inspiral<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		bool local_spline = false;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
		
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		bool local_spline = false;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
		
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_IMR<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for(int i = 0 ; i<length; i++){
				phase[i]*= (T)(-1.);
		}
	}

	return status ;
}
template int fourier_phase<double>(double *, int , double * , std::string, gen_params_base<double> *);
template int fourier_phase<adouble>(adouble *, int , adouble * , std::string,gen_params_base<adouble> *);
/*!\brief Function to produce the phase of the plus and cross mode of a quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 */
template<class T>
int fourier_phase(T *frequencies, /**<double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			T *phase_plus, /**<output array for the phase*/
			T *phase_cross, /**<output array for the phase*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params_base<T> *parameters
			)
{
	int status=1;
	bool NSflag1 = parameters->NSflag1;
	bool NSflag2 = parameters->NSflag2;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag1 || NSflag2)
	{
		cout<<"NS waveforms still under develpment - BH only"<<endl;
		return 0;
	}
	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters<T> params;
	//params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	params = params.populate_source_parameters(parameters);
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.incl_angle = parameters->incl_angle;
	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.sky_average = parameters->sky_average;
	params.dep_postmerger = parameters->dep_postmerger;

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<T> modeld;
		status = modeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_Inspiral<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		bool local_spline = false;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
		
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		bool local_spline = false;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
		
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD_log<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		T temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		for(int i = 0 ; i<length; i++){
			//phase_plus[i]*= (T)(-1.);
			phase_cross[i] = phase_plus[i]+ M_PI/2.;
		}

		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_IMR<T> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase_plus, &params);	
		//for(int i = 0 ; i<length; i++){
		//	phase_plus[i]*= (T)(-1.);
		//	phase_cross[i] = phase_plus[i]+ M_PI/2.;
		//}
	}
	else if(generation_method == "gIMRPhenomD")
	{
		params.delta_phi = parameters->delta_phi;
		params.delta_sigma = parameters->delta_sigma;
		params.delta_beta = parameters->delta_beta;
		params.delta_alpha = parameters->delta_alpha;
		params.phii = parameters->phii;
		params.sigmai = parameters->sigmai;
		params.betai = parameters->betai;
		params.alphai = parameters->alphai;
		params.Nmod_phi = parameters->Nmod_phi;
		params.Nmod_sigma = parameters->Nmod_sigma;
		params.Nmod_beta = parameters->Nmod_beta;
		params.Nmod_alpha = parameters->Nmod_alpha;
		gIMRPhenomD<T> gmodeld;
		status = gmodeld.construct_phase(frequencies, length, phase_plus, &params);	
		//for(int i = 0 ; i<length; i++){
		//	phase_plus[i]*= (T)(-1.);
		//	phase_cross[i] = phase_plus[i]+ M_PI/2.;
		//}
	}
	else if(generation_method == "IMRPhenomPv2")
	{
		//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);

		IMRPhenomPv2<T> modeld;
		//Calculate Waveform
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		T *phase_plus_temp = new T[length];
		T *phase_cross_temp = new T[length];
		status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
		unwrap_array(phase_plus_temp, phase_plus, length);
		unwrap_array(phase_cross_temp, phase_cross, length);
		delete [] phase_plus_temp;
		delete [] phase_cross_temp;
		//for(int i = 0 ; i<length; i++){
		//	phase_plus[i]*= (T)(-1.);
		//	phase_cross[i]*= (T)(-1.);
		//}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_Inspiral")
	{
		//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;

		ppE_IMRPhenomPv2_Inspiral<T> modeld;
		//Calculate Waveform
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		T *phase_plus_temp = new T[length];
		T *phase_cross_temp = new T[length];
		status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
		unwrap_array(phase_plus_temp, phase_plus, length);
		unwrap_array(phase_cross_temp, phase_cross, length);
		delete [] phase_plus_temp;
		delete [] phase_cross_temp;
		//for(int i = 0 ; i<length; i++){
		//	phase_plus[i]*= (T)(-1.);
		//	phase_cross[i]*= (T)(-1.);
		//}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_IMR")
	{
		//std::complex<T> ci = std::complex<T>(cos(params.incl_angle),0);
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;

		ppE_IMRPhenomPv2_IMR<T> modeld;
		//Calculate Waveform
		if((parameters->chip +1)>DOUBLE_COMP_THRESH){
			params.chip = parameters->chip;
			params.spin1z = parameters->spin1[2];
			params.spin2z = parameters->spin2[2];
			params.phip = parameters->phip;
			modeld.PhenomPv2_Param_Transform_reduced(&params);
		}
		else {
			modeld.PhenomPv2_Param_Transform(&params);
		}
		T *phase_plus_temp = new T[length];
		T *phase_cross_temp = new T[length];
		status = modeld.construct_phase(frequencies, length, phase_plus_temp, phase_cross_temp, &params);
		unwrap_array(phase_plus_temp, phase_plus, length);
		unwrap_array(phase_cross_temp, phase_cross, length);
		delete [] phase_plus_temp;
		delete [] phase_cross_temp;
		//for(int i = 0 ; i<length; i++){
		//	phase_plus[i]*= (T)(-1.);
		//	phase_cross[i]*= (T)(-1.);
		//}
	}

	return status ;
}



template int fourier_waveform<double>(double *, int, std::complex<double> *,std::complex<double> *, std::string, gen_params_base<double> *);
template int fourier_waveform<adouble>(adouble *, int, std::complex<adouble> *,std::complex<adouble> *, std::string, gen_params_base<adouble> *);
template int fourier_phase<double>(double *, int, double *,double *, std::string, gen_params_base<double> *);
template int fourier_phase<adouble>(adouble *, int, adouble *,adouble *, std::string, gen_params_base<adouble> *);
