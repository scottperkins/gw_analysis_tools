#include <iostream>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
#include "ppE_IMRPhenomP.h"
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
int fourier_waveform(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			std::complex<double> *waveform_plus, /**< [out] complex array for the output plus polarization waveform*/
			std::complex<double> *waveform_cross, /**< [out] complex array for the output cross polarization waveform*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters/**<structure containing all the source parameters*/
			)
{
	int status=1;
	bool NSflag = parameters->NSflag;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag)
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
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.incl_angle = parameters->incl_angle;
	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	double ci = cos(params.incl_angle);
	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i] * .5 *(1+ci*ci);
		}
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* .5 *(1+ci*ci);
		}
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

		
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* .5 *(1+ci*ci);
		}
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

		
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* .5 *(1+ci*ci);
		}
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
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* .5 *(1+ci*ci);
		}
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
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i]* .5 *(1+ci*ci);
		}
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++){
			waveform_cross[i] = ci*std::complex<double>(0,1) * waveform_plus[i];
			waveform_plus[i] = waveform_plus[i] * .5 *(1+ci*ci);
		}
	}
	else if(generation_method == "IMRPhenomPv2")
	{
		IMRPhenomPv2<double> modeld;
		//Initialize Pv2 specific params	

		params.thetaJN = parameters->thetaJN;
		params.alpha0 = parameters->alpha0;
		params.zeta_polariz = parameters->zeta_polariz;
		params.phi_aligned = parameters->phi_aligned;
		params.chil = parameters->chil;
		params.chip = parameters->chip;
		//Check to see if thetaJN was used
		//If not, its calculated
		if(params.thetaJN==-1 )
			//Compute transform with spins and L inclination
			modeld.PhenomPv2_Param_Transform(&params);
		else {
			//compute transform with spins and J inclination
			//Not supported at the moment
			modeld.PhenomPv2_Param_Transform_J(&params);
		}
		//Calculate Waveform
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<double> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = cos(2.*params.zeta_polariz)*tempPlus
					+sin(2.*params.zeta_polariz)*tempCross;
			waveform_cross[i] = (2.*params.zeta_polariz)*tempCross
					-sin(2.*params.zeta_polariz)*tempPlus;
		}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_Inspiral")
	{
		ppE_IMRPhenomPv2_Inspiral<double> modeld;
		//Initialize Pv2 specific params	

		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		//########################################
		params.thetaJN = parameters->thetaJN;
		params.alpha0 = parameters->alpha0;
		params.zeta_polariz = parameters->zeta_polariz;
		params.phi_aligned = parameters->phi_aligned;
		params.chil = parameters->chil;
		params.chip = parameters->chip;
		//Check to see if thetaJN was used
		//If not, its calculated
		if(params.thetaJN==-1 )
			//Compute transform with spins and L inclination
			modeld.PhenomPv2_Param_Transform(&params);
		else {
			//compute transform with spins and J inclination
			//Not supported at the moment
			modeld.PhenomPv2_Param_Transform_J(&params);
		}
		//Calculate Waveform
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<double> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = cos(2.*params.zeta_polariz)*tempPlus
					+sin(2.*params.zeta_polariz)*tempCross;
			waveform_cross[i] = (2.*params.zeta_polariz)*tempCross
					-sin(2.*params.zeta_polariz)*tempPlus;
		}
	}
	else if(generation_method == "ppE_IMRPhenomPv2_IMR")
	{
		ppE_IMRPhenomPv2_IMR<double> modeld;
		//Initialize Pv2 specific params	

		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		//########################################
		params.thetaJN = parameters->thetaJN;
		params.alpha0 = parameters->alpha0;
		params.zeta_polariz = parameters->zeta_polariz;
		params.phi_aligned = parameters->phi_aligned;
		params.chil = parameters->chil;
		params.chip = parameters->chip;
		//Check to see if thetaJN was used
		//If not, its calculated
		if(params.thetaJN==-1 )
			//Compute transform with spins and L inclination
			modeld.PhenomPv2_Param_Transform(&params);
		else {
			//compute transform with spins and J inclination
			//Not supported at the moment
			modeld.PhenomPv2_Param_Transform_J(&params);
		}
		//Calculate Waveform
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<double> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = cos(2.*params.zeta_polariz)*tempPlus
					+sin(2.*params.zeta_polariz)*tempCross;
			waveform_cross[i] = (2.*params.zeta_polariz)*tempCross
					-sin(2.*params.zeta_polariz)*tempPlus;
		}
	}
	else if(generation_method == "dCS_IMRPhenomPv2")
	{
		//########################################
		//convert betappe for dCS (alpha**2) to the full betappe
		//dCS only supports one modification
		dCS_IMRPhenomD<double> dcs_phenomd;
		params.Nmod = 1;
		params.bppe = new int[1];
		params.bppe[0] = -1;
		params.betappe = new double[1];
		params.betappe[0] = parameters->betappe[0];
		params.betappe[0] = dcs_phenomd.dCS_phase_mod(&params);
		//########################################

		ppE_IMRPhenomPv2_Inspiral<double> model;
		//Initialize Pv2 specific params	

		//########################################
		params.thetaJN = parameters->thetaJN;
		params.alpha0 = parameters->alpha0;
		params.zeta_polariz = parameters->zeta_polariz;
		params.phi_aligned = parameters->phi_aligned;
		params.chil = parameters->chil;
		params.chip = parameters->chip;
		//Check to see if thetaJN was used
		//If not, its calculated
		if(params.thetaJN==-1 )
			//Compute transform with spins and L inclination
			model.PhenomPv2_Param_Transform(&params);
		else {
			//compute transform with spins and J inclination
			//Not supported at the moment
			model.PhenomPv2_Param_Transform_J(&params);
		}
		//Calculate Waveform
		status = model.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
		std::complex<double> tempPlus,tempCross;
		for (int i =0;i < length; i++)
		{
			tempPlus = waveform_plus[i];	
			tempCross = waveform_cross[i];	
			waveform_plus[i] = cos(2.*params.zeta_polariz)*tempPlus
					+sin(2.*params.zeta_polariz)*tempCross;
			waveform_cross[i] = (2.*params.zeta_polariz)*tempCross
					-sin(2.*params.zeta_polariz)*tempPlus;
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
	std::complex<double> *waveform_plus = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveform_cross = (std::complex<double> *)malloc(sizeof(std::complex<double>)*length);
	fourier_waveform(frequencies, length, waveform_plus,waveform_cross, generation_method, parameters);
	for(int i = 0 ;i<length;i++)
	{
		waveform_plus_real[i] = real(waveform_plus[i]);
		waveform_plus_imag[i] = imag(waveform_plus[i]);
		waveform_cross_real[i] = real(waveform_cross[i]);
		waveform_cross_imag[i] = imag(waveform_cross[i]);
	}
	free(waveform_plus);
	free(waveform_cross);
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
	bool NSflag = parameters->NSflag;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag)
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
	fourier_waveform(frequencies, length, waveform, generation_method, parameters);
	for(int i = 0 ;i<length;i++)
	{
		waveform_real[i] = real(waveform[i]);
		waveform_imag[i] = imag(waveform[i]);
	}
	free(waveform);
}

/*!\brief Function to produce the amplitude of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 */
int fourier_amplitude(double *frequencies, /**< double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *amplitude, /**< output array for the amplitude*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters
			)
{
	int status=1;
	bool NSflag = parameters->NSflag;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag)
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
	//params.f_ref = parameters->f_ref;
	//params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		dCS_IMRPhenomD_log<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD")
	{
		dCS_IMRPhenomD<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		EdGB_IMRPhenomD_log<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		EdGB_IMRPhenomD<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}

	return status ;
}
/*!\brief Function to produce the phase of the (2,2) mode of an quasi-circular binary
 *
 * By using the structure parameter, the function is allowed to be more flexible in using different 
 * method of waveform generation - not all methods use the same parameters
 */
int fourier_phase(double *frequencies, /**<double array of frequencies for the waveform to be evaluated at*/
			int length,/**<integer length of all the arrays*/
			double *phase, /**<output array for the phase*/
			string generation_method,/**<String that corresponds to the generation method - MUST BE SPELLED EXACTLY*/
			gen_params *parameters
			)
{
	int status=1;
	bool NSflag = parameters->NSflag;

	/*Eventually, this will be where NS specific quantities are defined*/	
	if (NSflag)
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
	params.cosmology = parameters->cosmology;

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_phase(frequencies, length, phase, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
	}
	else if(generation_method == "dCS_IMRPhenomD_log")
	{
		bool local_spline = false;
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-1};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD_log<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		
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
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		dCS_IMRPhenomD<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD_log")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD_log<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "EdGB_IMRPhenomD")
	{
		params.betappe = parameters->betappe;
		params.Nmod = 1;
		int tempbppe[params.Nmod] = {-7};
		params.bppe = tempbppe;
		double temp[params.Nmod] ;
		for( int i = 0; i < params.Nmod; i++)
			temp[i] = params.betappe[i];
		EdGB_IMRPhenomD<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
		for( int i = 0; i < params.Nmod; i++)
			parameters->betappe[i] = temp[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		params.Nmod = parameters->Nmod;
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
	}

	return status ;
}



