#include <iostream>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "ppE_IMRPhenomD.h"
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
	params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++)
			waveform_cross[i] = std::complex<double>(0,1) * waveform_plus[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);
		for (int i =0 ; i < length; i++)
			waveform_cross[i] = std::complex<double>(0,1) * waveform_plus[i];
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		status = ppemodeld.construct_waveform(frequencies, length, waveform_plus, &params);	
		for (int i =0 ; i < length; i++)
			waveform_cross[i] = std::complex<double>(0,1) * waveform_plus[i];
	}
	else if(generation_method == "IMRPhenomPv2")
	{
		IMRPhenomPv2<double> modeld;
		//Initialize Pv2 specific params	
		modeld.PhenomPv2_Param_Transform(&params);
		//Calculate Waveform
		status = modeld.construct_waveform(frequencies, length, waveform_plus, waveform_cross, &params);
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
	params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);

	params.phi = parameters->phi;
	params.theta = parameters->theta;
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
		status = ppemodeld.construct_waveform(frequencies, length, waveform, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
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
	params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		status = ppemodeld.construct_amplitude(frequencies, length, amplitude, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
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
	params = params.populate_source_parameters(mass1, mass2, Luminosity_Distance, spin1, spin2, phi_c,t_c);

	if(generation_method == "IMRPhenomD")
	{
		IMRPhenomD<double> modeld;
		status = modeld.construct_phase(frequencies, length, phase, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_Inspiral")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		ppE_IMRPhenomD_Inspiral<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
	}
	else if(generation_method == "ppE_IMRPhenomD_IMR")
	{
		params.betappe = parameters->betappe;
		params.bppe = parameters->bppe;
		ppE_IMRPhenomD_IMR<double> ppemodeld;
		status = ppemodeld.construct_phase(frequencies, length, phase, &params);	
	}

	return status ;
}

//int fourier_waveform_polarizations(double *frequencies,
//				int length,
//				std::complex<double> hcross,
//				std::complex<double> hplus,
//				string generation_method,
//				gen_params *parameters
//				)
//{
//	int status = 1;
//	std::complex<double> *waveform = 
//			(std::complex<double> *)malloc(sizeof(std::complex<double>*length);	
//	if(generation_method == "IMRPhenomD")
//	{
//		std::complex<double> *waveform = 
//			(std::complex<double> *)malloc(sizeof(std::complex<double>*length);	
//		IMRPhenomD<double> modeld;
//		status = modeld.construct_phase(frequencies, length, waveform, &params);	
//		for (int i = 0; i<length;i++)
//		{
//			hcross[i]=waveform[i];
//			hplus[i]=waveform[i]*std::complex<double>(0.,1.);
//		}
//		free(waveform);
//	}
//	return status
//}


