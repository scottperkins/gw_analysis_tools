#include <math.h>
#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ppE_utilities.h"
#include "gwat/error.h"
#include <iostream>
#include <iomanip>

int main(int argc, char *argv[]){

	/*
	The params come from test_LISA_fishers in test_fishers.cpp
	Except I have substituted the detector to LIGO
	*/


    //std::cout.precision(15);
	//Create injection structure
	gen_params params;	
	//params.mass1 = 1.9 *MSOL_SEC;
	params.mass1 = 2.0 *MSOL_SEC;
	params.mass2 = 1.4 *MSOL_SEC;
	params.spin1[2] = .001;
	params.spin2[2] = .002 ;
	params.Luminosity_Distance = 100;
	params.incl_angle = 3*M_PI/4;

	params.NSflag1 = true;
	params.NSflag2 =true;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2);//*MSOL_SEC;
	double eta = calculate_eta(params.mass1,params.mass2);//*MSOL_SEC;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	params.tidal_love=true;
	params.tidal_s=200;
	
	//params.tidal_love=false;
	//params.tidal1=200;
	//params.tidal2=100;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;

	

	//double fmin = 5;
	//double fmax = 2048;

	double fmin = 7.2;
	double fmax = 1820;
	
	//double **psd = new double*[Ndetect];
	//double fmin = .006508;
	//double fmax = .0067506;
	//double fmin = .01208;
	//double fmax = 1.00;
	//double T =(t_0PN(fmin,chirpmass)- t_0PN(fmax,chirpmass));
	double T = 16;
	//std::cout<<"TIME: "<<T/T_year<<std::endl;

	params.tc = 3.*T/4.;
	int length = 2000;
	double *frequency = new double[length];
	int Ndetect = 3;
	double **psd = new double*[Ndetect];
	
	bool AD = true;
	bool GL = false;
	double *weights = new double[length];
	if(AD && GL){
	//if(false){
		gauleg(log10(fmin), log10(fmax),frequency,weights,length);
		for(int i = 0 ; i<length; i++){
			frequency[i] = pow(10,frequency[i]);	
		}
	}
	else{
		double deltaF = (fmax-fmin)/length;	
		for(int i = 0 ; i<length; i++){
			frequency[i] = fmin + deltaF*i;
		}
	}

	std::cout<<"Freq populated"<<std::endl;

	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdLIGOMidHigh"}; //"AdVIRGOPlus1"};
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}
	
	std::cout<<"frequency[10]:"<<frequency[10]<<std::endl;

	int dim = 11;
	double* output = new double[dim];
	std::string method = "IMRPhenomD";
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	//std::string detector = "Hanford";

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
		}
	}

	double snr; 
	double total_snr = 0;

	//###############################################
	//Calculate Fishers
	//###############################################
	/*
	for(int i = 0 ;i < Ndetect; i++){
		        total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
			//debugger_print(__FILE__,__LINE__,total_snr);
			fisher_numerical(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, 2,NULL,NULL, psd[i]);
		
		for(int k = 0 ; k<dim; k++){
			std::cout<<i<<": "<<std::endl;
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
				std::cout<<output_AD[i][j]<<" ";
			}
			std::cout<<std::endl;
		}
	}
	std::cout<<"Total SNR: "<<sqrt(total_snr)<<std::endl;
*/
	//Get Waveforms
	
	std::complex<double> hpg[length];
	std::complex<double> hcg[length];
	std::complex<double> hpppE[length];
	std::complex<double> hcppE[length];
	waveform_polarizations<double> wp;
	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, "IMRPhenomD",&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, "ppE_IMRPhenomD_IMR",&params);


	//calculate_systematic_error(frequency, hcg, hcppE, length, method, detectors, detectors[0], output, dim, &params, 2, psd[1]);

	calculate_statistical_error(frequency, length, method, detectors, detectors[0], output, dim, &params, 2, psd);
	std::cout<<"SNR: "<<sqrt(output[0])<<std::endl;

	

	delete [] frequency;
	delete [] weights;
	delete [] psd;
	
	return 0;
}