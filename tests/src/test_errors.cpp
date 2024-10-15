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
#include <fstream>
#include <iomanip>
#include <vector>

int main(int argc, char *argv[]){

	/*
	The params come from test_LISA_fishers in test_fishers.cpp
	Except I have substituted the detector to LIGO
	*/


    //std::cout.precision(15);
	//Create injection structure
	gen_params params;	
	//params.mass1 = 1.9 *MSOL_SEC;
	params.mass1 = 1.44;
	params.mass2 = 1.29399;
	params.spin1[2] = .003;
	params.spin2[2] = -.002 ;
	params.Luminosity_Distance = 63;
	params.incl_angle = 2.532207345558998;

	params.NSflag1 = true;
	params.NSflag2 =true;

	params.phiRef = 2.;
	params.RA = 3.42;
	params.DEC = -.37;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2);//*MSOL_SEC;
	double eta = calculate_eta(params.mass1,params.mass2);//*MSOL_SEC;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	params.tidal_love=true;
	params.tidal_s=242;

	double beta = 3.;
	int b = -7.;
	std::cout<<"Beta: "<<beta<<std::endl;
	std::cout<<"b: "<<b<<std::endl;

	
	
	params.Nmod = 1;
	params.bppe = new double[1];
	params.bppe[0] = b;
	params.betappe = new double[1];
	params.betappe[0] = beta;

	params.alphappe = new double[1];
	params.appe = new double[1];
	params.alphappe[0] = 0.;
	params.appe[0] = -2.;

	
	
	//params.tidal_love=false;
	//params.tidal1=200;
	//params.tidal2=100;
	params.psi = 2.;
	double gps = 1187008882.4;
	params.gmst = gps_to_GMST_radian(gps);
	params.sky_average = false;

	

	//double fmin = 5;
	//double fmax = 2048;

	double fmin = 10;
	double fmax = 2048;
	
	//double **psd = new double*[Ndetect];
	//double fmin = .006508;
	//double fmax = .0067506;
	//double fmin = .01208;
	//double fmax = 1.00;
	//double T =(t_0PN(fmin,chirpmass)- t_0PN(fmax,chirpmass));
	double T = 16;
	//std::cout<<"TIME: "<<T/T_year<<std::endl;

	double Tsignal = 4;
	//double Tsignal = 128; 
	double deltaF = 1./Tsignal;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= Tsignal*3./4.;
	int length = (int)((fmax-fmin)/deltaF);
	params.tc=Tsignal-T_merger;
	double *frequency = new double[length];
	int Ndetect = 3;
	double **psd = new double*[Ndetect];
	
	bool AD = false;
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
		//double deltaF = (fmax-fmin)/length;	
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
		output[i] = 0;
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
	
	for(int i = 0 ;i < Ndetect; i++){
		if(AD){
			if(GL){
				total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "GAUSSLEG", weights, true), 2);
				fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
				//debugger_print(__FILE__,__LINE__, total_snr);

			}
			else{
				total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
				fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "SIMPSONS",weights,false, psd[i],NULL,NULL);
				//debugger_print(__FILE__,__LINE__, total_snr);
			}
		}
		else{
		        total_snr += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
			//debugger_print(__FILE__,__LINE__,total_snr);
			fisher_numerical(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, 2,NULL,NULL, psd[i]);
		}
		for(int k = 0 ; k<dim; k++){
			//std::cout<<i<<": "<<std::endl;
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
				//std::cout<<std::setprecision(5)<<output_AD[i][j]<<" ";
			}
			//std::cout<<std::endl;
		}
	}
	std::cout<<"Total SNR: "<<sqrt(total_snr)<<std::endl;

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

	double* output_sys = new double[dim];
	double* output_stat = new double[dim];

	

	//Generate a plot of systematic error vs SNR, based on a range of luminosity distances

	std::ofstream output_file("sys_err_data.csv");
	std::ofstream stat_output_file("stat_err_data.csv");
	std::ofstream waveform_output_file("waveform_data.csv");
	std::ofstream bgr_waveform_output_file("bgr_waveform_data.csv");


	
	int no_of_DL_steps = 20;
	int DL_step_size = 10.5;
	int DLeval = 63;

	double total_snr_temp = 0;

	for(int a = 0; a < no_of_DL_steps; a++){
	params.Luminosity_Distance = DLeval;
	total_snr_temp = 0;

	for(int i = 0 ; i<dim; i++){
		output_sys[i] = 0;
		output_stat[i] = 0;
	}
	for(int i = 0; i < Ndetect; i++){
	total_snr_temp += pow_int( calculate_snr(SN[i],detectors[i],method, &params, frequency, length, "SIMPSONS", weights, false), 2);
	}

	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, "IMRPhenomD",&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, "ppE_IMRPhenomD_IMR",&params);

	for(int i = 0; i < length; i++){
		waveform_output_file<<frequency[i]<<","<<hcg[i]<<std::endl;
		waveform_output_file<<frequency[i]<<","<<hcppE[i]<<std::endl;
	}

	calculate_systematic_error(frequency, hcg, hcppE, length, method, detectors, detectors[0], output_sys, dim, &params, 2, psd[1]);
	calculate_statistical_error(frequency, length, method, detectors, detectors[0], output_stat, dim, &params, 2, psd);

	output_file<<sqrt(total_snr_temp)<< ",";
	stat_output_file<<sqrt(total_snr_temp)<< ",";
	for(int i = 0; i < dim-1; i++){
		output_file<<output_sys[i]<<",";
		stat_output_file<<output_stat[i]<<",";
	}
	output_file<<output_sys[dim-1]<<std::endl;
	stat_output_file<<output_stat[dim-1]<<std::endl;

	DLeval+= DL_step_size;

	}
	output_file.close();
	std::vector<std::string> param_info = {"RA", "DEC", "psi", "phiRef", "tc", "iota_L", "ln DL", "ln chirpmass", "eta", "chi1", "chi2"};

	

	for(int i = 0; i < dim; i++){
	std::cout<<__LINE__<<" ? "<< param_info[i]<<" - Systematic Error "<<i<<": "<<output_sys[i]<<std::endl;
	}
	

	calculate_statistical_error(frequency, length, method, detectors, detectors[0], output_stat, dim, &params, 2, psd);
	//std::cout<<"SNR: "<<sqrt(output[0])<<std::endl;

	for(int i = 0; i < dim; i++){
	std::cout<<param_info[i]<<" - Statistical Error "<<i<<": "<<output_stat[i]<<std::endl;
	}

	for(int i = 0; i < dim; i++){
	std::cout<<param_info[i]<<" - Systematic Error / Statistical Error "<<i<<": "<<output_sys[i]/output_stat[i]<<std::endl;
	}

	

	delete [] frequency;
	delete [] weights;
	delete [] psd;
	
	return 0;
}