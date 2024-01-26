#include <iostream>
#include <gwat/waveform_util.h>
#include <gwat/waveform_generator.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_gw_extended.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/standardPriorLibrary.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/utilities.h>
#include <bayesship/bayesshipSampler.h>

// Example injection for an initially non-spinning
// black hole binary.

int main(int argc, char *argv[])
{
	//######################################
	//######################################
	//######################################
	//Injection Parameters
	//######################################
	//######################################
	//######################################

	gen_params injection;

	//Masses in M_SOL, plus useful other parameters
	injection.mass1 = 2.0e6;
	injection.mass2 = 1.0e6; 
	injection.spin1[2] = .9; //spin_1z
	injection.spin2[2] = -.8; //spin_2z
	injection.spin1[1] = 0.0; //spin_1y
	injection.spin2[1] = 0.0;//spin_2y
	injection.spin1[0] = 0.0;//spin_1x
	injection.spin2[0] = 0.0;//spin_2x
	injection.Luminosity_Distance = 1.0e8; //std::atof(argv[3]);// Mpc
	injection.psi = 2.3;
	injection.RA = 1.7;
	injection.DEC = 1.05;
	injection.incl_angle = M_PI/3.0;
	double gps = 1126259642.413;
	injection.gmst = gps_to_GMST_radian(gps);
	injection.phiRef = 0;
	injection.f_ref = 1.0; 
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	//PhenomD injection
	std::string injection_method = "IMRPhenomD";
	
	//Detector properties
	int detect_number = 3;
	std::string detectors[1] = {"LISA_TDIAET_full"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//Frequency array properties
	double fmin = 1e-4;
	double fmax = 1.;
	//Length of total signal -- leave at 4 seconds for now
	double Tsignal = 3600*24*365;
	//double deltaF = 1./Tsignal;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= 0.0;
	//int length = (int)((fmax-fmin)/deltaF);
	int length = 40000;
	double deltaF = (std::log(fmax) - std::log(fmin))/(length-1);
	
	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];

	//Set merger time of injection
	injection.tc=0.;
	
	
	
	//Populate freq, psd, and allocate array for data
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		//Populate freq array
		
		freq[i][0] = fmin;
		for(int j =1 ; j<length; j++){
			freq[i][j] = std::exp(std::log(fmin) + j*deltaF);		
		}
		//for(int j =0 ; j<length; j++){
		//	freq[i][j] = fmin + j*deltaF;
		//}
		
	}
	
	
	
	
	
	int out_rows_wp = 5;//Number of cols in output array
	double **data_out_wp = new double*[out_rows_wp]; //Declare output array
	for(int i = 0 ; i<out_rows_wp; i++){
		data_out_wp[i] = new double[length];
	}
	waveform_polarizations<double> *wp = new waveform_polarizations<double>();//Structure to store polarizations
	wp->allocate_memory(length);//Allocate memory for polarizations
	fourier_waveform(freq[0], length, wp, injection_method, &injection);//Generate polarizations
	for(int i = 0 ; i<length; i++){
		data_out_wp[0][i] = freq[0][i];
		data_out_wp[1][i] = std::real(wp->hplus[i]); //wp->hplus[i]=data_in[column1][i]+std::complex<double>(0,1)*data_in[column2][i]
		data_out_wp[2][i] = std::imag(wp->hplus[i]);
		data_out_wp[3][i] = std::real(wp->hcross[i]);
		data_out_wp[4][i] = std::imag(wp->hcross[i]);
	}
	wp->deallocate_memory();//Deallocate memory for wp structure
	delete wp;
	//Write out freq, real(hplus), imag(hplus)
	write_file("data/waveform_debugging.csv",data_out_wp, 5, length);
	
	
	for(int i = 0 ; i<out_rows_wp; i++){
		delete [] data_out_wp[i];
	}
	delete [] data_out_wp;
	
	
	

	//Set data_length array 
	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}

	//Create coherent gravitational wave detector between all detectors
	create_coherent_GW_detection(detectors,detect_number, freq, data_lengths,  true, &injection, injection_method, data);
	
	
	int data_out_rows = detect_number*2+1; //Number of cols in output array
	double **data_out_new = new double*[data_out_rows]; //Declare output array
	for(int i = 0 ; i<data_out_rows; i++){
		data_out_new[i] = new double[length];
	}
	//Write out data
	for(int i = 0 ; i<length; i++){
		data_out_new[0][i] = freq[0][i] ; //Write first freq set, since they're identical
		int ct = 0 ;
		for(int j = 1 ; j<data_out_rows; j+=2){ //change detector_number to out_rows
			data_out_new[j][i] = std::real(data[ct][i]);
			data_out_new[j+1][i] = std::imag(data[ct][i]);
			ct+=1;
		}
	}
	//Write out freq, real(response_1), imag(reasponse_1) ...
	write_file("data/responses_debugging.csv",data_out_new, data_out_rows, length);
	for(int i = 0 ; i<data_out_rows; i++){
		delete [] data_out_new[i];
	}
	delete [] data_out_new;
	
	
	
	
	
	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] injection.bppe;
	delete [] injection.betappe;
	delete [] freq;
	delete [] psd;
	delete [] data;

	return 0;
}
