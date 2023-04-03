#include <iostream>
#include <gwat/waveform_util.h>
#include <gwat/waveform_generator.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/utilities.h>
#include <gwat/fisher.h>
#include <gwat/ortho_basis.h>

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
	injection.mass1 = 18.0;
	injection.mass2 = 14.0;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;

	//Spins
	injection.spin1[2] = .0;
	injection.spin2[2] = .0;
	injection.spin1[1] = .0;
	injection.spin2[1] = .0;
	injection.spin1[0] = .0;
	injection.spin2[0] = .0;

	//Extrinsic
	injection.Luminosity_Distance = 600;
	injection.psi = .2;
        injection.RA = 3.42;
	injection.DEC = -.37;
	injection.incl_angle = 2.532207345558998;
	double gps = 1187008882.4;
	injection.gmst = gps_to_GMST_radian(gps);

	//Trivial constants
	injection.phiRef = 2.;
	injection.f_ref = 20.;

	//Specific flags -- probably don't modify
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	//PhenomD injection
	std::string injection_method = "IMRPhenomD";
	
	//Detector properties
	int detect_number = 3;
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//Frequency array properties
	double fmin = 10;
	double fmax =2048;
	//Length of total signal -- leave at 4 seconds for now
	double Tsignal = 4;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= Tsignal*3./4.;
	//int length = (int)((fmax-fmin)/deltaF);
	int length = 1000;
	double deltaF = (fmax-fmin)/length;
	
	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];
	double **weights = new double*[detect_number];

	//Set merger time of injection
	injection.tc=Tsignal-T_merger;
	bool log10f = true;
	
	//Populate freq, psd, and allocate array for data
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		weights[i]= new double[length];
	
		gauleg(log10(fmin),log10(fmax), freq[i], weights[i], length);
		
		//Populate freq array
		for(int j =0 ; j<length; j++){
			freq[i][j] =pow(10,freq[i][j]);
		}
		//Populate PSD array
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
	}

	//Set data_length array 

	int dim = 11;
	double **fisher = new double*[dim];
	double **fisher_temp = new double*[dim];
	for(int i = 0 ; i<dim ; i++){
		fisher[i] = new double[dim];
		fisher_temp[i] = new double[dim];
		for(int j = 0 ; j<dim ; j++){
			fisher[i][j] = 0;
		}
	}
	for(int i = 0 ; i<detect_number ; i++){
		fisher_autodiff(freq[i], length, "IMRPhenomD",detectors[i], detectors[0], fisher_temp, dim, &injection, "GAUSSLEG", weights[i], log10f, psd[i], (int *)nullptr, (int *)nullptr);
		for(int i = 0 ; i<dim ; i++){
			for(int j = 0 ; j<dim ; j++){
				fisher[i][j] += fisher_temp[i][j];
			}
		}
	}
	
	bayesship::writeCSVFile("data/fisher_output.csv",fisher, dim, dim);

	double injectionOut[dim] = {injection.RA,injection.DEC, injection.psi, injection.incl_angle,injection.phiRef, injection.tc, log(injection.Luminosity_Distance), log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	bayesship::writeCSVFile("data/injection.csv",injectionOut,dim);
	
	
	for (int i = 0 ; i<dim ; i++){
		delete [] fisher[i];
		delete [] fisher_temp[i];
	}
	delete [] fisher;
	delete [] fisher_temp;
	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] weights[i];	
	}
	delete [] freq;
	delete [] psd;
	delete [] weights;

	return 0;
}
