#include <iostream>
#include <gwat/waveform_util.h>
#include <gwat/waveform_generator.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/utilities.h>
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
	injection.mass2 = 4.0;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;

	//Spins
	injection.spin1[2] = .8;
	injection.spin2[2] = -.3;
	injection.spin1[1] = .1;
	injection.spin2[1] = -.7;
	injection.spin1[0] = .1;
	injection.spin2[0] = .3;

	//Extrinsic
	injection.Luminosity_Distance = 600;
	injection.psi = .8;
        injection.RA = 3.42;
	injection.DEC = -.37;
	//injection.incl_angle = 2.532207345558998;
	injection.incl_angle =2.;
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
	std::string injection_method = "IMRPhenomPv2";
	
	//Detector properties
	int detect_number = 3;
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	//Frequency array properties
	double fmin = 10;
	double fmax =2048;
	//Length of total signal -- leave at 4 seconds for now
	double Tsignal = 4;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	double T_merger= Tsignal*3./4.;
	//int length = (int)((fmax-fmin)/deltaF);
	int length = 5000;
	double deltaF = (fmax-fmin)/length;
	
	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	std::complex<double> *hplus = new std::complex<double>[length];
	std::complex<double> *hcross = new std::complex<double>[length];
	double *phaseplus = new double[length];
	double *phasecross = new double[length];
	double **freq = new double*[detect_number];
	int data_lengths[detect_number]; 

	//Set merger time of injection
	
	//Populate freq and allocate array for data
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i] = length;
		data[i]= new std::complex<double>[length];
		freq[i] = new double[length];
		for(int j=0 ; j<length; j++){
			freq[i][j] = fmin + j * deltaF;
		}
	}
	
	//injection.tc=Tsignal - T_merger;
	injection.tc=T_merger;
	
	
	create_coherent_GW_detection(detectors, detect_number, freq, data_lengths, true, &injection, injection_method, data);

	waveform_polarizations<double> *wp = new waveform_polarizations<double>;
	wp->hplus = hplus;
	wp->hcross = hcross;
	fourier_waveform(freq[0], length, wp, injection_method, &injection);
	fourier_phase(freq[0], length, phaseplus, phasecross,  injection_method, &injection);
	delete wp;

	double **outputResponse = new double*[detect_number * 2 + 1];
	for(int i = 0 ; i<detect_number * 2 + 1; i++){
		outputResponse[i] = new double[length];
	}
	double **outputPhase = new double*[2 + 1];
	for(int i = 0 ; i<3; i++){
		outputPhase[i] = new double[length];
	}
	double **outputPol = new double*[2*2 + 1];
	for(int i = 0 ; i<5; i++){
		outputPol[i] = new double[length];

	}
	for ( int i = 0 ; i< length; i++){
		outputResponse[0][i] = freq[0][i];	
		outputPhase[0][i] = freq[0][i];	
		outputPol[0][i] = freq[0][i];	
	}

	for ( int i = 0 ; i< detect_number; i++){
		for(int j = 0 ; j<length; j++){
			outputResponse[2*i+1][j] = std::real(data[i][j]);
			outputResponse[2*i+2][j] = std::imag(data[i][j]);
		}
	}
	for(int j = 0 ; j<length; j++){
		outputPhase[1][j] = phaseplus[j];
		outputPhase[2][j] = phasecross[j];
		outputPol[1][j] = std::real(hplus[j]);
		outputPol[2][j] = std::imag(hplus[j]);
		outputPol[3][j] = std::real(hcross[j]);
		outputPol[4][j] = std::imag(hcross[j]);
	}
	
	bayesship::writeCSVFile("data/Responses.csv",outputResponse, 2*detect_number+1, length);
	bayesship::writeCSVFile("data/Phases.csv",outputPhase, 3, length);
	bayesship::writeCSVFile("data/Polarizations.csv",outputPol, 5, length);

	for(int i = 0 ; i<detect_number * 2 + 1; i++){
		delete [] outputResponse[i] ;
	}
	delete [] outputResponse ;
	for(int i = 0 ; i<3; i++){
		delete [] outputPhase[i] ;
	}
	delete [] outputPhase ;
	for(int i = 0 ; i<5; i++){
		delete [] outputPol[i] ;

	}
	delete [] outputPol ;
	
//	bayesship::writeCSVFile("data/fisher_output.csv",fisher, dim, dim);
//
//	//double injectionOut[dim] = {injection.RA,injection.DEC, injection.psi, injection.incl_angle,injection.phiRef, injection.tc, log(injection.Luminosity_Distance), log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
//	bayesship::writeCSVFile("data/injection.csv",injectionOut,dim);
	
	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] data[i];	
	}
	delete [] freq;
	delete [] hplus;
	delete [] hcross;
	delete [] phaseplus;
	delete [] phasecross;
	delete [] data;	

	return 0;
}
