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


    std::cout.precision(15);
	gen_params params;	
	params.mass1 = 36.9e5;
	params.mass2 = 35.93e5;
	params.spin1[2] = .8;
	params.spin2[2] = .8 ;
	params.chip = .07;
	params.phip = 1.0;
	params.Luminosity_Distance = 48000;
	//params.Luminosity_Distance = DL_from_Z(10,"PLANCK15");
	//params.Luminosity_Distance = DL_from_Z(0.0019,"PLANCK15");
	//params.incl_angle = acos(.75);
	params.incl_angle = 0;
	params.phi_l = M_PI/4;
	params.theta_l = M_PI/2;

	params.NSflag1 = false;
	params.NSflag2 =false;

	params.phiRef = .0;
	params.RA = .234;
	params.DEC = asin(-.45);
	params.f_ref = 1;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	std::cout<<"chirpmass: "<<chirpmass/MSOL_SEC<<std::endl;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	//params.equatorial_orientation=true;
	params.equatorial_orientation=false;
	
	params.gmst = 2.;
	//params.sky_average = false;
	params.sky_average = true;
	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new double[1];
	params.bppe[0]=-13;
	params.betappe[0]=0;
	

	double fmin = 1e-5;
	double fmax = 1e-3;
	
	//double **psd = new double*[Ndetect];
	//double fmin = .006508;
	//double fmax = .0067506;
	//double fmin = .01208;
	//double fmax = 1.00;
	double T =(t_0PN(fmin,chirpmass)- t_0PN(fmax,chirpmass));
	std::cout<<"TIME: "<<T/T_year<<std::endl;

	params.tc = 3.*T/4.;
	int length = 5000;
	double *frequency = new double[length];
	double *weights = new double[length];
	double *psd = new double[length];
	
	gauleg(log10(fmin),log10(fmax), frequency, weights, length);
	for(int j = 0 ; j<length; j++){
		frequency[j] = pow(10,frequency[j]);	
	}
	string SN = "AdLIGOMidHigh";
	populate_noise(frequency, SN,psd, length, 48);
	for(int j = 0 ; j<length; j++){
		psd[j]*=psd[j];	
	}
	
	
	int dim = 8;
	double* output = new double[dim];
	std::string method = "gIMRPhenomD";
	std::string detector = "Hanford";

	//Get Waveforms
	std::complex<double> hpg[length];
	std::complex<double> hcg[length];
	std::complex<double> hpppE[length];
	std::complex<double> hcppE[length];
	waveform_polarizations<double> wp;
	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, "gIMRPhenomD",&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, "ppE_IMRPhenomD_IMR",&params);


	//calculate_systematic_error(frequency, hpg, hpppE, length, method, detector,detector, output, dim, &params, 2, psd);
	//calculate_systematic_error(frequency, hcg, hcppE, length, method, detector, detector, output, dim, &params, 2, psd);
	//fisher_autodiff(frequency, length, method, detector,detector, output, dim, &params, "GAUSSLEG",weights,true, psd,NULL,NULL);
	std::cout<<"SNR: "<<sqrt(output[0])<<std::endl;

	

	delete [] frequency;
	delete [] weights;
	delete [] psd;
	
	return 0;
}