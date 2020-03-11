#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/io_util.h"
#include <iostream>

int time(int argc, char *argv[]);
void RT_ERROR_MSG();

int main(int argc, char *argv[])
{
	std::cout<<"TESTING WAVEFORM CALCULATIONS"<<std::endl;
		
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = stoi(argv[1]);	
	if(runtime_opt == 0){
		return time(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int time(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phic = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	
	//params.mass1 = 5e6;
	//params.mass2 = 5e6;
	//params.f_ref = 20;
	//params.psi = .1;
	//double fmin = 20;
	//double fmax = 2048;
	//params.incl_angle = .1;
	//params.tc = 10;
	//params.equatorial_orientation = false;
	//double T = 32;
	params.mass1 = 9e1;
	params.mass2 = 4e1;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 0;
	params.equatorial_orientation = true;
	double fmin = 3e-2;
	double fmax = 1e-1;
	double T = T_year/12;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}

	std::string method = "IMRPhenomPv2";

	double *times_N = new double[length];
	double *times_AD = new double[length];
	double *times_0PN = new double[length];
	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	time_phase_corrected(times_N, length, frequency, &params, method, false);
	double chirpm = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	for(int i = 0 ; i<length; i++){
		times_0PN[i]=t_0PN(frequency[i],chirpm);
	}
		
	double **output = allocate_2D_array(length,2);
	for(int i = 0 ; i<length; i++){
		output[i][0] = frequency[i];
		output[i][1] = times_N[i];
	}
	write_file("data/times_N.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_AD[i];
	}
	write_file("data/times_AD.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_0PN[i];
	}
	write_file("data/times_0PN.csv",output,length,2);
	
	deallocate_2D_array(output,length,2);
	delete [] times_N;
	delete [] times_AD;
	delete [] times_0PN;
	
	delete [] frequency;
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Test phase-time relationship"<<std::endl;
}
