#include "gwat/util.h"
#include "gwat/detector_util.h"
#include "gwat/io_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
#include <iostream>



int time_comparison_NADPN(int argc, char *argv[]);
int time_comparison_NADPN_MBH(int argc, char *argv[]);
int integration_interval_SM(int argc, char *argv[]);
int integration_interval_MBH(int argc, char *argv[]);
void RT_ERROR_MSG();

int main(int argc, char *argv[])
{
	std::cout<<"TESTING WAVEFORM UTIL CALCULATIONS"<<std::endl;
		
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		return time_comparison_NADPN(argc,argv);
	}
	else if(runtime_opt == 1){
		return time_comparison_NADPN_MBH(argc,argv);
	}
	else if(runtime_opt == 2){
		return integration_interval_SM(argc,argv);
	}
	else if(runtime_opt == 3){
		return integration_interval_MBH(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int integration_interval_SM(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = -.3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 1000;
	params.phic = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
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

	std::string method = "IMRPhenomPv2";

	double bounds[2];
	double fmin= 1e-5;
	double fmax= 1;
	int length = 1e5;
	double *freqs = new double[length];
	double *psd = new double[length];
	double *integrand = new double[length];
	std::complex<double> *response = new std::complex<double>[length];
	double *times = new double[length];
	double deltaf = (fmax-fmin)/length;
	for (int i =0 ; i<length ; i++){
		freqs[i]=fmin + i*deltaf;
	}
	populate_noise(freqs, "LISA_CONF",psd, length, 48);
	for(int i =0 ; i<length ; i++){
		psd[i]*=psd[i];
	}
	
	time_phase_corrected_autodiff(times, length, freqs, &params, method, false, (int *)NULL);
	fourier_detector_response(freqs,length, response, "LISA", method,&params, times);
	
	std::cout<<"OUTPUTTING"<<std::endl;
	for(int i = 0 ; i<length ; i++){
		integrand[i]=std::real( response[i] * std::conj(response[i]) ) /psd[i];
	}
	double **output = allocate_2D_array(length, 3);
	
	for(int i =0 ; i<length ; i++){
		output[i][0]=freqs[i];
		output[i][1]=integrand[i];
		output[i][2]=psd[i];
	}
	write_file("data/integration_bounds.csv",output, length, 3);
	deallocate_2D_array(output, length, 3);
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	integration_bounds(&params, method, "LISA","LISA_CONF",fmin,fmax,.1,0.01,bounds);
	std::cout<<bounds[0]<<" "<<bounds[1]<<std::endl;
	delete [] freqs;
	delete [] psd;
	delete [] integrand;
	delete [] response;
	delete [] times;
	return 0;
}
int integration_interval_MBH(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = -.3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 5000;
	params.phic = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
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
	params.mass1 = 9e6;
	params.mass2 = 4e6;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 0;
	params.equatorial_orientation = true;

	std::string method = "IMRPhenomPv2";

	double bounds[2];
	double fmin= 1e-5;
	double fmax= 1;
	int length = 1e5;
	double *freqs = new double[length];
	double *psd = new double[length];
	double *integrand = new double[length];
	std::complex<double> *response = new std::complex<double>[length];
	double *times = new double[length];
	double deltaf = (fmax-fmin)/length;
	for (int i =0 ; i<length ; i++){
		freqs[i]=fmin + i*deltaf;
	}
	populate_noise(freqs, "LISA_CONF",psd, length, 48);
	for(int i =0 ; i<length ; i++){
		psd[i]*=psd[i];
	}
	
	time_phase_corrected_autodiff(times, length, freqs, &params, method, false, (int *)NULL);
	fourier_detector_response(freqs,length, response, "LISA", method,&params, times);
	
	std::cout<<"OUTPUTTING"<<std::endl;
	for(int i = 0 ; i<length ; i++){
		integrand[i]=std::real( response[i] * std::conj(response[i]) ) /psd[i];
	}
	double **output = allocate_2D_array(length, 3);
	
	for(int i =0 ; i<length ; i++){
		output[i][0]=freqs[i];
		output[i][1]=integrand[i];
		output[i][2]=psd[i];
	}
	write_file("data/integration_bounds.csv",output, length, 3);
	deallocate_2D_array(output, length, 3);
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	integration_bounds(&params, method, "LISA","LISA_CONF",fmin,fmax,.1,0.01,bounds);
	std::cout<<bounds[0]<<" "<<bounds[1]<<std::endl;
	delete [] freqs;
	delete [] psd;
	delete [] integrand;
	delete [] response;
	delete [] times;
	return 0;
}
int time_comparison_NADPN_MBH(int argc, char *argv[])
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
	params.shift_time=true;
	params.shift_phase=true;
	
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
	params.mass1 = 9e6;
	params.mass2 = 4e6;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 0;
	params.equatorial_orientation = true;
	double fmin = 3e-5;
	double fmax = 1e-3;
	double T = T_year/2;

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
	write_file("data/times_N_MBH.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_AD[i];
	}
	write_file("data/times_AD_MBH.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_0PN[i];
	}
	write_file("data/times_0PN_MBH.csv",output,length,2);
	
	deallocate_2D_array(output,length,2);
	delete [] times_N;
	delete [] times_AD;
	delete [] times_0PN;
	
	delete [] frequency;
	return 0;
}
int time_comparison_NADPN(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .1;
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
	params.shift_time=true;
	params.shift_phase=true;
	
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
	double fmin = 3e-3;
	double fmax = 1e-2;
	double T = T_year/2;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}

	std::string method = "IMRPhenomD";

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
	std::cout<<"0 --- Time comparison N AD PN"<<std::endl;
	std::cout<<"1 --- Time comparison N AD PN for MBH (wiht merger)"<<std::endl;
	std::cout<<"2 --- Integration Interval Stellar Mass"<<std::endl;
	std::cout<<"3 --- Integration Interval MBH"<<std::endl;
}


