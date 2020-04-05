#include "gwat/util.h"
#include "gwat/detector_util.h"
#include "gwat/io_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ortho_basis.h"
#include <iostream>



int time_comparison_NADPN(int argc, char *argv[]);
int time_comparison_NADPN_MBH(int argc, char *argv[]);
int integration_interval_SM(int argc, char *argv[]);
int integration_interval_MBH(int argc, char *argv[]);
int observation_interval_SM(int argc, char *argv[]);
int observation_interval_MBH(int argc, char *argv[]);
int threshold_times_MBH(int argc, char *argv[]);
int threshold_times_SM(int argc, char *argv[]);
int tbm_testing(int argc, char *argv[]);
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
	else if(runtime_opt == 4){
		return tbm_testing(argc,argv);
	}
	else if(runtime_opt == 5){
		return threshold_times_SM(argc,argv);
	}
	else if(runtime_opt == 6){
		return threshold_times_MBH(argc,argv);
	}
	else if(runtime_opt == 7){
		return observation_interval_SM(argc,argv);
	}
	else if(runtime_opt == 8){
		return observation_interval_MBH(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
/*
 * Runs the regular, AD derivative to get t(f) --
 *
 */
int tbm_testing(int argc, char *argv[])
{
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = -.3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	params.sky_average=true;
	
	//params.mass1 = 36;
	//params.mass2 = 29;
	params.mass1 = 36e5;
	params.mass2 = 29e4;
	params.theta_l = 1;
	params.phi_l = 1;
	params.phiRef= 10;
	params.equatorial_orientation = false;
	params.psi = 0;
	params.incl_angle = 0;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;

	std::string method = "IMRPhenomD";
	//std::string method = "ppE_IMRPhenomPv2_IMR";
	//params.Nmod = 1;
	//params.betappe = new double[params.Nmod];
	//params.betappe[0]=100;
	//params.bppe = new int[params.Nmod];
	//params.bppe[0]=-1;


	double fpeak, frd, fdamp;
	postmerger_params(&params, method, &fpeak,&fdamp, &frd);
	std::cout<<fpeak<<std::endl;

	double fmin= 1e-6;
	double fmax= 5e-3;
	//double fmin= 1e-2;
	//double fmax= 1;
	int length = 1e3;
	double *freqs = new double[length];
	double *times = new double[length];
	double deltaf = (fmax-fmin)/length;
	for (int i =0 ; i<length ; i++){
		freqs[i]=fmin + i*deltaf;
	}
	
	params.tc = 0.;
	//params.tc = 3./(deltaf*4);

	time_phase_corrected_autodiff(times, length, freqs, &params, method, false, (int *)NULL);
	
	double **output = allocate_2D_array(length, 2);
	
	for(int i =0 ; i<length ; i++){
		output[i][0]=freqs[i];
		output[i][1]=times[i];
	}
	write_file("data/Tbm_times_base.csv",output, length, 2);

	params.tc = 3./(deltaf*4);
	
	double **output2=allocate_2D_array(length,2);
	bool autodiff=false;
	double tol = 1e-10;
	autodiff=true;
	int max_iteration = 100;
	std::cout.precision(15);
	clock_t start,end;
	double ave=0;
	bool relative_time =true;
	for(int i =0 ; i<length; i++){
		output2[i][1]=output[i][1];
		start = clock();	
		Tbm_to_freq(&params, method,output2[i][1], &(output2[i][0]) ,tol,autodiff,max_iteration,relative_time);
		end = clock();
		ave += (double)(end - start)/CLOCKS_PER_SEC;
	}
	std::cout<<"Average time for AD: "<<ave/length<<std::endl;
	ave=0;
	write_file("data/Tbm_times_inv_AD.csv",output2, length, 2);

	autodiff=false;
	for(int i =0 ; i<length; i++){
		start = clock();	
		Tbm_to_freq(&params, method,output2[i][1], &(output2[i][0]) ,tol,autodiff,max_iteration,relative_time);
		end = clock();
		ave += (double)(end - start)/CLOCKS_PER_SEC;
	}
	std::cout<<"Average time for N: "<<ave/length<<std::endl;
	write_file("data/Tbm_times_inv_N.csv",output2, length, 2);

	for(int i =0 ; i<length; i++){
		start = clock();	
		output2[i][0] = f_0PN(output2[i][1],chirpmass);
		end = clock();
		ave += (double)(end - start)/CLOCKS_PER_SEC;
	}
	std::cout<<"Average time for PN: "<<ave/length<<std::endl;
	write_file("data/Tbm_times_inv_PN.csv",output2, length, 2);

	deallocate_2D_array(output2,length,2);
	deallocate_2D_array(output, length, 2);
	delete [] freqs;
	delete [] times;
	//delete [] params.betappe;
	//delete [] params.bppe;
	return 0;
}
int threshold_times_SM(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.spin1[2] = .0;
	params.spin2[2] = -.0;
	params.chip = .5;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 0.;
	params.DEC = -0.1;
	params.f_ref = 1e0;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 36;
	params.mass2 = 29;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 10*T_year;
	params.equatorial_orientation = false;
	params.incl_angle = 0;
	params.sky_average=true;

	std::string method = "IMRPhenomD";

	double bounds[2];
	double fmin= 1e-5;
	double fmax= 1;
	double Tobs = 4*T_year;
	double Twait = 20*T_year;
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	bool autodiff = false;
	int np = 1000;
	double tol = 1e-10;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);

	clock_t start = clock();
	int status =threshold_times_gsl(&params, method, Tobs,Twait,fmin,fmax,"LISA_SADC_CONF",8.,bounds,tol,w,np);
	std::cout<<"TIME: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	std::cout<<bounds[0]/T_year<<" "<<bounds[1]/T_year<<std::endl;
	gsl_integration_workspace_free(w);


	int snr_pts = 500;
	double freqs[snr_pts];
	double weights[snr_pts];
	int iterations = 1000;
	double delta_t = pow(200./.1,1./iterations);
	double **snr_out = new double*[iterations];
	double fbounds[2];
	autodiff = true;
	int max_iterations = 50;
	double fpeak, frd, fdamp;
	postmerger_params(&params, method, &fpeak,&fdamp, &frd);
	bool relative_time = true;
	for(int i = 0 ; i<iterations; i++){
		snr_out[i]=new double[2];
		snr_out[i][0]= .1*T_year*pow_int(delta_t,i);	
		Tbm_to_freq(&params, method,snr_out[i][0], &(fbounds[0]) ,tol,autodiff,max_iterations,relative_time);
		if(snr_out[i][0]-Tobs >0){
			Tbm_to_freq(&params, method,(snr_out[i][0]-Tobs), &(fbounds[1]) ,tol,autodiff,max_iterations,relative_time);
		}
		else{
			if(fpeak<fmax){
				fbounds[1] = 2.*fpeak;
			}
			else{
				fbounds[1] = fmax;
			}
			
		}
		gauleg(log10(fbounds[0]), log10(fbounds[1]), freqs, weights, snr_pts);
		for(int i  = 0 ; i<snr_pts; i++){
			freqs[i]  = pow(10.,freqs[i]);
		}
		snr_out[i][1] = calculate_snr("LISA_SADC_CONF","LISA",method, &params, freqs, snr_pts, "GAUSSLEG",weights, true);

	}
	write_file("data/threshold_output_SM.csv", snr_out, iterations, 2);
	

	for(int i = 0 ; i<iterations; i++){
		delete [] snr_out[i];
	}
	delete [] snr_out;
	
	//delete [] freqs;
	//delete [] psd;
	//delete [] integrand;
	//delete [] response;
	//delete [] times;
	return 0;
}
int threshold_times_MBH(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.chip = .5;
	params.phip = 0.1;
	//params.spin1[2] = 0.7022*cos( 2.572);
	//params.spin2[2] = 0.7873*cos(1.6073);
	params.spin1[2] = 0.407*cos( 2.9);
	params.spin2[2] = 0.998*cos(2.17);
	params.Luminosity_Distance = DL_from_Z(3.704,"Planck15");
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=true;
	params.sky_average=true;
	
	//params.mass1 =1.0142e+04;
	//params.mass2 = 0.3773e+03;
	params.mass1 = 5.1277e7;
	params.mass2 = 7.117e5;
	params.theta_l = 1;
	params.phi_l = 2;
	params.equatorial_orientation = false;
	params.incl_angle = 0;
	params.psi = 0;

	std::string method = "IMRPhenomD";

	double bounds[2];
	double fmin= 1e-5;
	double fmax= 1;
	double Tobs = 10*T_year;
	double Twait = 200000000*T_year;
	//params.tc = 3.*Tobs/4.;
	params.tc = 0.*Tobs/4.;
	std::cout<<"tc in years: "<<params.tc/T_year<<std::endl;
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	bool autodiff = false;
	int np = 1000;
	double tol = 1e-10;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);
	clock_t start = clock();
	int status =threshold_times_gsl(&params, method, Tobs,Twait,fmin,fmax,"LISA_SADC_CONF",8.,bounds,tol,w,np);
	std::cout<<"STATUS: "<<status<<std::endl;
	std::cout<<"TIME: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;
	std::cout<<bounds[0]/T_year<<" "<<bounds[1]/T_year<<std::endl;
	gsl_integration_workspace_free(w);


	int snr_pts = 500;
	double freqs[snr_pts];
	double weights[snr_pts];
	int iterations = 100;
	double delta_t = pow(20./.1,1./iterations);
	double **snr_out = new double*[iterations];
	double fbounds[2];
	autodiff = true;
	int max_iterations = 50;
	double fpeak, frd, fdamp;
	postmerger_params(&params, method, &fpeak,&fdamp, &frd);
	std::cout<<"Fpeak: "<<fpeak<<std::endl;
	bool relative_time=true;
	double chirpmass=calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	tol = 1e-8;
	for(int i = 0 ; i<iterations; i++){
		snr_out[i]=new double[2];
		snr_out[i][0]= .1*T_year*pow_int(delta_t,i);	
		int status = Tbm_to_freq(&params, method,snr_out[i][0], &(fbounds[0]) ,tol,autodiff,max_iterations,relative_time);
		if(snr_out[i][0]-Tobs >0){
			status = Tbm_to_freq(&params, method,(snr_out[i][0]-Tobs), &(fbounds[1]) ,tol,autodiff,max_iterations,relative_time);
		}
		else{
			if(2*fpeak<fmax){
				fbounds[1] = fpeak*2.;
			}
			else{
				fbounds[1] = fmax;
			}
			
		}
		//double time[2];
		//time_phase_corrected_autodiff(time,2,fbounds, &params,"IMRPhenomD",false);
		//std::cout<<fbounds[0]<<" "<<fbounds[1]<<" "<<(-time[1]+time[0])/T_year<<std::endl;
		
		gauleg(log10(fbounds[0]), log10(fbounds[1]), freqs, weights, snr_pts);
		for(int i  = 0 ; i<snr_pts; i++){
			freqs[i]  = pow(10.,freqs[i]);
		}
		snr_out[i][1] = calculate_snr("LISA_SADC_CONF","LISA",method, &params, freqs, snr_pts, "GAUSSLEG",weights, true);

	}
	write_file("data/threshold_output_MBH.csv", snr_out, iterations, 2);
	

	for(int i = 0 ; i<iterations; i++){
		delete [] snr_out[i];
	}
	delete [] snr_out;
	
	return 0;
}
int observation_interval_SM(int argc, char *argv[])
{
	std::cout<<" TESTING stellar mass observation interval"<<std::endl;
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = -.3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 1000;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 9e1;
	params.mass2 = 4e1;
	params.theta_l = 1;
	params.phi_l = 2;
	//params.tc = T_year;
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
	double **output = allocate_2D_array(length, 4);
	
	for(int i =0 ; i<length ; i++){
		output[i][0]=freqs[i];
		output[i][1]=integrand[i];
		output[i][2]=psd[i];
		output[i][3]=times[i];
	}
	write_file("data/observation_bounds.csv",output, length, 4);
	deallocate_2D_array(output, length, 4);
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	bool autodiff = true;
	observation_bounds(2.,4*T_year, "LISA","LISA_CONF",method,&params,bounds,autodiff);
	std::cout<<bounds[0]<<" "<<bounds[1]<<std::endl;
	delete [] freqs;
	delete [] psd;
	delete [] integrand;
	delete [] response;
	delete [] times;
	return 0;
}
int observation_interval_MBH(int argc, char *argv[])
{
	std::cout<<" TESTING massive mass observation interval"<<std::endl;
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = -.3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 6600;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 4.6e7;
	params.mass2 = 4.5e7;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = T_year;
	//params.tc = 0;
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
	double **output = allocate_2D_array(length, 4);
	
	for(int i =0 ; i<length ; i++){
		output[i][0]=freqs[i];
		output[i][1]=integrand[i];
		output[i][2]=psd[i];
		output[i][3]=times[i];
	}
	write_file("data/observation_bounds.csv",output, length, 4);
	deallocate_2D_array(output, length, 4);
	
	std::cout<<"BOUNDS CALC"<<std::endl;
	bool autodiff = true;
	observation_bounds(2.,4*T_year, "LISA","LISA_CONF",method,&params,bounds,autodiff);
	std::cout<<bounds[0]<<" "<<bounds[1]<<std::endl;
	delete [] freqs;
	delete [] psd;
	delete [] integrand;
	delete [] response;
	delete [] times;
	return 0;
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
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 4e1;
	params.mass2 = 2e1;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = T_year;
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
	bool autodiff = true;
	integration_bounds(&params, method, "LISA","LISA_CONF",fmin,fmax,.1,0.01,bounds,autodiff);
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
	params.spin1[2] = .5;
	params.spin2[2] = -.3;
	params.chip = .5;
	params.phip = 0.1;
	params.Luminosity_Distance = 5000;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 9e6;
	params.mass2 = 4e5;
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
	bool autodiff = false;
	integration_bounds(&params, method, "LISA","LISA_CONF",fmin,fmax,.1,0.01,bounds,autodiff);
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
	params.Luminosity_Distance = 4608.28;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 =5.3357e+08*(1+0.44288);
	params.mass2 =932820*(1+0.44288);
	params.theta_l = 1;
	params.phi_l = 2;
	//params.tc = T_year*3./4.;
	params.tc = 0;
	params.equatorial_orientation = true;
	double fmin = 3e-7;
	double fmax = 5e-4;
	double T = 4*T_year;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}

	std::string method = "IMRPhenomD";

	double fpeak, frd, fdamp;
	postmerger_params(&params, method, &fpeak,&fdamp, &frd);
	std::cout<<fpeak<<std::endl;

	double *times_N = new double[length];
	double *times_AD = new double[length];
	double *times_0PN = new double[length];
	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	time_phase_corrected(times_N, length, frequency, &params, method, false,1);
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
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 3.6e1;
	params.mass2 = 2.9e1;
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
	std::cout<<"1 --- Time comparison N AD PN for MBH (with merger)"<<std::endl;
	std::cout<<"2 --- Integration Interval Stellar Mass"<<std::endl;
	std::cout<<"3 --- Integration Interval MBH"<<std::endl;
	std::cout<<"4 --- Tbm (time before merger) testing"<<std::endl;
	std::cout<<"5 --- threhold times for LISA detection Stellar Mass"<<std::endl;
	std::cout<<"6 --- threhold times for LISA detection MBH"<<std::endl;
	std::cout<<"7 --- observation bounds for LISA detection Stellar Mass"<<std::endl;
	std::cout<<"8 --- observation bounds for LISA detection MBH"<<std::endl;
}


