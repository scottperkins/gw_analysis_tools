#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include <string>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "mcmc_routines.h"
#include "noise_util.h"
#include "util.h"
#include "waveform_util.h"
#include <adolc/adouble.h>
#include "fisher.h"
#include "ppE_IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "waveform_generator_C.h"
#include "mcmc_sampler.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_rng.h"
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/taping.h"
#include "limits"
#include "GWATConfig.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


using namespace std;

bool show_progress = true;
void test1();
void test2();
void test3();
void test4();
void test5();
void test6();
void test7();
void test8();
void test9();
void test10();
void test11();
void test12();
void test13();
void test14();
void test15();
void test16();
void test17();
void test18();
void test19();
void test20();
double test_ll(double *pos, int dim);
double test_lp(double *pos, int dim);
double test_lp_nts(double *pos, int dim, int chain_id);
double test_lp_GW(double *pos, int dim, int chain_id);
double test_lp_GW_Pv2(double *pos, int dim, int chain_id);
double test_lp_GW_7dim(double *pos, int dim, int chain_id);
double test_lp_GW_DFull(double *pos, int dim, int chain_id);
double test_lp_GW_dCS(double *pos, int dim, int chain_id);
double test_lp_GW_dCS_log(double *pos, int dim, int chain_id);
double test_lp_GW_ppE(double *pos, int dim, int chain_id);
void test_fisher(double *pos, int dim, double **fisher);
double log_student_t (double *x,int dim);
double log_neil_proj3 (double *x,int dim);
double log_neil_proj3_nts (double *x,int dim, int chain_id);
double log_neil_proj32 (double *c,int dim);
void fisher_neil_proj3 (double *x,int dim, double **fish);
adouble dist(adouble *pos, int dimension);

const gsl_rng_type* Y;
gsl_rng * g;
static std::complex<double> *waveformout11=NULL;
static double *freq=NULL;
static double *psd=NULL;

int main(){

	test20();	
	return 0;
}
void test20()
{
	std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(50000,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(5,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, -5};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0,4* log(50000/(3e8))};
	//double initial_pos[dimension]={.9, 2, 1.,log(300),log(10), .24,- .0,-.0,-0};
	double initial_pos[dimension]={.9, 2, -1.,log(300),log(30), .24,- .0,-.0,-0};
	double *seeding_var = NULL;
	std::cout<<initial_pos[8]<<std::endl;
	int n_steps = 1000000;
	int chain_N=10 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.3;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 1;
	int *bppe = new int[Nmod];
	bppe[0] = -1;
	int numThreads = 8;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "ppE_IMRPhenomD_Inspiral";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD";
	
	
	std::string autocorrfile = "";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_ppE.csv";
	std::string chainfile = "testing/data/mcmc_output_ppE.csv";
	std::string statfilename = "testing/data/mcmc_statistics_ppE.txt";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_ppE,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[0][j][0];
		output_transform[j][1]=output[0][j][1];
		output_transform[j][2]=output[0][j][2];
		output_transform[j][3]=std::exp(output[0][j][3]);
		output_transform[j][4]=std::exp(output[0][j][4]);
		output_transform[j][5]=output[0][j][5];
		output_transform[j][6]=output[0][j][6];
		output_transform[j][7]=output[0][j][7];
		//output_transform[j][8]=std::exp(output[0][j][8]);
		output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	//chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	chainfile = "testing/data/mcmc_output_ppE_hot.csv";
	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[chain_N-1][j][0];
		output_transform[j][1]=output[chain_N-1][j][1];
		output_transform[j][2]=output[chain_N-1][j][2];
		output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
		output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
		output_transform[j][5]=output[chain_N-1][j][5];
		output_transform[j][6]=output[chain_N-1][j][6];
		output_transform[j][7]=output[chain_N-1][j][7];
		//output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
		output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	delete [] bppe;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test19()
{
	gen_params params;
	double RA=1, DEC=.2;
	double gps_time = 1185389807.3;//TESTING -- gw170729
	params.mass1 = 80;
	params.mass2 = 40;
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin1[2] = .1;
	params.spin2[2] = .1;
	params.Nmod = 1;
	params.betappe= new double[params.Nmod];
	params.betappe[0] = -1;
	params.sky_average = false;
	params.NSflag = false;
	params.incl_angle = M_PI/3.;
	params.theta = 0;
	params.phi = 0;
	params.Luminosity_Distance = 400;
	params.phic = 0;
	params.tc = 0;
	std::string gen_method = "MCMC_dCS_IMRPhenomD_log_Full";
	std::string detector = "Hanford";
	
	int length = 4000;
	double freqs[length];
	double fhigh = 1024;
	double flow = 10;
	double fstep = (fhigh-flow)/length;
	for(int i =0; i<length; i++){
		freqs[i] = flow + i*fstep;
	}
	double noise[length];
	populate_noise(freqs,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	int dim = 9;
	double **fishout = allocate_2D_array(dim, dim);
	std::complex<double> response[length];
	std::complex<double> responsedcs[length];
	double tc=2;
	params.tc = tc;
	double phic;
	double phicmax, tcmax;
	fftw_outline plan;
	initiate_likelihood_function(&plan, length);

	//fisher(freqs, length, gen_method, detector,fishout,dim, &params,NULL, NULL, noise);
	fourier_detector_response(freqs, length, response, detector, "IMRPhenomD",&params);
	fourier_detector_response(freqs, length, responsedcs, detector, "dCS_IMRPhenomD_log",&params);


	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	double llmax = maximized_coal_Log_Likelihood(response, noise, freqs, length, &params,
			"Hanford", "dCS_IMRPhenomD_log",&plan, &tcmax,&phicmax);
	std::cout<<tcmax<<" "<<phicmax<<std::endl;
	double temptheta = params.theta;
	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	double ll = Log_Likelihood(response, noise, freqs, length, &params,
			"Livingston", "dCS_IMRPhenomD_log", &plan);
	std::cout.precision(15);
	std::cout<<llmax<<std::endl;
	std::cout<<ll<<std::endl;
	std::cout<<ll+llmax<<std::endl;
	for(int i =0; i<length/4;i++){
	//for(int i =0; i<dim;i++){
		for(int j =0; j<dim; j++){
			//std::cout<<i<<" "<<j<<" "<<fishout[i][j]<<std::endl;
			//std::cout<<i<<" "<<j<<" "<<fishout[i][j]<<std::endl;
		}
		//std::cout<<responsedcs[i]-response[i]<<std::endl;
	}
	deallocate_2D_array(fishout, dim , dim);
	delete [] params.betappe;
	deactivate_likelihood_function(&plan);
}
void test18()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	std::string psd_file = "testing/data/GWTC1_GW151226_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	int num_detectors = 2, psd_length = 8032, length;
	//int num_detectors = 3, psd_length = 4016, length;
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1185389792-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 9;
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(50000,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, log(pow(MPC_SEC,4)*pow(5,4))};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0, -5};
	//double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0,4* log(50000/(3e8))};
	double initial_pos[dimension]={.9, 2, 1.,log(300),log(10), .24,- .0,-.0,-40};
	//double initial_pos[dimension]={.9, 2, 1.,log(3000),log(50), .24,- .0,-.0,5e-18};
	double *seeding_var = NULL;
	std::cout<<initial_pos[8]<<std::endl;
	int n_steps = 15000;
	int chain_N=6 ;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 1;
	int *bppe = NULL;
	int numThreads = 10;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "dCS_IMRPhenomD_log";
	//std::string generation_method = "EdGB_IMRPhenomD_log";
	//std::string generation_method = "dCS_IMRPhenomD";
	
	
	//std::string autocorrfile = "";
	std::string autocorrfile = "testing/data/auto_corr_mcmc_dCS.csv";
	std::string chainfile = "testing/data/mcmc_output_dCS.csv";
	std::string statfilename = "testing/data/mcmc_statistics_dCS.txt";
	//std::string autocorrfile = "testing/data/auto_corr_mcmc_EdGB.csv";
	//std::string chainfile = "testing/data/mcmc_output_EdGB.csv";
	//std::string statfilename = "testing/data/mcmc_statistics_EdGB.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_dCS_log,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile);	

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[0][j][0];
		output_transform[j][1]=output[0][j][1];
		output_transform[j][2]=output[0][j][2];
		output_transform[j][3]=std::exp(output[0][j][3]);
		output_transform[j][4]=std::exp(output[0][j][4]);
		output_transform[j][5]=output[0][j][5];
		output_transform[j][6]=output[0][j][6];
		output_transform[j][7]=output[0][j][7];
		output_transform[j][8]=std::exp(output[0][j][8]);
		//output_transform[j][8]=output[0][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	//chainfile = "testing/data/mcmc_output_EdGB_hot.csv";
	chainfile = "testing/data/mcmc_output_dCS_hot.csv";
	for(int j = 0; j<n_steps;j++){
		output_transform[j][0]=output[chain_N-1][j][0];
		output_transform[j][1]=output[chain_N-1][j][1];
		output_transform[j][2]=output[chain_N-1][j][2];
		output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
		output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
		output_transform[j][5]=output[chain_N-1][j][5];
		output_transform[j][6]=output[chain_N-1][j][6];
		output_transform[j][7]=output[chain_N-1][j][7];
		output_transform[j][8]=std::exp(output[chain_N-1][j][8]);
		//output_transform[j][8]=output[chain_N-1][j][8];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test17()
{
	std::cout<<Z_from_DL(1000.,"PLANCK13")<<std::endl;
	std::cout<<Z_from_DL(10000.,"PLANCK13")<<std::endl;
	std::cout<<Z_from_DL(100.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(10.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(1.,"PLANCK13")<<std::endl;
	std::cout<<DL_from_Z(.1,"PLANCK13")<<std::endl;
	//clock_t start7,end7;
	//int threads =10;
	//start7 = clock();
	//gsl_spline **z_d_spline = (gsl_spline **)malloc(sizeof(gsl_spline*)*threads);	
	//gsl_interp_accel **z_d_accel= (gsl_interp_accel **)malloc(sizeof(gsl_interp_accel *)*threads);	
	//for(int i =0; i<threads; i++){
	//	initiate_LumD_Z_interp(&z_d_accel[i], &z_d_spline[i]);
	//}
	//int iterations = 100;
	//omp_set_num_threads(threads);
	//#pragma omp parallel for 
	//for (int j =0; j<iterations; j++){
	//	gsl_spline *z_d_spline;
	//	gsl_interp_accel *z_d_accel;
	//	initiate_LumD_Z_interp(&z_d_accel, &z_d_spline);
	//	std::cout<<j<<std::endl;
	//	std::cout<<Z_from_DL_interp(j+10, z_d_accel, z_d_spline)<<std::endl;;
	//	free_LumD_Z_interp(&z_d_accel, &z_d_spline);
	//}

	////for (int i =0; i<threads;i++)
	////	free_LumD_Z_interp(&z_d_accel[i], &z_d_spline[i]);
	//end7 = clock();

	//cout<<"TIMING interp: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
}

void test16()
{
	//std::string psd_file = "testing/data/GWTC1_GW150914_PSDs.dat.txt";
	//std::string data_file = "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//std::string psd_file = "testing/data/GWTC1_GW170729_PSDs.dat.txt";
	std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW170608/GWTC1_GW170608_PSDs.dat.txt";
	//std::string psd_file = "/Users/sperkins/Downloads/LOSC_data/GW150914/GWTC1_GW150914_PSDs.dat.txt";
	//int rows = 8032;
	//int cols = 3;
	int datalength = 131075;
	//double **psd = allocate_2D_array(rows, cols);
	//read_LOSC_PSD_file(psd_file, psd, rows, cols);
	//double data_start_time, duration, fs;
	//int num_detectors = 2, psd_length = 8032, length;
	int num_detectors = 2, psd_length = 16064, length;
	//int num_detectors = 3, psd_length = 4016, length;
	//double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1185389807.3;//TESTING -- gw170729
	double gps_time = 1180922494.5;//TESTING -- gw170608
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	//detectors[2] = "Virgo";
	std::string *detector_files = new std::string[num_detectors];
	//detector_files[0] =  "testing/data/H-H1_GWOSC_4KHZ_R1-1135136335-32.txt";
	//detector_files[1] =  "testing/data/L-L1_GWOSC_4KHZ_R1-1135136335-32.txt";
	detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/H-H1_GWOSC_4KHZ_R1-1180922479-32.txt";
	detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW170608/L-L1_GWOSC_4KHZ_R1-1180922479-32.txt";
	//detector_files[0] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[1] =  "/Users/sperkins/Downloads/LOSC_data/GW150914/L-L1_GWOSC_4KHZ_R1-1126259447-32.txt";
	//detector_files[2] =  "testing/data/V-V1_GWOSC_4KHZ_R1-1185389792-32.txt";
 	//double trigger_time= 1135136350.6;
 	double trigger_time = gps_time;
	double **psd = allocate_2D_array(num_detectors,psd_length);
	double **freqs = allocate_2D_array(num_detectors,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*num_detectors);
	for(int i =0; i<num_detectors; i++)
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);

	allocate_LOSC_data(detector_files, psd_file, num_detectors, psd_length, datalength, trigger_time, data, psd, freqs);


	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =psd_length;
	data_length[1] =psd_length;
	//data_length[2] =psd_length;

	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={.3, 2., -0.2,log(400),log(30), .24,- .0,-.0};
	double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0};
	//double initial_pos[dimension]={-.9, 2, -1.2,log(410),log(30), .24,-.4,.3};
	//double initial_pos[dimension]={-.0, 0, -0,log(500),log(50), .2,-.0,.0};
	//double initial_pos[dimension]={-.99, 2, -1.2,log(410),log(30.78), .24,-.4,.3};
	double *seeding_var = NULL;
	int n_steps = 3000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 10;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 5;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_DFull.csv";
	std::string statfilename = "testing/data/mcmc_statistics_DFull.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,freqs, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_DFull_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	delete [] detectors;
	free(data_length);
	//free_LOSC_data(data, psd,freqs, num_detectors, length);
	deallocate_2D_array(psd,num_detectors, psd_length);
	deallocate_2D_array(freqs,num_detectors, psd_length);
	for(int i =0; i<num_detectors; i++)
		free(data[i]);
	free(data);
	delete [] detector_files;
	//deallocate_2D_array(psd, rows, cols);
	//free(data);
	
	//int length = 2000;
	//double x[length];
	//double y[length];
	//double xlim = 10.;
	//double xstart=0;
	//double xstep = (xlim-xstart)/length;
	//for (int i =0; i<length; i++){
	//	x[i]= i*xstep;
	//	y[i]= x[i]*x[i];
	//}
	//double sum = simpsons_sum(xstep, length, y);
	//std::cout<<sum<<std::endl;

}
void test15()
{
	int length = 5000;
	int num_detectors =3;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	//double RA = 5.;
	//double DEC = 1.;
	double RA = 1.5;
	double DEC = .2;
	double chirpm = 20.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = .5;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	params.phic = .0;
	double tc = 2;
	params.Luminosity_Distance = 410.;
	params.NSflag = false;
	params.incl_angle = M_PI/3.;
	params.sky_average=false;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =1000;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	
	fourier_detector_response(freq,length, waveformout,"Hanford","IMRPhenomD", &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	//double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	double snr2 = pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Hanford", "IMRPhenomD",&params),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	std::cout<<"Time at Livingston of injection: "<<params.tc<<std::endl;
	fourier_detector_response(freq,length, waveformout,"Livingston","IMRPhenomD", &params);

	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Livingston", "IMRPhenomD",&params),2);
	celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	std::cout<<"Time at Virgo of injection: "<<params.tc<<std::endl;
	fourier_detector_response(freq,length, waveformout,"Virgo","IMRPhenomD", &params);

	for(int j = 0; j<data_length[0]; j++){
		frequencies[2][j] = freq[j];	
		psd[2][j] = (noise[j]);	
		data[2][j] = waveformout[j];	
	}
	//#########################################################
	//snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	snr2 += pow(data_snr_maximized_extrinsic(freq,length,waveformout, noise,"Virgo", "IMRPhenomD",&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	//snr = data_snr_maximized_extrinsic(freq,length, waveformout,"Hanford_O1_fitted","IMRPhenomD",params );
	//std::cout<<"SNR of injection calculated as data: "<<snr<<std::endl;

	
	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={.3, 2., -0.2,log(400),log(40), .24,- .0,-.0};
	double initial_pos[dimension]={.0, 1, 0.,log(300),log(10), .2,- .0,-.0};
	//double initial_pos[dimension]={-.9, 2, -1.2,log(410),log(30), .24,-.4,.3};
	//double initial_pos[dimension]={-.0, 0, -0,log(500),log(50), .2,-.0,.0};
	//double initial_pos[dimension]={-.99, 2, -1.2,log(410),log(30.78), .24,-.4,.3};
	double *seeding_var = NULL;
	int n_steps = 20000;
	int chain_N= 5;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 3;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 3;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_injection.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_injection.csv";
	std::string statfilename = "testing/data/mcmc_statistics_injection.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool,show_progress,
			num_detectors, 
			data, psd,frequencies, data_length,gps_time, detectors,Nmod, bppe,
			generation_method,statfilename,"",autocorrfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_injection_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test14()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	
	int raw_length=  4096;
	//int raw_length=  8032;
	int cutoff = 50;
	int length = raw_length-cutoff;

	int num_detectors =2;
	//int num_detectors =1;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	//std::string filebase = "testing/data/gw170608_";
	std::string filebase = "testing/data/gw151226_";
	//std::string filebase = "testing/data/gw150914_";
	//std::string filebase = "testing/data/gw_150914_";
	//double gps_time = 1126259462;//TESTING -- gw150914
	double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	//read_file(filebase+"data_H_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H_new.csv",temp_psd);
	//read_file(filebase+"freq_H_new.csv",temp_freq);
	//read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L.csv",temp_psd);
	//read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd.csv",temp_psd);
	//read_file(filebase+"freq.csv",temp_freq);
	
	double delta_f = temp_freq[length/2]-temp_freq[length/2-1];
	//double delta_f = 1;
	//double delta_f = 1./(temp_freq[length/2]-temp_freq[length/2-1]);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] = 1./delta_f * (temp_psd[j+cutoff]);	
		data[0][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data_L_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L_new.csv",temp_psd);
	//read_file(filebase+"freq_L_new.csv",temp_freq);
	//read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H.csv",temp_psd);
	//read_file(filebase+"freq_H.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] = 1./delta_f *(temp_psd[j+cutoff]);	
		data[1][j] = 1./delta_f *std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);

	//######### SNR #####################	
	gen_params params;
	//150914
	//double RA = 1.69;
	//double DEC = -1.21;
	//double chirpm = 30.78;
	//double eta =.24;
	//params.mass1 = 38.56;
	//params.mass2 = 33.18;
	//chirpm = calculate_chirpmass(params.mass1,params.mass2);
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = -.0;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .0;
	//params.phic = .0;
	//double tc = 2;
	//params.Luminosity_Distance = 431.;
	////params.Luminosity_Distance = 800.;
	//params.NSflag = false;
	//params.incl_angle = M_PI;
	//params.sky_average=false;
	//params.phi = 0;
	//params.theta = 0;
	//params.tc = tc; 
	
	//151226
	//double RA = 3.34;
	//double DEC = -.5;
	//double chirpm = 30.78;
	//params.mass1 = 15.26;
	//params.mass2 = 8.29;
	//chirpm = calculate_chirpmass(params.mass1,params.mass2);
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = .5;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .0;
	//params.phic = .0;
	//double tc = 2;
	//params.Luminosity_Distance = 425.;
	//params.NSflag = false;
	//params.incl_angle = M_PI/3.;
	//params.sky_average=false;


	//170806
	double RA = 2.1;
	double DEC = .5;
	double chirpm = .78;
	params.mass1 = 11.9;
	params.mass2 = 8.09;
	chirpm = calculate_chirpmass(params.mass1,params.mass2);
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = .3;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	params.phic = .0;
	double tc = 2;
	params.Luminosity_Distance = 425.;
	params.NSflag = false;
	params.incl_angle = 0;
	params.sky_average=false;

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, detectors[0], &params.phi, &params.theta);
	double snr2 = pow(data_snr_maximized_extrinsic(frequencies[0],length,data[0], psd[0],detectors[0], "IMRPhenomD",&params),2);
	celestial_horizon_transform(RA, DEC, gps_time, detectors[1], &params.phi, &params.theta);
	snr2 += pow(data_snr_maximized_extrinsic(frequencies[1],length,data[1], psd[1],detectors[1], "IMRPhenomD",&params),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;

	//#########################################################
	//mcmc options
	int dimension = 8;
	//double initial_pos[dimension]={-.9, 1.,-1.,log(400),log(30), .24,- .0,-.0};
	double initial_pos[dimension]={cos(params.incl_angle)*.9, RA,DEC,
		log(params.Luminosity_Distance),log(chirpm), .24,- .0,-.0};
	//double initial_pos[dimension]={-.0, 0.,-0.,log(400),log(10), .24,- .0,-.0};
	//double initial_pos[dimension]={log(8), .24,- .0,-.0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int n_steps = 50000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 10;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.4;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 8;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_DFull_Neil.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_DFull_Neil.csv";
	std::string statfilename = "testing/data/mcmc_statistics_DFull_Neil.txt";
	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_DFull,numThreads, pool, show_progress, 
			num_detectors, 
			data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[0][j][0];
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=std::exp(output[0][j][4]);
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
			output_transform[j][7]=output[0][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_DFull_hot_Neil.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0]=output[chain_N-1][j][0];
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=std::exp(output[chain_N-1][j][3]);
			output_transform[j][4]=std::exp(output[chain_N-1][j][4]);
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
			output_transform[j][7]=output[chain_N-1][j][7];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test13()
{
	double RV, RH, RL;
	RV = radius_at_lat(V_LAT, V_elevation);
	RH = radius_at_lat(H_LAT, H_elevation);
	RL = radius_at_lat(L_LAT, L_elevation);
	double dtoa = DTOA(80,80,"Hanford","Virgo");
	
	cout.precision(15);
	std::cout<<"Distance to center of the earth for Virgo: "<<RV<<std::endl;
	std::cout<<"Distance to center of the earth for Livingston: "<<RL<<std::endl;
	std::cout<<"Distance to center of the earth for Hanford: "<<RH<<std::endl;
	std::cout<<"DTOA: "<<dtoa<<std::endl;
}
void test12()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	
	int raw_length=  4096;
	int cutoff = 50;
	int length = raw_length-cutoff;

	int num_detectors =2;
	//int num_detectors =1;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	//std::string filebase = "testing/data/gw170608_";
	//std::string filebase = "testing/data/gw151226_";
	std::string filebase = "testing/data/gw150914_";
	//std::string filebase = "testing/data/gw_150914_";
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_h.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_h.csv",temp_psd);
	read_file(filebase+"freq_h.csv",temp_freq);
	//read_file(filebase+"data.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd.csv",temp_psd);
	//read_file(filebase+"freq.csv",temp_freq);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] = (temp_psd[j+cutoff]);	
		data[0][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_l.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_l.csv",temp_psd);
	read_file(filebase+"freq_l.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] = (temp_psd[j+cutoff]);	
		data[1][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	
	//#########################################################
	//mcmc options
	int dimension = 7;
	double initial_pos[dimension]={0,log(30), .24, .0,.0,0,0};
	double *seeding_var = NULL;
	int n_steps = 100;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	int swp_freq = 5;
	double chain_temps[chain_N];
	chain_temps[0]=1.;
	double c = 1.1;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
	
	int Nmod = 0;
	int *bppe = NULL;
	int numThreads = 20;
	bool pool = false;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomPv2";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_Pv2.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output_Pv2.csv";
	std::string statfilename = "testing/data/mcmc_statistics_Pv2.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_Pv2,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0] = output[0][j][0];
			output_transform[j][1]=std::exp(output[0][j][1]);
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
			output_transform[j][4]=output[0][j][4];
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_Pv2_hot.csv";
	for(int j = 0; j<n_steps;j++){
			output_transform[j][0] = output[chain_N-1][j][0];
			output_transform[j][1]=std::exp(output[chain_N-1][j][1]);
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=output[chain_N-1][j][3];
			output_transform[j][4]=output[chain_N-1][j][4];
			output_transform[j][5]=output[chain_N-1][j][5];
			output_transform[j][6]=output[chain_N-1][j][6];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test11()
{
	int length = 4000;
	int loops = 1000;
	
	//std::complex<double> **out = (std::complex<double> **)malloc(sizeof(std::complex<double>)*loops);

	//for (int i =0 ; i< loops; i++)
	//{
	//	out[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*length);
	//}
	double fhigh =200;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	freq = (double *)malloc(sizeof(double) * length);
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	//Synthetic data
	double gps_time = 1135136350.6;//TESTING -- gw151226
	gen_params params_data;
	double chirpm = 50.78;
	double eta =.24;
	double Ra = 1;
	double Dec = -1;
	params_data.mass1 = calculate_mass1(chirpm,eta);
	params_data.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//complex<double> waveformout[length];
	params_data.spin1[0] = 0;
	params_data.spin1[1] = 0;
	params_data.spin1[2] = -.4;
	params_data.spin2[0] = 0;
	params_data.spin2[1] = 0;
	params_data.spin2[2] = .3;
	params_data.phic = .0;
	params_data.tc = 5;
	params_data.Luminosity_Distance = 400.;
	params_data.NSflag = false;
	celestial_horizon_transform(Ra, Dec, gps_time, "Hanford", &params_data.phi, &params_data.theta);	
	//params_data.phi = 0;
	//params_data.theta = 0;
	params_data.incl_angle = M_PI/3.;
	params_data.sky_average=false;
	waveformout11 = (std::complex<double>*) malloc(sizeof(std::complex<double>)*length);
	std::complex<double> *waveformout11L = (std::complex<double>*) malloc(sizeof(std::complex<double>)*length);
	
	//fourier_detector_response(freq, length, waveformout11,method,&params_data);
	fourier_detector_response(freq, length, waveformout11,"Hanford",method,&params_data);
	double temptheta=params_data.theta;
	celestial_horizon_transform(Ra, Dec, gps_time, "Livingston", &params_data.phi, &params_data.theta);	

	params_data.tc = params_data.tc + DTOA(temptheta, params_data.theta,"Hanford", "Livingston");
	fourier_detector_response(freq, length, waveformout11L,"Livingston",method,&params_data);

	//double psd[length];
	psd = (double *) malloc(sizeof(double)*length);
	populate_noise(freq,"Hanford_O1_fitted", psd,length);
	for (int i =0; i<length;i++){
		psd[i] = psd[i]*psd[i];
	}
	//fftw_outline *plans = (fftw_outline *)malloc(sizeof(fftw_outline)*loops);	
	//for(int i =0 ; i<loops; i ++)
	//{
	//	initiate_likelihood_function(&plans[i] , length);
	//}
	fftw_outline plan;
	initiate_likelihood_function(&plan, length);
	double ll_true =maximized_Log_Likelihood(waveformout11, psd, freq,
			length, &params_data, "Livingston", "IMRPhenomD", &plan);
	std::cout<<"TRUE ll: "<<ll_true<<std::endl;	
	omp_set_num_threads(10);

	double *DLs = (double *)malloc(sizeof(double)*loops);
	double *lls1 = (double *)malloc(sizeof(double)*loops);
	double *lls2 = (double *)malloc(sizeof(double)*loops);
	double *lls3 = (double *)malloc(sizeof(double)*loops);
	double DLmin = 50;
	double DLmax = 1000;
	double DLstep = (DLmax - DLmin)/loops;
	for (int i =0; i<loops; i++)
	{
		DLs[i] = DLstep * i + DLmin;
	}
	//#pragma omp parallel 
	{
		//#pragma omp for
		for (int i =0; i<loops; i ++)
		{
			gen_params params;
			double chirpm = 50.78;
			double eta =.24;
			params.mass1 = calculate_mass1(chirpm,eta);
			params.mass2 = calculate_mass2(chirpm,eta);
			string method= "IMRPhenomD";
			params.spin1[0] = 0;
			params.spin1[1] = 0;
			params.spin1[2] = -.4;
			params.spin2[0] = 0;
			params.spin2[1] = 0;
			params.spin2[2] = .3;
			params.phic = .0;
			params.tc = 0;
			params.Luminosity_Distance = DLs[i];
			params.NSflag = false;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			params.incl_angle = 0;
			params.sky_average=false;

			double tc,phic;
			lls1[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			//lls1[i]=0;
			//tc=5;
			//phic = 0;

			//std::cout<<"LL 1: "<<lls1[i]<<std::endl;
			double phi1 = params.phi;
			double theta1 = params.theta;
			//std::cout<<"HANFORD ANGLES: "<<params.theta<<" "<<params.phi<<std::endl;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			//std::cout<<"LIV ANGLES: "<<params.theta<<" "<<params.phi<<std::endl;
			//std::cout<<"HANFORD TIME: "<<tc<<std::endl;
			params.tc = tc + DTOA(theta1,params.theta, "Hanford", "Livingston");
			//std::cout<<"LIV TIME: "<<params.tc<<std::endl;
			//std::cout<<"LIV TIME REAL: "<<params_data.tc<<std::endl;
			params.phic = phic;
			//params.tc = params_data.tc;
			lls1[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);
			//params.tc = 0;
			//params.phic =0;
			//lls1[i] +=maximized_coal_Log_Likelihood(waveformout11L, psd, freq,
			//		length, &params, "Livingston", "IMRPhenomD", &plan, &tc, &phic);
			//std::cout<<"LL 2: "<<lls1[i]<<std::endl;

			params.incl_angle = M_PI/3.;
			params.tc = 0;
			params.phic = 0;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			lls2[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			//std::cout<<"HANFORD TIME: "<<tc<<std::endl;
			phi1 = params.phi;
			theta1 = params.theta;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			params.tc = tc + DTOA(theta1,params.theta, "Hanford", "Livingston");
			//std::cout<<"LIV TIME: "<<params.tc<<std::endl;
			params.phic = phic;
			lls2[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);

			params.incl_angle = -M_PI/4.;
			params.tc = 0;
			params.phic = 0;
			celestial_horizon_transform(Ra, Dec, gps_time, "Hanford",&params.phi,&params.theta);	
			lls3[i] =maximized_coal_Log_Likelihood(waveformout11, psd, freq,
					length, &params, "Hanford", "IMRPhenomD", &plan, &tc, &phic);
			phi1 = params.phi;
			theta1 = params.theta;
			celestial_horizon_transform(Ra, Dec, gps_time, "Livingston",&params.phi,&params.theta);	
			
			params.tc =tc+ DTOA(theta1,params.theta, "Hanford", "Livingston");
			params.phic = phic;
			lls3[i] +=Log_Likelihood(waveformout11L, psd, freq,
					length, &params, "Livingston", "IMRPhenomD", &plan);
			//if(ll != ll_true)
			//{
			//	std::cout<<ll<<std::endl;
			//}
		}
	}
	double **output = allocate_2D_array(loops, 4);
	for (int i=0;i<loops; i++){
		output[i][0] = DLs[i];
		output[i][1] = lls1[i];
		output[i][2] = lls2[i];
		output[i][3] = lls3[i];
	}
	std::string filename = "testing/data/dl_ll.csv";
	write_file(filename, output, loops, 4);
	//for(int i =1;i < loops; i++)
	//{
	//	for(int j =0; j <length; j++){
	//		if(out[i-1][j]!=out[i][j]){
	//			std::cout<<"YIKES"<<std::endl;
	//		}
	//	}
	//}

	//for(int i =0;i < loops; i++)
	//{
	//	free(out[i]);
	//}
	//free(out);
	//for(int i =0 ; i<loops; i ++)
	//{
	//	deactivate_likelihood_function(&plans[i]);
	//}
	//free(plans);
	deactivate_likelihood_function(&plan);
	deallocate_2D_array(output, loops,2);
	free(freq);
	free(DLs);
	free(lls1);
	free(lls2);
	free(lls3);
	free(psd);
	free(waveformout11);
	free(waveformout11L);
}
void test10()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	int cutoff = 40;
	int raw_length=  4096;
	int length = raw_length-cutoff;
	//int length = raw_length-cutoff;
	int num_detectors =2;

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	std::string filebase = "testing/data/gw150914_";
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	//bool check = true;
	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] = (temp_psd[j+cutoff]);	
		data[0][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] = (temp_psd[j+cutoff]);	
		data[1][j] = std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	//#########################################################
	//MCMC options
	int dimension = 7;
	double initial_pos[dimension]={log(400),2,2,log(30), .24, 0,0};
	//double initial_pos[dimension]={log(200*MPC_SEC),log(20*MSOL_SEC), .15, 0,0};
	double *seeding_var = NULL;
	int N_steps = 50000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 20;
	bool pool = false;
	//#########################################################
	//GW options
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc_7dim.csv";
	std::string chainfile = "testing/data/mcmc_output_7dim.csv";
	std::string statfilename = "testing/data/mcmc_statistics_7dim.txt";

	MCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW_7dim,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length,gps_time, detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=std::exp(output[0][j][3]);
			output_transform[j][4]=output[0][j][4];
			output_transform[j][5]=output[0][j][5];
			output_transform[j][6]=output[0][j][6];
	}
	write_file(chainfile, output_transform, N_steps, dimension);
	//ofstream mcmc_out;
	//mcmc_out.open("testing/data/mcmc_output.csv");
	//mcmc_out.precision(15);
	////for(int i = 0;i<chain_N;i++){
	//for(int j = 0; j<N_steps;j++){
	//	//for(int k = 0; k<dimension; k++){
	//		mcmc_out<<std::exp(output[0][j][0])/MPC_SEC<<" , "<<std::exp(output[0][j][1])/MSOL_SEC<<" , "<<output[0][j][2]<<" , "<<output[0][j][3]<<" , "<<output[0][j][4]<<endl;
	//	//}
	//}
	////}
	//mcmc_out.close();

	deallocate_3D_array(output, chain_N, N_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< N_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	//free(detectors);
	free(data_length);
	
}
void test9()
{
	//int raw_length = 28673;
	//int cutoff = 600;
	//int high_cut = 11000;
	//int length = high_cut-cutoff;
	
	int raw_length=  4096;
	//int raw_length=  8032;
	int cutoff = 50;
	int length = raw_length-cutoff;

	int num_detectors =2;
	//int num_detectors =1;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";

	double **temp_data = allocate_2D_array(raw_length,2);
	double *temp_psd = (double *)malloc(sizeof(double)*raw_length);
	double *temp_freq = (double *)malloc(sizeof(double)*raw_length);
	//std::string filebase = "testing/data/gw170608_";
	//std::string filebase = "testing/data/gw151226_";
	std::string filebase = "testing/data/gw150914_";
	//std::string filebase = "testing/data/gw_150914_";
	double gps_time = 1126259462;//TESTING -- gw150914
	//double gps_time = 1135136350.6;//TESTING -- gw151226
	//double gps_time = 1180922494.5; //TESTING -- gw170608

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	read_file(filebase+"data_H.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_H.csv",temp_psd);
	read_file(filebase+"freq_H.csv",temp_freq);
	//read_file(filebase+"data_H_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_H_new.csv",temp_psd);
	//read_file(filebase+"freq_H_new.csv",temp_freq);
	//read_file(filebase+"data.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd.csv",temp_psd);
	//read_file(filebase+"freq.csv",temp_freq);
	//double delta_f = temp_freq[length/2]-temp_freq[length/2-1];
	double delta_f = 1;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = temp_freq[j+cutoff];	
		psd[0][j] =  (temp_psd[j+cutoff]);	
		data[0][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}
	read_file(filebase+"data_L.csv",temp_data, raw_length,2);
	read_file(filebase+"psd_L.csv",temp_psd);
	read_file(filebase+"freq_L.csv",temp_freq);
	//read_file(filebase+"data_L_new.csv",temp_data, raw_length,2);
	//read_file(filebase+"psd_L_new.csv",temp_psd);
	//read_file(filebase+"freq_L_new.csv",temp_freq);
	for(int j = 0; j<data_length[1]; j++){
		frequencies[1][j] = temp_freq[j+cutoff];	
		psd[1][j] =  (temp_psd[j+cutoff]);	
		data[1][j] = 1./delta_f * std::complex<double>(temp_data[j+cutoff][0],temp_data[j+cutoff][1]);	
	}

	deallocate_2D_array(temp_data,raw_length,2);
	free(temp_psd);
	free(temp_freq);
	
	//gen_params params;
	//imrphenomd<double> modeld;
	//params.mass1 = 63.187;
	//params.mass2 = 17.331;
	//string method= "imrphenomd";
	//params.spin1[0] = 0;
	//params.spin1[1] = 0;
	//params.spin1[2] = -.16;
	//params.spin2[0] = 0;
	//params.spin2[1] = 0;
	//params.spin2[2] = .01388;
	//params.phic = .0;
	//params.tc = 0;
	//params.luminosity_distance = 100.;
	//params.nsflag = false;
	//params.phi = 0;
	//params.theta = 0;
	//params.incl_angle = 0;
	//params.sky_average=false;
	//fftw_outline plan;
	//initiate_likelihood_function(&plan, length);
	//double ll =maximized_log_likelihood(data[0], psd[0], frequencies[0],
	//		length, &params, "hanford", "imrphenomd", &plan);
	//double snr = data_snr_maximized_extrinsic(frequencies[0],length, data[0],"Hanford_O1_fitted","imrphenomd",params );
	//std::cout<<"ll "<<ll<<std::endl;
	//std::cout<<"snr "<<snr<<std::endl;
	//deactivate_likelihood_function(&plan);
	//exit(1);
	//#########################################################
	//mcmc options
	int dimension = 4;
	//double initial_pos[dimension]={log(30), .24,- .0,-.0};
	double initial_pos[dimension]={log(8), .12,- .0,-.0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int n_steps = 40000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, n_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 5;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	int Nmod = 0;
	int *bppe = NULL;
	
	int numThreads = 20;
	bool pool = true;
	//#########################################################
	//gw options
	std::string generation_method = "IMRPhenomD";
	
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	//std::string autocorrfile = "";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";

	MCMC_MH_GW(output, dimension, n_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW,numThreads, pool,show_progress, num_detectors, data, psd, 
			frequencies, data_length,gps_time, detectors,Nmod, bppe, generation_method,
			statfilename,"",autocorrfile);	
	std::cout<<"ended"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*n_steps);
	for (int j =0; j<n_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
	}
	write_file(chainfile, output_transform, n_steps, dimension);
	//output hottest chain too
	chainfile = "testing/data/mcmc_output_hot.csv";
	for(int j = 0; j<n_steps;j++){
			//output_transform[j][0]=std::exp(output[0][j][0])/mpc_sec;
			output_transform[j][0]=std::exp(output[chain_N-1][j][0]);
			output_transform[j][1]=output[chain_N-1][j][1];
			output_transform[j][2]=output[chain_N-1][j][2];
			output_transform[j][3]=output[chain_N-1][j][3];
	}
	write_file(chainfile, output_transform, n_steps, dimension);

	deallocate_3D_array(output, chain_N, n_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< n_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	free(frequencies );
	delete [] detectors;
	free(data_length);
}
void test8()
{
	int length = 5000;
	int num_detectors =3;
	std::string *detectors = new std::string[num_detectors];//(std::string*)malloc(sizeof(std::string)*50*num_detectors);
	detectors[0] = "Hanford";
	detectors[1] = "Livingston";
	detectors[2] = "Virgo";

	double gps_time = 1135136350.6;//TESTING -- gw151226

	std::complex<double> **data= (std::complex<double>**)malloc(
			sizeof(std::complex<double>*)*num_detectors);
	double **psd = (double **)malloc(sizeof(double *)*num_detectors);
	double **frequencies = (double **)malloc(sizeof(double *)*num_detectors);
	int *data_length= (int*)malloc(sizeof(int)*num_detectors);
	data_length[0] =length;
	data_length[1] =length;
	data_length[2] =length;

	for (int i =0; i<num_detectors; i++){
		data[i] = (std::complex<double> *)malloc(
			sizeof(std::complex<double>)*data_length[i]);
		
		psd[i] = (double *)malloc(sizeof(double)*data_length[i]);
		frequencies[i] = (double *)malloc(sizeof(double)*data_length[i]);
	}
	//#########################################################
	//Make trial data
	gen_params params;
	std::string generation_method = "IMRPhenomD";
	double RA = 2.;
	double DEC = -1.2;
	double chirpm = 20.78;
	double eta =.249;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.4;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	params.phic = .0;
	double tc = 1;
	params.Luminosity_Distance = 810.;
	params.NSflag = false;
	params.incl_angle = 0;
	params.sky_average=false;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	//############################################################
	double fhigh =400;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	//double *freq = (double *)malloc(sizeof(double) * length);
	double freq[length];

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	//############################################################
	//############################################################
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++){
		noise[i] = noise[i]*noise[i];
	}
	//############################################################

	params.phi = 0;
	params.theta = 0;
	params.tc = tc; 
	celestial_horizon_transform(RA, DEC, gps_time, "Hanford", &params.phi, &params.theta);
	fourier_detector_response(freq,length, waveformout,"Hanford","IMRPhenomD", &params);
	double temptheta = params.theta;
	for(int j = 0; j<data_length[0]; j++){
		frequencies[0][j] = freq[j];	
		psd[0][j] = noise[j];	
		data[0][j] = waveformout[j];	
	}
	double snr2 = pow(calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);

	celestial_horizon_transform(RA, DEC, gps_time, "Livingston", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Livingston"); 
	fourier_detector_response(freq,length, waveformout,"Livingston","IMRPhenomD", &params);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[1][j] = freq[j];	
		psd[1][j] = (noise[j]);	
		data[1][j] = waveformout[j];	
	}
	snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	celestial_horizon_transform(RA, DEC, gps_time, "Virgo", &params.phi, &params.theta);
	params.tc = tc+ DTOA(temptheta, params.theta, "Hanford","Virgo"); 
	fourier_detector_response(freq,length, waveformout,"Virgo","IMRPhenomD", &params);
	for(int j = 0; j<data_length[0]; j++){
		frequencies[2][j] = freq[j];	
		psd[2][j] = (noise[j]);	
		data[2][j] = waveformout[j];	
	}
	//#########################################################
	snr2+=pow( calculate_snr("Hanford_O1_fitted",waveformout, freq, length),2);
	std::cout<<"SNR of injection: "<<sqrt(snr2)<<std::endl;
	
	//#########################################################
	//MCMC options
	int dimension = 4;
	//double initial_pos[dimension]={log(chirpm*MSOL_SEC), eta, params.spin1[2],params.spin2[2]};
	double initial_pos[dimension]={log(10), .2, 0,0};
	//double initial_pos[dimension]={log(200),log(20), .15, 0,0};
	double *seeding_var = NULL;
	int N_steps = 15000;
	int chain_N= 8;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 10;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//double temp_step = 20./(chain_N);
	chain_temps[0]=1.;
	double c = 1.2;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] = c*chain_temps[i-1];
		//chain_temps[i] = (1.+i*temp_step);
	
	int Nmod = 0;
	int *bppe = NULL;
	
	std::string autocorrfile = "testing/data/auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/mcmc_output.csv";
	std::string statfilename = "testing/data/mcmc_statistics.txt";

	
	int numThreads = 20;
	bool pool = true;

	MCMC_MH_GW(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, 
			swp_freq, test_lp_GW,numThreads, pool,show_progress, 
			num_detectors, 
			data, psd, 
			frequencies, data_length, gps_time,detectors,Nmod, bppe, generation_method,
			statfilename, "" ,autocorrfile);	
	std::cout<<"ENDED"<<std::endl;

	double **output_transform=(double **)malloc(sizeof(double*)*N_steps);
	for (int j =0; j<N_steps; j++)
		output_transform[j] = (double *)malloc(sizeof(double)*dimension);

	for(int j = 0; j<N_steps;j++){
			output_transform[j][0]=std::exp(output[0][j][0]);
			output_transform[j][1]=output[0][j][1];
			output_transform[j][2]=output[0][j][2];
			output_transform[j][3]=output[0][j][3];
	}
	write_file(chainfile, output_transform, N_steps, dimension);
	

	deallocate_3D_array(output, chain_N, N_steps, dimension);
	for(int i =0; i< num_detectors; i++){
		free(data[i]);
		free(psd[i]);
		free(frequencies[i]);
	}
	for(int i =0; i< N_steps; i++){
		free(output_transform[i]);
	}
	free(output_transform);
	free(data);
	free(psd);
	//free(freq);
	free(frequencies );
	//free(detectors);
	delete [] detectors;
	free(data_length);
}
void test7()
{
	int dimension = 2;
	double initial_pos[2]={0,0.};
	double *seeding_var = NULL;

	
	int N_steps = 10000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 1;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	double temp_step = 500./(chain_N);
	for(int i =0; i < chain_N;  i ++)
		//chain_temps[i]=1.;
		chain_temps[i] = 1+ temp_step * i;
	//double chain_temps[chain_N] ={1};
	std::string autocorrfile = "testing/data/neil_auto_corr_mcmc.csv";
	std::string chainfile = "testing/data/neil_mcmc_output.csv";
	std::string statfilename = "testing/data/neil_mcmc_statistics.txt";
	int numThreads = 10;
	bool pool = true;
	
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj3,fisher_neil_proj3,statfilename,chainfile,autocorrfile );	
	//auto lambda = [](double *x, int dim){return log_neil_proj3(x,dim);};
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj3,NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile );	
	MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, test_lp, log_neil_proj3,NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile );	
	std::cout<<"ENDED"<<std::endl;

	//autocorrfile = "testing/data/neil_auto_corr_mcmc2.csv";
	//chainfile = "testing/data/neil_mcmc_output2.csv";
	//statfilename = "testing/data/neil_mcmc_statistics2.txt";
	//MCMC_MH(output, dimension, N_steps, chain_N, initial_pos,chain_temps, swp_freq, test_lp, log_neil_proj32,NULL,statfilename,chainfile,autocorrfile );	
	//std::cout<<"ENDED"<<std::endl;
	


	//write_file("testing/data/mcmc_output.csv", output[0],N_steps, dimension);
	//ofstream mcmc_out;
	//mcmc_out.open("testing/data/mcmc_output.csv");
	//mcmc_out.precision(15);
	////for(int i = 0;i<chain_N;i++){
	//for(int j = 0; j<N_steps;j++){
	//	//for(int k = 0; k<dimension; k++){
	//		mcmc_out<<output[0][j][0]<<" , "<<output[0][j][1]<<endl;
	//	//}
	//}
	////}
	//mcmc_out.close();

	deallocate_3D_array(output, chain_N, N_steps, dimension);
}

void test6()
{

	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 8;
	double chirpm = 49.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	//params.betappe = new double[1] ;
	//params.betappe[0]=-100.;
	//params.bppe  =new int[1];
	//params.bppe[0] = -3;
	//params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	//params.f_ref = 100;
	//params.phiRef = 1.0;
	
	double fhigh =300;
	double flow =15;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	int dimension = 7;
	int dimensionmcmc = 5;

	clock_t start7,end7;

	double **output = (double **)malloc(dimension * sizeof(**output));	
	double **output2 = (double **)malloc(dimension * sizeof(**output2));	
	double **output3 = (double **)malloc(dimensionmcmc * sizeof(**output3));	

	for (int i = 0;i<dimension;i++){
		output[i] = (double *)malloc(dimension*sizeof(double));
		output2[i] = (double *)malloc(dimension*sizeof(double));
	}
	for (int i = 0;i<dimensionmcmc;i++){
	
		output3[i] = (double*)malloc(dimensionmcmc*sizeof(double));
	}
	start7 = clock();
	fisher(freq, length, "MCMC_IMRPhenomD_single_detect","Hanford_O1_fitted", 
			output3, dimensionmcmc, &params );

	end7 = clock();

	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;

	cout.precision(5);
	for (int i = 0;i <dimensionmcmc;i++)
	{
		for (int j=0;j <dimensionmcmc; j++)
			cout<<output3[i][j]<<"   ";
		cout<<endl;
	}
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(5);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	start7 = clock();
	fisher_autodiff(freq, length, "IMRPhenomD","Hanford_O1_fitted", output2, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER autodiff: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output2[i][j]<<"   ";
		cout<<endl;
	}
	std::cout<<"fractional DIFF: "<<std::endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<(output2[i][j]-output[i][j])/output2[i][j]<<"   ";
		cout<<endl;
	}
	for(int i =0; i<dimension; i++)
	{
		free(output[i]);
		free(output2[i]);
	}
	for(int i =0; i<dimensionmcmc; i++)
	{
		free(output3[i]);
	}
	free(output);
	free(output2);
	free(output3);
	free(freq);
}
void test5()
{

	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 1000;
	params.mass1 = 200;
	params.mass2 = 50;
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=10.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.f_ref = 100;
	params.phiRef = 1.0;
	params.sky_average=true;
	
	double fhigh =100;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];

	int dimension =8;
	clock_t start7,end7;
	double **output = (double **)malloc(dimension * sizeof(**output));	
	for (int i = 0;i<dimension;i++)
		output[i] = (double *)malloc(dimension*sizeof(double));
	
	start7 = clock();
	fisher(freq, length, "ppE_IMRPhenomD_Inspiral","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	free(output);
	free(freq);
}
void test4()
{

	cout.precision(15);

	gen_params params;

	int length = 16000;
	double fhigh =20;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	double *freqnew = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++){
		freq[i]=flow+i*df;
		//freqnew[i] = freq[i]-.1;
		freqnew[i] = freq[i];
	}
	//for(int i=0;i<length;i++){
	//	cout<<freqnew[i]<<" "<<freq[i] <<endl;
	//}
	
	double chirpmass = 20;
	double eta = .2;
	params.mass1 = calculate_mass1(chirpmass,eta);
	params.mass2 = calculate_mass2(chirpmass,eta);
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	string method= "IMRPhenomPv2";
	complex<double> *waveformout = (complex<double> *)malloc(sizeof(complex<double>) * length);
	complex<double> *waveformoutcross = (complex<double> *)malloc(sizeof(complex<double>) * length);
	params.spin1[0] = 0.01;
	params.spin1[1] = 0;
	params.spin1[2] = .1;
	params.spin2[0] = 0.01;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 1.0;
	params.tc = .0;
	params.Luminosity_Distance = 100.2;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 1.2;
	params.theta =3.4;
	params.incl_angle = 1.3;
	//params.f_ref=10;
	//params.phiRef=0.;
	
	
	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,waveformoutcross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;

	double *noise = (double *)malloc(sizeof(double)*length);
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];
	
	//params.mass1 = 100;
	//params.mass2 = 5;
	fftw_outline plan;
	initiate_likelihood_function(&plan,length);
	
	int masslen = 10;
	//double chirp = calculate_chirpmass(params.mass1,params.mass2);
	//double eta = calculate_eta(params.mass1,params.mass2);
	double masses[masslen];
	for (int i =0; i <masslen;i++)
		masses[i] = (i+1.)*.1*chirpmass;
	double mass2 = params.mass2;
	cout<<"Mass2: "<<mass2<<std::endl;
	


	double *real = (double *)malloc(sizeof(double)*length);
	double *imag = (double *)malloc(sizeof(double)*length);
	for ( int i =0; i<length;i++)
	{
		real[i]=(waveformout[i]).real();
		imag[i]=(waveformout[i]).imag();
	}
	
	complex<double> *hplus_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *detector_response = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hcross_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hplus_old = (complex<double> *)malloc(sizeof(complex<double>)*length);
	double llnew, llold,sum;
	complex<double> q;
	for (int i =0;i<masslen;i++)
	{
		params.mass1 = calculate_mass1(masses[i],eta);
		params.mass2 = calculate_mass2(masses[i],eta);

		
		start = clock();
		llnew = maximized_Log_Likelihood(waveformout, noise,freqnew,length, 
					&params,"Hanford","IMRPhenomPv2",&plan);
		//llnew = maximized_Log_Likelihood(real,imag, noise,freqnew,length, 
					//&params,"Hanford","IMRPhenomD",&plan);
		end = clock();
		cout<<"logl new  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		start = clock();
		fourier_detector_response(freqnew,length, detector_response,"Hanford",
					"IMRPhenomPv2",&params);
		end = clock();
		cout<<"waveform Pv2  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		start = clock();
		llold = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, real,imag, 
					noise,  calculate_chirpmass(params.mass1,params.mass2), 
					calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
		end = clock();
		
		start = clock();
		q = Q(params.theta,params.phi,params.incl_angle);
		fourier_waveform(freq,length, hplus_old,
					"IMRPhenomD",&params);
		for (int i = 0; i<length;i++)
			hplus_old[i] = q * hplus_old[i];
		end = clock();
		cout<<"waveform old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		//cout<<"LOGLnew: "<<llnew<<endl;
		//cout<<"LOGLold: "<<llold<<endl;
		cout<<"chirpmass: "<<masses[i]<<" new: "<<llnew<<" old: "<<llold<<" diff ll: "<<(llold-llnew)/llold<<endl;
		for (int i =0 ; i< length;i ++)
			sum += abs((detector_response[i]-hplus_old[i])/hplus_old[i]);
		cout<<"Average Diff waveform plus: "<<sum/length<<endl;
	}



	free(hplus_new);
	free(hcross_new);
	free(hplus_old);
	free(real);
	free(imag);
	free(waveformout);
	free(waveformoutcross);
	free(freq);
	free(freqnew);
	free(detector_response);
	free(noise);
	delete [] params.betappe;
	delete [] params.bppe;

	deactivate_likelihood_function(&plan);	
}
void test3()
{
	gen_params params;
	IMRPhenomPv2<double> modeld;
	IMRPhenomPv2<adouble> modela;
	int length = 900;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomPv2";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout_plus[length];
	complex<double> waveformout_cross[length];
	params.spin1[0] = .0;
	params.spin1[1] = .0;
	params.spin1[2] = -.2;
	params.spin2[0] = .0;
	params.spin2[1] = 0.0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	//params.phi = M_PI/3.;
	//params.theta = M_PI/3;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	
	double freq[length];
	for(int i=0;i<length;i++)
		freq[i]=10.+i*1e-1;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout_plus,waveformout_cross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	delete [] params.betappe;
	delete [] params.bppe;
	
}
void test2()
{
	double alpha, epsilon;
	IMRPhenomPv2<double> modelP;
	alpha = modelP.alpha(2,3,1./2,.75);
	epsilon = modelP.epsilon(2,3,1./2,.75);
	cout<<alpha<<endl; 
	cout<<epsilon<<endl; 
	long factorial_num = factorial(15);
	cout<<factorial_num<<endl;

	double d = modelP.d(2,1,0,.4);
	cout<<d<<endl;
}
void test1()
{
	gsl_spline *Z_DL_spline_ptr;
	gsl_interp_accel *Z_DL_accel_ptr;
	initiate_LumD_Z_interp(&Z_DL_accel_ptr,&Z_DL_spline_ptr);
	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 5000;
	double chirpm = 49.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	params.betappe = new double[1] ;
	params.betappe[0]=4*log(  1000/3.e5);
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	params.Z_DL_accel_ptr = Z_DL_accel_ptr;
	params.Z_DL_spline_ptr = Z_DL_spline_ptr;
	//params.f_ref = 30.5011;
	//params.phiRef =58.944425/2.;
	
	double fhigh =200;
	double flow =10;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	end=clock();
	
	clock_t  start2, end2;
	start2 = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	end2=clock();

	clock_t  start3, end3;
	start3 = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	end3=clock();

	ofstream ampfile;
	ampfile.open("testing/data/amplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	ofstream phasefile;
	phasefile.open("testing/data/phase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++)
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
	phasefile.close();

	ofstream wavefilereal;
	wavefilereal.open("testing/data/real_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	ofstream wavefileimag;
	wavefileimag.open("testing/data/imag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING amp: "<<(double)(end2-start2)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING phase: "<<(double)(end3-start3)/CLOCKS_PER_SEC<<endl;
	
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];

	double snr = calculate_snr("Hanford_O1_fitted", waveformout, freq, length);
	cout<<"SNR: "<<snr<<endl;
	//snr = data_snr_maximized_extrinsic(freq,length,waveformout,"Hanford_O1_fitted","IMRPhenomD",
	//		&params);
	//cout<<" FULL SNR: "<<snr<<endl;

	double logl;
	fftw_outline plan;
	clock_t  start4, end4;
	initiate_likelihood_function(&plan,length);
	start4 = clock();
	logl = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	end4 = clock();
	double logl2 = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, 1.2*calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	//params.mass1=300;
	//params.mass2=100;
	//params.spin1[2] = .9;
	//params.spin2[2] = -.2;
	start = clock();
	double logl3 = maximized_Log_Likelihood(waveformout, noise,freq,length, 
				&params,"Hanford","IMRPhenomD",&plan);
	end = clock();
	cout<<"logl new TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	start = clock();
	double logl4 = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, waveformout, 
				noise,  calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
	end = clock();
	cout.precision(15);
	cout<<"logl old old TIMING: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<logl<<endl;
	cout<<logl2<<endl;
	cout<<"LOGLnew: "<<logl3<<endl;
	cout<<"LOGLold: "<<logl4<<endl;
	deactivate_likelihood_function(&plan);	



	double real_data[length];
	double imag_data[length];
	double loglpy;
	for (int i = 0; i<length; i++)
	{
		real_data[i] = real(waveformout[i]);
		imag_data[i] = imag(waveformout[i]);
	}
	start4 = clock();
	loglpy = maximized_coal_log_likelihood_IMRPhenomD(freq, length, real_data, imag_data, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false);
	end4 = clock();
	cout<<"logl TIMING with setup: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<loglpy<<endl;

//###################################################################################################
	
	//method = "ppE_IMRPhenomD_Inspiral";
	method = "dCS_IMRPhenomD_log";
	//method = "EdGB_IMRPhenomD_log";
	clock_t  startppe, endppe;
	startppe = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	endppe=clock();
	cout<<"TIMING waveform ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;
	
	startppe = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	endppe=clock();
	cout<<"TIMING amplitude ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;

	startppe = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	endppe=clock();
	cout<<"TIMING phase ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;
	params.betappe[0]=2.;

	//ofstream ampfile;
	ampfile.open("testing/data/ppeamplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	//ofstream phasefile;
	phasefile.open("testing/data/ppephase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++){
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
		//std::cout<<phaseout[i]<<std::endl;
	}
	phasefile.close();

	//ofstream wavefilereal;
	wavefilereal.open("testing/data/ppereal_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	//ofstream wavefileimag;
	wavefileimag.open("testing/data/ppeimag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	
	
	
	
	
	int dimension = 7;
	
	double parameters[dimension] = {params.mass1,params.mass2,params.Luminosity_Distance,spin1[2],spin2[2],params.phic,params.tc};
	double **amp_derivative = (double**) malloc(dimension * sizeof(**amp_derivative));
	for (int i = 0; i<dimension;i++)
		amp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **phase_derivative = (double**) malloc(dimension * sizeof(**phase_derivative));
	for (int i = 0; i<dimension;i++)
		phase_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double spin1vec[3] = {0,0,parameters[3]};
	double spin2vec[3] = {0,0,parameters[4]};
	source_parameters<double> source_params;
	source_params = source_params.populate_source_parameters_old(parameters[0],
			parameters[1],parameters[2],spin1vec,spin2vec,parameters[5],
			parameters[6],false);

	lambda_parameters<double> lambda;
	modeld.assign_lambda_param(&source_params, &lambda);
	modeld.post_merger_variables(&source_params);
	source_params.f1 = 0.014/(source_params.M);
	source_params.f3 = modeld.fpeak(&source_params, &lambda);
	source_params.f1_phase = 0.018/(source_params.M);
	source_params.f2_phase = source_params.fRD/2.;
	source_params.bppe = new int[1];
	source_params.bppe[0] =-1;
	source_params.betappe = new double[1];
	source_params.betappe[0] = 10;
	source_params.Nmod=1;

	double A0 = source_params.A0; 
	double tc = source_params.tc;
	double phic = source_params.phic;
	double chirpmass = source_params.chirpmass;
	double symm = source_params.eta;
	double chi_s = source_params.chi_s;
	double chi_a = source_params.chi_a;


	IMRPhenomD<adouble> model;
	int amptapes[3] = {10,11,12};
	int phasetapes[3] = {13,14,15};
	model.amplitude_tape(&source_params, amptapes);
	model.phase_tape(&source_params, phasetapes);



	clock_t start5,end5;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modela.construct_amplitude_derivative(freq,length,dimension,amp_derivative, &source_params); 
	
	modela.construct_phase_derivative(freq,length,dimension,phase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream derivA;
	derivA.open("testing/data/deriv_amp.csv");
	derivA.precision(15);
	for(int i = 0;i<length;i++)
		derivA<<freq[i]<<','<<A0*amp_derivative[0][i]<<','<<amp_derivative[1][i]<<','<<amp_derivative[2][i]<<','<<chirpmass*amp_derivative[3][i]<<','<<symm*amp_derivative[4][i]<<','<<amp_derivative[5][i]<<','<<amp_derivative[6][i]<<endl;
	derivA.close();
	ofstream derivp;
	derivp.open("testing/data/deriv_phase.csv");
	derivp.precision(15);
	for(int i = 0;i<length;i++)
		derivp<<freq[i]<<','<<A0*phase_derivative[0][i]<<','<<phase_derivative[1][i]<<','<<phase_derivative[2][i]<<','<<chirpmass*phase_derivative[3][i]<<','<<symm*phase_derivative[4][i]<<','<<phase_derivative[5][i]<<','<<phase_derivative[6][i]<<endl;
	derivp.close();

	
	clock_t start7,end7;
	double **output = (double **)malloc(dimension * sizeof(**output));	
	for (int i = 0;i<dimension;i++)
		output[i] = (double *)malloc(dimension*sizeof(double));
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(4);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}

	int dimensionppe = dimension +1;
	double **outputppe = (double **)malloc(dimensionppe * sizeof(**outputppe));	
	for (int i = 0;i<dimensionppe;i++)
		outputppe[i] = (double *)malloc(dimensionppe*sizeof(double));
	start7 = clock();
	//fisher(freq, length, "ppE_IMRPhenomD_IMR","Hanford_O1_fitted", outputppe, dimensionppe, 
	//			&params );

	//end7 = clock();
	//cout<<"TIMING: FISHER ppE: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	//for (int i = 0;i <dimensionppe;i++)
	//{
	//	for (int j=0;j <dimensionppe; j++)
	//		cout<<outputppe[i][j]<<"   ";
	//	cout<<endl;
	//}


	double **ppeamp_derivative = (double**) malloc(dimensionppe * sizeof(**ppeamp_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppeamp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **ppephase_derivative = (double**) malloc(dimensionppe * sizeof(**ppephase_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppephase_derivative[i] = (double *)malloc(length * sizeof(double)); 


	ppE_IMRPhenomD_Inspiral<double> modelppe;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modelppe.construct_amplitude_derivative(freq,length,dimensionppe,ppeamp_derivative, &source_params); 
	
	modelppe.construct_phase_derivative(freq,length,dimensionppe,ppephase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad ppE: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream ppederivA;
	ppederivA.open("testing/data/ppederiv_amp.csv");
	ppederivA.precision(15);
	for(int i = 0;i<length;i++)
		ppederivA<<freq[i]<<','<<A0*ppeamp_derivative[0][i]<<','<<ppeamp_derivative[1][i]<<','<<ppeamp_derivative[2][i]<<','<<chirpmass*ppeamp_derivative[3][i]<<','<<symm*ppeamp_derivative[4][i]<<','<<ppeamp_derivative[5][i]<<','<<ppeamp_derivative[6][i]<<','<<ppeamp_derivative[7][i]<<endl;
	ppederivA.close();
	ofstream ppederivp;
	ppederivp.open("testing/data/ppederiv_phase.csv");
	ppederivp.precision(15);
	for(int i = 0;i<length;i++)
		ppederivp<<freq[i]<<','<<A0*ppephase_derivative[0][i]<<','<<ppephase_derivative[1][i]<<','<<ppephase_derivative[2][i]<<','<<chirpmass*ppephase_derivative[3][i]<<','<<symm*ppephase_derivative[4][i]<<','<<ppephase_derivative[5][i]<<','<<ppephase_derivative[6][i]<<','<<ppephase_derivative[7][i]<<endl;
	ppederivp.close();
	
	
	
	
	for (int i =0;i<dimension;i++)
	{
		free( amp_derivative[i]);
		free( phase_derivative[i]);
		free(output[i]);
	}
	for (int i =0;i<dimensionppe;i++)
	{
		free(outputppe[i]);
		free( ppeamp_derivative[i]);
		free( ppephase_derivative[i]);
	}
	free(amp_derivative);
	free(freq);
	free(phase_derivative);
	free(ppeamp_derivative);
	free(ppephase_derivative);
	free(output);
	free(outputppe);
	delete [] params.betappe;
	delete [] params.bppe;
	delete [] source_params.betappe;
	delete [] source_params.bppe;
	free_LumD_Z_interp(&Z_DL_accel_ptr, &Z_DL_spline_ptr);
}

void fisher_neil_proj3 (double *pos,int dimension, double **fisher)
{
	//int alpha = (int)(gsl_rng_uniform(g)*1e7);
 	adouble* x = new adouble[dimension];
 	adouble y = 1;  
 	double out =1;
 	trace_on(1);
 	for (int i =0; i< dimension; i++){
 	        x[i]<<= pos[i];
 	}
 	y =-1* log(dist(x, dimension));
 	y>>=out;
 	delete[] x;
 	trace_off();
 	hessian(1,dimension,pos,fisher);
	for (int i = 0 ; i<dimension; i++){
        	for (int j=0;j<i;j++){
        	        if (i!=j) fisher[j][i] =fisher[i][j];
        	}
	}

}
adouble dist(adouble *pos, int dimension){
        adouble x = pos[0];
        adouble y = pos[1];
        adouble exponent_1 = - pow(x,2) - pow(9 + 4*pow(x,2) + 8*y , 2);
        adouble exponent_2 = - 8*pow(x,2) - 8*pow(y - 2, 2);
        adouble out =( 16/(3 * M_PI) ) * ( exp(exponent_1) + 0.5 * exp(exponent_2) ); 
 
        return out;
}
double log_neil_proj3_nts (double *c,int dim, int chainid)
{
	return log_neil_proj3(c,dim);
}
double log_neil_proj3 (double *c,int dim)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_neil_proj32 (double *c,int dim)
{
	double x = c[0]*c[1];
	double y = c[0]/c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_student_t (double *x,int dim){

	double  mu=1, nu=3,  sigma=1;
        double g1 = gsl_sf_gamma( (nu + 1) / 2 ) ;
        double g2 = gsl_sf_gamma( (nu/2) );
        double parenth = 1 + (1/nu) *pow( (x[0] - mu) / sigma, 2 );
        return log(g1 / (g2 * sqrt(nu * M_PI ) * sigma ) * pow(parenth,-(nu+1)/2    ));
}
void test_fisher(double *pos, int dim, double **fisher)
{
	fisher[0][0] = .5;
}
double test_ll(double *pos, int dim)
{
	//std::cout<<"LL"<<std::endl;
	//std::cout<<"Pos in LL: "<<pos[0]<<std::endl;
	return -pos[0]*pos[0]/(4.);
	//return  0;
}
double test_lp_nts(double *pos, int dim, int chain_id){
	return test_lp(pos,dim);
}
double test_lp(double *pos, int dim)
{
	//return 0;
	return -pos[0]*pos[0]/(10.)- pos[1]*pos[1]/20.;
}	
double test_lp_GW(double *pos, int dim, int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	//if (std::exp(pos[0])/MPC_SEC<50 || std::exp(pos[0])/MPC_SEC>1000){return a;}
	if (std::exp(pos[0])<2 || std::exp(pos[0])>100){return a;}
	else if ((pos[1])<.1 || (pos[1])>.249999){return a;}
	else if ((pos[2])<-.9 || (pos[2])>.9){return a;}
	else if ((pos[3])<-.9 || (pos[3])>.9){return a;}
	//else {return 0.;}
	else {return pos[0] ;}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10);}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10)-(pos[1]-.24)*(pos[1]-.24)/(2*.010);}
}
double test_lp_GW_Pv2(double *pos, int dim, int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	if (std::exp(pos[1])<2 || std::exp(pos[1])>100){return a;}
	else if ((pos[0])<-1 || (pos[0])>1){return a;}
	else if ((pos[2])<.1 || (pos[2])>.249999){return a;}
	else if ((pos[3])<0 || (pos[3])>.9){return a;}
	else if ((pos[4])<0 || (pos[4])>.9){return a;}
	else if ((pos[5])<-1 || (pos[5])>1){return a;}
	else if ((pos[6])<-1 || (pos[6])>1){return a;}
	//else {return 0.;}
	else {return pos[0] ;}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10);}
	//else {return log(std::exp(pos[0])/MSOL_SEC)-(std::exp(pos[0])/MSOL_SEC-30)*(std::exp(pos[0])/MSOL_SEC-30)/(2*10)-(pos[1]-.24)*(pos[1]-.24)/(2*.010);}
}
double test_lp_GW_7dim(double *pos, int dim, int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	if (std::exp(pos[0])<50 || std::exp(pos[0])>1000){return a;}
	if ((pos[1])<0 || (pos[1])>4){return a;}
	if ((pos[2])<0 || (pos[2])>2*M_PI){return a;}
	if (std::exp(pos[3])<2 || std::exp(pos[3])>100){return a;}
	if ((pos[4])<.1 || (pos[4])>.245){return a;}
	if ((pos[5])<-.9 || (pos[5])>.9){return a;}
	if ((pos[6])<-.9 || (pos[6])>.9){return a;}
	//else {return 0.;}
	//else {return log(std::exp(pos[3])*std::exp(pos[0])*std::exp(pos[0])*std::exp(pos[0]));}
	else {return pos[3]+4*pos[0];}
}
double test_lp_GW_DFull(double *pos, int dim, int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//Flat priors across physical regions
	//if ((pos[0])<0 || (pos[0])>M_PI){return a;}
	if ((pos[0])<-1 || (pos[0])>1){return a;}//cos \iota
	if ((pos[1])<0 || (pos[1])>2*M_PI){return a;}//RA
	if ((pos[2])<-M_PI/2. || (pos[2])>M_PI/2.){return a;}//DEC
	if (std::exp(pos[3])<100 || std::exp(pos[3])>10000){return a;}//DL
	if (std::exp(pos[4])<2 || std::exp(pos[4])>100){return a;}//chirpmass
	if ((pos[5])<.1 || (pos[5])>.249999){return a;}//eta
	if ((pos[6])<-.9 || (pos[6])>.9){return a;}//chi1 
	if ((pos[7])<-.9 || (pos[7])>.9){return a;}//chi2
	//else {return log(-sin(pos[0]))+pos[4]+pos[3];}
	else {return pos[4]+3*pos[3];}
	//else {return pos[4]+1*pos[3];}
}
double test_lp_GW_dCS(double *pos, int dim , int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	double alphasqmin = 0;
	double alphasqmax = 1e4;//>10^6 km in \sqrt{\alpha}
	//if((std::exp(pos[8])<alphamin) || std::exp(pos[8])>alphamax){return a;}
	//if((std::exp(pos[8])<alphamin) ){return a;}
	//if(((pos[8])<lnalphamin) ){return a;}
	if(((pos[8])<alphasqmin )|| (pos[8]>alphasqmax  )){return a;}
	
	//Uniform prior on \alpha^2, not \alpha
	else { return test_lp_GW_DFull(pos,dim,chain_id)-.75*std::log(pos[8]);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_GW_dCS_log(double *pos, int dim , int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//double lnalphamin = -500;//Corresponds to a limit on root alpha of 10^-46 m
	//double lnalphamin = -100;//Corresponds to a limit on root alpha of 4 mm
	double lnalphamin = -70;//Corresponds to a limit on root alpha of 7.5 m
	double lnalphamax = 0;
	//if((std::exp(pos[8])<alphamin) || std::exp(pos[8])>alphamax){return a;}
	//if((std::exp(pos[8])<alphamin) ){return a;}
	//if(((pos[8])<lnalphamin) ){return a;}
	if(((pos[8])<lnalphamin)|| (pos[8]>lnalphamax) ){return a;}
	
	//Uniform prior on \alpha^.5, not \alpha
	//else { return test_lp_GW_DFull(pos,dim,chain_id)+.25*pos[8];}
	//Uniform in ln \alpha^2
	else { return test_lp_GW_DFull(pos,dim,chain_id);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
double test_lp_GW_ppE(double *pos, int dim , int chain_id)
{
	double a = -std::numeric_limits<double>::infinity();
	//if(((pos[8])<lnalphamin)|| (pos[8]>lnalphamax) ){return a;}
	//Uniform prior on \alpha^.5, not \alpha
	//else { return test_lp_GW_DFull(pos,dim,chain_id)+.25*pos[8];}
	//Uniform in ln \alpha^2
	{ return test_lp_GW_DFull(pos,dim,chain_id);}
	//else { return test_lp_GW_DFull(pos,dim,chain_id);}
}
