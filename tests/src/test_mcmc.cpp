#include <gwat/waveform_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/mcmc_gw.h>
#include <gwat/detector_util.h>
#include <iostream>


void RT_ERROR_MSG();
int mcmc_standard_test(int argc, char *argv[]);
int mcmc_injection(int argc, char *argv[]);
double log_test (double *c,int dim,void *parameters);
double log_test_prior (double *c,int dim,void *parameters);
void fisher_test(double *c,int dim, double **fisher,  void *parameters);
double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters);
double chirpmass_eta_jac(double chirpmass, double eta);
double T_mcmc_gw_tool ;

int main(int argc, char *argv[])
{
	std::cout<<"TESTING MCMC CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"MCMC student t"<<std::endl;
		return mcmc_standard_test(argc,argv);
	}
	else if(runtime_opt == 1){
		std::cout<<"MCMC Injected GW"<<std::endl;
		return mcmc_injection(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int mcmc_injection(int argc, char *argv[])
{
	gen_params injection;
	injection.mass1 = 36;
	injection.mass2 = 29;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	injection.Luminosity_Distance = 450;
	injection.psi = 1.;
	injection.phiRef = 2.;
	injection.f_ref = 20.;
	double tc_ref = 3.;
	injection.RA = 1.5;
	injection.DEC = -1.4;
	injection.spin1[2] = .01;
	injection.spin2[2] = -.01;
	injection.incl_angle = -M_PI +.01;
	//injection.incl_angle = -M_PI +.51;
	double gps = 1126259462.4;
	injection.gmst = gps_to_GMST_radian(gps);
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;
	
	int detect_number = 2;
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	std::string SN[3] = {"Hanford_O1_fitted","Hanford_O1_fitted","Hanford_O1_fitted"};
	std::string injection_method = "IMRPhenomD";
	double fmin = 5;
	double fmax =1024;
	T_mcmc_gw_tool= 4;
	double deltaf = 1./T_mcmc_gw_tool;
	int length = (fmax-fmin)/deltaf;
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];
	double total_snr = 0;
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		for(int j =0 ; j<length; j++){
			freq[i][j] = fmin + j*deltaf;
		}
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
		double deltat = DTOA_DETECTOR(injection.RA,injection.DEC,injection.gmst,detectors[0],detectors[i]);
		injection.tc = tc_ref + deltat;
		//injection.tc = tc_ref - deltat;
		fourier_detector_response(freq[i],length, data[i],detectors[i],injection_method, &injection, (double *)NULL);
		total_snr += pow_int( calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false), 2);
	}
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;

	int dim = 11;
	std::string recovery_method = "IMRPhenomD";
	int chains = 20;
	double temps[chains];
	double c = 1.3;
	temps[0]=1;
	for(int i = 1 ; i<chains; i++){
		temps[i]=temps[i-1]*c;
	}
	
	
	double initial_position[dim]= {injection.RA, sin(injection.DEC),injection.psi, cos(injection.incl_angle), injection.phiRef, tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	//double initial_position[dim]= {injection.RA, sin(injection.DEC), cos(injection.incl_angle), injection.phiRef, tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	write_file("data/injections.csv",initial_position,dim);
	initial_position[3] = -initial_position[3];
	initial_position[1] = 0;
	double *seeding = NULL;
	int swap_freq = 3;
	int threads = 10;
	bool pool = true;
	bool show_progress = true;
	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}
	int Nmod = 0;
	int *bppe = NULL;
	std::string stat_file = "data/injection_stat.txt";
	std::string output_file = "data/injection_output.csv";
	std::string ll_file = "data/injection_ll.csv";
	std::string checkpoint_file = "data/injection_checkpoint.csv";
	std::string ac_file = "data/injection_ac.csv";
	//std::string ac_file = "";
	
	int samples = 100000;
	double ***output  = allocate_3D_array(chains, samples, dim);
	PTMCMC_MH_GW(output, dim , samples, chains, initial_position, seeding, temps, swap_freq, standard_log_prior_D, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, Nmod, bppe, recovery_method, stat_file, output_file, ac_file,ll_file, checkpoint_file);
	deallocate_3D_array(output, chains, samples, dim);
	//int samples = 2000;
	//double **output  = allocate_2D_array( samples, dim);
	//double t0 = 1000;
	//double nu = 500;
	//double corr_threshold = .01;
	//double corr_segments = 10;
	//double corr_converge_thresh = .1;
	//double corr_target_ac = .01;
	//std::string chain_distribution="double";
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dim , samples, chains, chains,initial_position, seeding, temps, swap_freq, t0,nu,corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_distribution,standard_log_prior_D, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, Nmod, bppe, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//deallocate_2D_array(output,  samples, dim);

	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] freq;
	delete [] psd;
	delete [] data;

	return 0;
}
int mcmc_standard_test(int argc, char *argv[])
{
	int dimension = 2;
	double initial_pos[2]={1,0.};
	double *seeding_var = NULL;
	int N_steps = 100000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.8;
	for(int i =1; i < chain_N/2;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	chain_temps[chain_N/2] = 1;
	for(int i =chain_N/2+1; i < chain_N;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output.csv";
	std::string statfilename = "data/mcmc_statistics.txt";
	std::string checkpointfile = "data/mcmc_checkpoint.csv";
	
	int numThreads = 10;
	bool pool = true;
	bool show_progress = true;
	
	PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;


	deallocate_3D_array(output, chain_N, N_steps, dimension);
		
	return 0;

}
double log_test (double *c,int dim,void *parameters)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_test_prior (double *c,int dim,void *parameters)
{
	return 1.;
}

double fisher11(double *c)
{
	double x = c[0];
	double y = c[1];

	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(x,2) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(-2 - 128*pow(x,2) - 16*(9 + 4*pow(x,2) + 8*y)) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y),2));
}
double fisher12(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*x + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*x*(-2 + y) - 
     16*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(9 + 4*pow(x,2) + 8*y)*(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y)));

}
double fisher22(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) - 128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2)) + 
     128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(-2 + y,2) + 256*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(9 + 4*pow(x,2) + 8*y,2));

}
void fisher_test(double *c,int dim, double **fisher,  void *parameters)
{
	
	fisher[0][0] = fisher11(c);
	fisher[1][1] = fisher22(c);
	fisher[1][0] = fisher12(c);
	fisher[0][1] = fisher12(c);
	return;
} 
	
double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>60 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<-.95 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<-.95 || (pos[10])>.95){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}

}
//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}

void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Standard"<<std::endl;
	std::cout<<"1 --- GW injection"<<std::endl;
}
