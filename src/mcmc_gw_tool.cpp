#include "mcmc_gw.h"
#include "io_util.h"
#include "util.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <limits>

/*! \file 
 *
 * Command line tool for analyzing LOSC GW data
 *
 * Runs the ``uncorrelated'' MCMC sampler with a generic prior
 *
 * See mcmc_sampler.cpp documentation for more explanation on the sampler
 *
 * See mcmc_gw.cpp for GW specific explanation
 *
 * Usage:
 * 	
 * 	mcmc_gw_tool /PATH/TO/PARAM/FILE
 *
 * See data/mcmc_gw_tool_param_template.dat for an example parameter file
 *
 * See data/sample_init_pos.csv for an example initial position file
 */

double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters);
double standard_log_prior_Pv2(double *pos, int dim, int chain_id,void *parameters);
int main(int argc, char *argv[])
{
	if(argc !=2){
		std::cout<<"ERROR -- A parameter file is required"<<std::endl;
		return 1;
	}
	std::string param_file(argv[1]);

	std::unordered_map<std::string, int> int_dict;
	std::unordered_map<std::string, double> dbl_dict;
	std::unordered_map<std::string, float> flt_dict;
	std::unordered_map<std::string, bool> bool_dict;
	std::unordered_map<std::string, std::string> str_dict;
	std::string input_param_file(argv[1]);
	int status = unpack_input_io_file(input_param_file, &int_dict, &str_dict, &dbl_dict, &flt_dict, &bool_dict);

	//Global
	int detector_N = int_dict["detector number"];
	std::cout<<"DETECTOR NUMBER: "<<detector_N<<std::endl;
	int samples = int_dict["samples"];
	std::cout<<"Samples: "<<samples<<std::endl;
	int max_thermo_chain_N = int_dict["max thermo chain number"];
	std::cout<<"Max thermo ensemble chain number: "<<max_thermo_chain_N<<std::endl;
	int t0 = int_dict["t0"];
	std::cout<<"t0: "<<t0<<std::endl;
	int nu = int_dict["nu"];
	std::cout<<"nu: "<<nu<<std::endl;
	
	std::string detectors[detector_N];
	std::string detector_files[detector_N];
	std::string psd_file = str_dict["PSD filepath"];
	std::cout<<"PSD file: "<<psd_file<<std::endl;
	for(int i = 0 ; i<detector_N; i++){
		detectors[i]= str_dict["detector name "+std::to_string(i)];
		detector_files[i]= str_dict["data file "+std::to_string(i)];
	}
	std::string generation_method = str_dict["generation method"];
	
	double data_time_length = dbl_dict["data length"];
	std::cout<<"data length: "<<data_time_length<<std::endl;
	int data_length=131075;
	//if((int)data_time_length == 32){
	//	data_length = 131072;
	//}	
	count_lines_LOSC_data_file(detector_files[0], &data_length);
	std::cout<<"data length (lines): "<<data_length<<std::endl;
	double gps_time = dbl_dict["gps"];
	std::cout<<"GPS time : "<<gps_time<<std::endl;
	std::string allocation_scheme = str_dict["allocation scheme"];
	std::cout<<"Allocation: "<<allocation_scheme<<std::endl;
	int swap_freq = int_dict["swap frequency"];
	std::cout<<"swap_freq: "<<swap_freq<<std::endl;
	int threads = int_dict["thread number"];
	std::cout<<"threads: "<<threads<<std::endl;
	int correlation_thresh = int_dict["correlation threshold"];
	int correlation_segs = int_dict["correlation segments"];
	double correlation_convergence_thresh = dbl_dict["correlation convergence threshold"];
	double ac_target = dbl_dict["autocorrelation target"];
	std::string output_file = str_dict["output data file"];
	std::string stat_file = str_dict["output stat file"];
	std::string check_file = str_dict["checkpoint file"];
	
	std::string initial_position_file, initial_checkpoint_file;
	int chain_N;
	bool continue_from_checkpoint=false;
	if(str_dict.find("initial checkpoint file") == str_dict.end()){
		//Original
		initial_position_file = str_dict["initial position file"];
		chain_N = int_dict["chain number"];
		std::cout<<"Chain number: "<<chain_N<<std::endl;
	}
	else{
		//Continue	
		continue_from_checkpoint=true;
		initial_checkpoint_file = str_dict["initial checkpoint file"];
		std::cout<<"INITIAL checkpoint file: "<<initial_checkpoint_file<<std::endl;
		
	}
	

	
	int psd_length ;
	count_lines_LOSC_PSD_file(psd_file, &psd_length);
	//########################################################################

	double **psd = allocate_2D_array(detector_N,psd_length);
	double **freqs = allocate_2D_array(detector_N,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*detector_N);
	for(int i =0; i<detector_N; i++){
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);
	}

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time, data, psd, freqs);
	std::cout<<"DATA loaded"<<std::endl;


	int *data_lengths= (int*)malloc(sizeof(int)*detector_N);
	for(int i = 0 ; i<detector_N; i++){
		data_lengths[i] =psd_length;
	}
	//########################################################################
	int dimension = 11;
	double **output;
	output = allocate_2D_array(samples, dimension );
	double chain_temps[chain_N];
	
	int Nmod = 0;
	int *bppe = NULL;
	bool pool = true;
	bool show_progress = true;
	//#########################################################

	if(continue_from_checkpoint){
		continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(initial_checkpoint_file,output, samples,  
				max_thermo_chain_N, chain_temps, 
				swap_freq, t0, nu, correlation_thresh, correlation_segs,
				correlation_convergence_thresh , ac_target,allocation_scheme, 
				standard_log_prior_D,threads, pool,show_progress,detector_N, 
				data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
				generation_method,stat_file,output_file, "",check_file);	

	}
	else{
		double **initial_position = new double*[1];
		initial_position[0] = new double[dimension];
		read_file(initial_position_file, initial_position,1,dimension);
		double *init_pos = initial_position[0];
		double *seeding_var = NULL;
		PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dimension, samples, chain_N, 
				max_thermo_chain_N, init_pos,seeding_var,chain_temps, 
				swap_freq, t0, nu, correlation_thresh, correlation_segs,
				correlation_convergence_thresh , ac_target,allocation_scheme, 
				standard_log_prior_D,threads, pool,show_progress,detector_N, 
				data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
				generation_method,stat_file,output_file, "",check_file);	
		delete [] initial_position[0]; delete [] initial_position;
	}


	//write_file(chainfile, output[0], n_steps, dimension);

	deallocate_2D_array(output, samples, dimension);
	free(data_lengths);
	//free_LOSC_data(data, psd,freqs, detector_N, length);
	deallocate_2D_array(psd,detector_N, psd_length);
	deallocate_2D_array(freqs,detector_N, psd_length);
	for(int i =0; i<detector_N; i++)
		free(data[i]);
	free(data);

	
	return 0;

}
double standard_log_prior_D(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if ((pos[5])<0 || (pos[5])>10){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<-.9 || (pos[9])>.9){return a;}//chi1 
	if ((pos[10])<-.9 || (pos[10])>.9){return a;}//chi2
	else {return pos[7]+3*pos[6] ;}

}
double standard_log_prior_Pv2(double *pos, int dim, int chain_id,void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	double q =m1/m2;
	double W = (3*q +4)/ ( 4*q*q +3*q);
	//Max values
	double chi1l = pos[9];
	double chi2l = pos[10];
	double chi1p = std::sqrt(1- chi1l*chi1l);
	double chi2p = std::sqrt(1- chi2l*chi2l);
	double chi_thresh=W*chi2p ;
	if(chi1p > W*chi2p){ chi_thresh =chi1p;}
	if(pos[11] > chi_thresh){ return a;}

	//Flat priors across physical regions
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota

	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if ((pos[5])<0 || (pos[5])>10){return a;}//PhiRef
	if (std::exp(pos[6])<10 || std::exp(pos[6])>10000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>100 || std::isnan(pos[4])){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<-.9 || (pos[9])>.9){return a;}//chi1 
	if ((pos[10])<-.9 || (pos[10])>.9){return a;}//chi2
	if ((pos[11])<0 || (pos[11])>.9){return a;}//chip
	if ((pos[12])<0 || (pos[12])>2*M_PI){return a;}//phip
	else {return pos[7]+3*pos[6];}
}
