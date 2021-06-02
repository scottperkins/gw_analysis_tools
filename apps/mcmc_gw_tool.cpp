#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "ppE_utilities.h"
#include "io_util.h"
#include "util.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <limits>

double T_mcmc_gw_tool;
double T_merger; 
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
 * See data/sample_config_files/mcmc_gw_tool_param_template.dat for an example parameter file
 *
 * See data/sample_init_pos.csv for an example initial position file
 *
 * NOTE: SkySearch generation method requires an initialization file with 11 dimensions. 
 * The extra parameters correspond to the injected intrinsic parameters -- exactly the format of PhenomD
 */

double **mod_priors;
double chirpmass_prior[2];
double mass1_prior[2];
double mass2_prior[2];
double tidal1_prior[2];
double tidal2_prior[2];
double DL_prior[2];
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_NRT(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2(double *pos,  mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_D_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_Pv2_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_skysearch(double *pos, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_Pv2_mod(double *pos, mcmc_data_interface *interface,void *parameters);
double chirpmass_eta_jac(double m1,double m2);
double chirpmass_q_jac(double chirpmass, double q);
double aligned_spin_prior(double chi);
int main(int argc, char *argv[])
{
	std::cout.precision(15);
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
		std::cout<<"Detector stream file: "<<detectors[i]<<" : "<<detector_files[i]<<std::endl;
	}
	std::string generation_method = str_dict["generation method"];
	
	double data_time_length = dbl_dict["data length"];
	std::cout<<"data length: "<<data_time_length<<std::endl;
	int data_length=131072;
	if((int)data_time_length == 32){
		data_length = 131072;
	}	
	else if((int)data_time_length == 4096){
		data_length = 16777216;
	}	
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
	int correlation_thresh =10; 
	int correlation_segs = 10;
	double correlation_convergence_thresh =0.2;
	double ac_target = 0.01;
	std::string output_file = str_dict["output data file"];
	std::string stat_file = str_dict["output stat file"];
	std::string check_file = str_dict["checkpoint file"];
	int dimension = int_dict["dimension"];
	int max_chunk_size = int_dict["Max chunk size"];
	std::cout<<"Max chunk size: "<<max_chunk_size<<std::endl;
	
	std::string initial_position_file="", initial_checkpoint_file="",initial_ensemble_position_file="";
	int chain_N;
	bool continue_from_checkpoint=false;
	if(str_dict.find("initial checkpoint file") == str_dict.end()){
		//Original
		if(str_dict.find("initial position file") != str_dict.end()){
			initial_position_file = str_dict["initial position file"];
			chain_N = int_dict["chain number"];
			std::cout<<"Chain number: "<<chain_N<<std::endl;
			std::cout<<"Initial position file: "<<initial_position_file<<std::endl;
		}
		if (str_dict.find("initial ensemble position file") != str_dict.end()){
			initial_ensemble_position_file = str_dict["initial ensemble position file"];
			chain_N = int_dict["chain number"];
			std::cout<<"Chain number: "<<chain_N<<std::endl;
			std::cout<<"Initial ensemble position file: "<<initial_ensemble_position_file<<std::endl;
		}
		if ( initial_position_file =="" && initial_ensemble_position_file ==""){
			std::cout<<"ERROR -- you need an initial checkpoint file, an initial position file, or an initial ensemble file"<<std::endl;
			exit(1);
		}
	}
	else{
		//Continue	
		continue_from_checkpoint=true;
		initial_checkpoint_file = str_dict["initial checkpoint file"];
		std::cout<<"INITIAL checkpoint file: "<<initial_checkpoint_file<<std::endl;
		chain_number_from_checkpoint_file(initial_checkpoint_file,&chain_N);
		std::cout<<"Chain number: "<<chain_N<<std::endl;
		
	}

	if(dbl_dict.find("Mass1 minimum") == dbl_dict.end()){
		mass1_prior[0]=1;
		mass1_prior[1]=80;
	}
	else{
		mass1_prior[0]=dbl_dict["Mass1 minimum"];
		mass1_prior[1]=dbl_dict["Mass1 maximum"];
	}
	if(dbl_dict.find("Mass2 minimum") == dbl_dict.end()){
		mass2_prior[0]=1;
		mass2_prior[1]=80;
	}
	else{
		mass2_prior[0]=dbl_dict["Mass2 minimum"];
		mass2_prior[1]=dbl_dict["Mass2 maximum"];
	}
	if(dbl_dict.find("Luminosity distance minimum") == dbl_dict.end()){
		DL_prior[0]=10;
		DL_prior[1]=5000;
	}
	else{
		DL_prior[0]=dbl_dict["Luminosity distance minimum"];
		DL_prior[1]=dbl_dict["Luminosity distance maximum"];
	}
	std::cout<<"Range of Mass1: "<<mass1_prior[0]<<" - "<<mass1_prior[1]<<std::endl;
	std::cout<<"Range of Mass2: "<<mass2_prior[0]<<" - "<<mass2_prior[1]<<std::endl;
	std::cout<<"Range of DL: "<<DL_prior[0]<<" - "<<DL_prior[1]<<std::endl;
	if(dbl_dict.find("tidal1 minimum") == dbl_dict.end()){
		tidal1_prior[0]=1;
		tidal1_prior[1]=500;
	}
	else{
		tidal1_prior[0]=dbl_dict["tidal1 minimum"];
		tidal1_prior[1]=dbl_dict["tidal1 maximum"];
	}
	if(dbl_dict.find("tidal2 minimum") == dbl_dict.end()){
		tidal2_prior[0]=1;
		tidal2_prior[1]=500;
	}
	else{
		tidal2_prior[0]=dbl_dict["tidal2 minimum"];
		tidal2_prior[1]=dbl_dict["tidal2 maximum"];
	}
	if(generation_method.find("NRT") != std::string::npos){
		std::cout<<"Range of tidal1: "<<tidal1_prior[0]<<" - "<<tidal1_prior[1]<<std::endl;
		std::cout<<"Range of tidal2: "<<tidal2_prior[0]<<" - "<<tidal2_prior[1]<<std::endl;
	}
	

	
	int psd_length ;
	count_lines_LOSC_PSD_file(psd_file, &psd_length);
	std::cout<<"Length of PSD: "<<psd_length<<std::endl;

	double **psd = allocate_2D_array(detector_N,psd_length);
	double **freqs = allocate_2D_array(detector_N,psd_length);
	std::complex<double> **data = (std::complex<double> **)malloc(sizeof(std::complex<double> *)*detector_N);
	for(int i =0; i<detector_N; i++){
		data[i] = (std::complex<double>*)malloc(sizeof(std::complex<double>)*psd_length);
	}
	
	double post_merger_duration = 2;
	if( dbl_dict.find("Post merger signal duration") != dbl_dict.end()){
		post_merger_duration = dbl_dict["Post merger signal duration"];
	}
	std::cout<<"Poster merger signal duration: "<<post_merger_duration<<std::endl;

	allocate_LOSC_data(detector_files, psd_file, detector_N, psd_length, data_length, gps_time,post_merger_duration, data, psd, freqs);
	//########################################################################
	//TEST
	//for(int i = 0 ; i<detector_N; i++){
	//	for(int j =0 ; j<psd_length; j++){
	//		psd[i][j]=unpack1[j][i+1];
	//		data[i][j]=std::complex<double>(unpack2[j][2*i+1],unpack2[j][2*i+2]);
	//		freqs[i][j]= unpack1[j][0];
	//		//std::cout<<unpack1[j][1]<<std::endl;
	//		//std::cout<<unpack2[j][4]<<std::endl;
	//		//std::cout<<unpack1[1][j]<<std::endl;
	//	}	
	//}
	//freqs[0][0]=.0001;
	//freqs[1][0]=.0001;
	//########################################################################
	std::cout<<"DATA loaded"<<std::endl;

	T_mcmc_gw_tool = 1./(freqs[0][2]-freqs[0][1]);
	T_merger = T_mcmc_gw_tool - post_merger_duration;
	double df = 1./T_mcmc_gw_tool;
	//std::cout<<df<<std::endl;
	//for(int i = 0 ; i<detector_N; i++){
	//	for(int j =0 ; j<psd_length; j++){
	//		psd[i][j]*=df;
	//	}	
	//}

	std::cout<<"Total time: "<<T_mcmc_gw_tool<<std::endl;
	int *data_lengths= (int*)malloc(sizeof(int)*detector_N);
	for(int i = 0 ; i<detector_N; i++){
		data_lengths[i] =psd_length;
	}


	//double **output_test = allocate_2D_array(psd_length, 10 );
	//
	//for(int i = 0 ; i<psd_length; i++)
	//{
	//	output_test[i][0] = freqs[0][i];	
	//	output_test[i][1] = psd[0][i];	
	//	output_test[i][2] = psd[1][i];	
	//	output_test[i][3] = psd[2][i];	
	//	output_test[i][4] =real( data[0][i]);	
	//	output_test[i][5] =imag( data[0][i]);	
	//	output_test[i][6] =real( data[1][i]);	
	//	output_test[i][7] =imag( data[1][i]);	
	//	output_test[i][8] =real( data[2][i]);	
	//	output_test[i][9] =imag( data[2][i]);	
	//}
	//write_file("data/processed_data.csv",output_test,psd_length,10);
	//deallocate_2D_array(output_test,psd_length,10);
	

	//#############################################################
	//#############################################################
	
	//double **whitened = allocate_2D_array(data_lengths[0],10);
	//for(int i = 0 ; i<data_lengths[0]; i++){
	//	whitened[i][0]= freqs[0][i];
	//	whitened[i][1]=psd[0][i];
	//	whitened[i][2]=psd[1][i];
	//	whitened[i][3]=psd[2][i];
	//	whitened[i][4]=real(data[0][i])/std::pow(psd[0][i],.5);
	//	whitened[i][5]=imag(data[0][i])/std::pow(psd[0][i],.5);
	//	whitened[i][6]=real(data[1][i])/std::pow(psd[1][i],.5);
	//	whitened[i][7]=imag(data[1][i])/std::pow(psd[1][i],.5);
	//	whitened[i][8]=real(data[2][i])/std::pow(psd[2][i],.5);
	//	whitened[i][9]=imag(data[2][i])/std::pow(psd[2][i],.5);
	//}	
	//write_file("data/whitened_data.csv",whitened,data_lengths[0],10);
	//deallocate_2D_array(whitened, data_lengths[0],10);
	
	//#############################################################
	//#############################################################

	//########################################################################
	double **output;
	output = allocate_2D_array(samples, dimension );
	double chain_temps[chain_N];
	
	int Nmod = 0;
	int gNmod_phi = 0;
	int gNmod_sigma = 0;
	int gNmod_beta = 0;
	int gNmod_alpha = 0;
	int *gIMR_phii = NULL;
	int *gIMR_sigmai = NULL;
	int *gIMR_betai = NULL;
	int *gIMR_alphai = NULL;
	double *bppe = NULL;
	std::cout<<"Generation method: "<<generation_method<<std::endl;
	if(generation_method.find("ppE") != std::string::npos || check_theory_support(generation_method)){
		Nmod = int_dict["Number of modifications"];
		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
		bppe= new double[Nmod];
		mod_priors = new double*[Nmod];
		for(int i =0; i<Nmod ; i++){
			bppe[i] = dbl_dict["ppE b "+std::to_string(i)];
			mod_priors[i]= new double[2];
			std::cout<<i<<" : "<<bppe[i]<<std::endl;
			mod_priors[i][0] = dbl_dict["ppE beta "+std::to_string(i)+" minimum"];
			mod_priors[i][1] = dbl_dict["ppE beta "+std::to_string(i)+" maximum"];
			std::cout<<"Min"<<" : "<<mod_priors[i][0]<<std::endl;
			std::cout<<"Max"<<" : "<<mod_priors[i][1]<<std::endl;
		}
		
	}
	//if(generation_method.find("dCS") != std::string::npos
	//|| generation_method.find("EdGB") != std::string::npos){
	//	if(generation_method.find("EdGB_GHO") != std::string::npos){
	//		Nmod = 2;
	//		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
	//		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
	//		bppe= new double[Nmod];
	//		mod_priors = new double*[Nmod];
	//		bppe[0] = -7;
	//		bppe[1] = -5;
	//		mod_priors[0]= new double[2];
	//		mod_priors[1]= new double[2];
	//		mod_priors[0][0] = dbl_dict["ppE beta "+std::to_string(0)+" minimum"];
	//		mod_priors[0][1] = dbl_dict["ppE beta "+std::to_string(0)+" maximum"];
	//		mod_priors[1][0] = dbl_dict["ppE beta "+std::to_string(1)+" minimum"];
	//		mod_priors[1][1] = dbl_dict["ppE beta "+std::to_string(1)+" maximum"];
	//		std::cout<<0<<" : "<<bppe[0]<<std::endl;
	//		std::cout<<"Min"<<" : "<<mod_priors[0][0]<<std::endl;
	//		std::cout<<"Max"<<" : "<<mod_priors[0][1]<<std::endl;
	//		std::cout<<0<<" : "<<bppe[1]<<std::endl;
	//		std::cout<<"Min"<<" : "<<mod_priors[1][0]<<std::endl;
	//		std::cout<<"Max"<<" : "<<mod_priors[1][1]<<std::endl;
	//	
	//	}
	//	else{
	//		Nmod = 1;
	//		std::cout<<"Number of ppE modifications: "<<Nmod<<std::endl;
	//		std::cout<<"ppE b parmeters: "<<Nmod<<std::endl;
	//		bppe= new double[Nmod];
	//		mod_priors = new double*[Nmod];
	//		if(generation_method.find("dCS") != std::string::npos){
	//			bppe[0] = -1;
	//		}
	//		if(generation_method.find("EdGB") != std::string::npos){
	//			bppe[0] = -7;
	//		}
	//		mod_priors[0]= new double[2];
	//		std::cout<<0<<" : "<<bppe[0]<<std::endl;
	//		mod_priors[0][0] = dbl_dict["ppE beta "+std::to_string(0)+" minimum"];
	//		mod_priors[0][1] = dbl_dict["ppE beta "+std::to_string(0)+" maximum"];
	//	}
	//	
	//}
	if(generation_method.find("gIMR") != std::string::npos){
		gNmod_phi = int_dict["Number of phi modifications"];
		gNmod_sigma = int_dict["Number of sigma modifications"];
		gNmod_beta = int_dict["Number of beta modifications"];
		gNmod_alpha = int_dict["Number of alpha modifications"];
		int Nmod_tot = gNmod_phi+ gNmod_sigma+gNmod_beta+gNmod_alpha;
		std::cout<<"Number of general modifications: "<<Nmod_tot<<std::endl;
		std::cout<<"delta phi i: "<<Nmod<<std::endl;
		mod_priors = new double*[Nmod_tot];
		int ct = 0;
		if(gNmod_phi != 0){
			gIMR_phii= new int[gNmod_phi];
			for(int i =0; i<gNmod_phi ; i++){
				gIMR_phii[i] = int_dict["delta phi "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_phii[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta phi "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta phi "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_sigma != 0){
			gIMR_sigmai= new int[gNmod_sigma];
			for(int i =0; i<gNmod_sigma ; i++){
				gIMR_sigmai[i] = int_dict["delta sigma "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_sigmai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta sigma "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta sigma "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_beta != 0){
			gIMR_betai= new int[gNmod_beta];
			for(int i =0; i<gNmod_beta ; i++){
				gIMR_betai[i] = int_dict["delta beta "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_betai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta beta "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta beta "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		if(gNmod_alpha != 0){
			gIMR_alphai= new int[gNmod_alpha];
			for(int i =0; i<gNmod_alpha ; i++){
				gIMR_alphai[i] = int_dict["delta alpha "+std::to_string(i)+" i"];
				mod_priors[ct]= new double[2];
				std::cout<<i<<" : "<<gIMR_alphai[i]<<std::endl;
				mod_priors[ct][0] = dbl_dict["delta alpha "+std::to_string(i)+" minimum"];
				mod_priors[ct][1] = dbl_dict["delta alpha "+std::to_string(i)+" maximum"];
				ct++;
			}
		}
		
	}
	int total_mods = Nmod+gNmod_phi+gNmod_sigma+gNmod_beta+gNmod_alpha;
	bool pool = true;
	if(pool){
		debugger_print(__FILE__,__LINE__,"POOLING");
	}
	else{
		debugger_print(__FILE__,__LINE__,"NOT POOLING");
	}
	bool show_progress = true;
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = Nmod;
	mod_struct.bppe = bppe;
	mod_struct.gIMR_Nmod_phi = gNmod_phi;
	mod_struct.gIMR_phii = gIMR_phii;
	mod_struct.gIMR_Nmod_sigma = gNmod_sigma;
	mod_struct.gIMR_sigmai = gIMR_sigmai;
	mod_struct.gIMR_Nmod_beta = gNmod_beta;
	mod_struct.gIMR_betai = gIMR_betai;
	mod_struct.gIMR_Nmod_alpha = gNmod_alpha;
	mod_struct.gIMR_alphai = gIMR_alphai;
	if(bool_dict.find("NS Flag 0") != bool_dict.end()){
		mod_struct.NSflag1 = bool_dict["NS Flag 0"];
		if(mod_struct.NSflag1){
			std::cout<<"Object 0 is a NS"<<std::endl;
		}
	}
	if(bool_dict.find("NS Flag 1") != bool_dict.end()){
		mod_struct.NSflag2 = bool_dict["NS Flag 1"];
		if(mod_struct.NSflag2){
			std::cout<<"Object 1 is a NS"<<std::endl;
		}
	}

	//#########################################################
	if(generation_method.find("SkySearch") != std::string::npos){
		std::complex<double> *hplus=new std::complex<double>[data_lengths[0]];
		std::complex<double> *hcross=new std::complex<double>[data_lengths[0]];

		double **initial_position = new double*[1];
		initial_position[0] = new double[11];
		read_file(initial_position_file, initial_position,1,11);

		gen_params params ;
		for(int i = 0 ; i<11; i++){
			std::cout<<initial_position[0][i]<<std::endl;
		}
		params.mass1 = calculate_mass1(exp(initial_position[0][7]),initial_position[0][8]);
		params.mass2 = calculate_mass2(exp(initial_position[0][7]),initial_position[0][8]);
		params.spin1[2]=initial_position[0][9];
		params.spin2[2]=initial_position[0][10];
		params.phiRef = 0;
		params.tc = 0;
		params.f_ref =20;
		params.NSflag1=false;
		params.NSflag2=false;
		params.shift_phase=true;
		params.shift_time=true;
		params.equatorial_orientation=false;
		params.sky_average=false;
		params.Luminosity_Distance = exp(initial_position[0][6]);
		params.incl_angle = acos(initial_position[0][3]);

		//fourier_waveform(freqs[0],data_lengths[0],hplus, hcross, "IMRPhenomD",&params);
		waveform_polarizations<double> wp;
		wp.hplus = hplus;
		wp.hcross = hcross;
		fourier_waveform(freqs[0],data_lengths[0],&wp, "IMRPhenomD",&params);
	
	
		double *seeding_var = NULL;
		double **ensemble_initial_position=NULL;
		mcmc_sampler_output sampler_output(chain_N, dimension);
		SkySearch_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dimension, samples, chain_N, 
				max_thermo_chain_N, initial_position[0],seeding_var,ensemble_initial_position,chain_temps, 
				swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
				standard_log_prior_skysearch,threads, pool,show_progress,detector_N, 
				data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
				hplus,hcross,stat_file,output_file, "",check_file);	
		delete []  hplus;
		delete []  hcross;
		delete [] initial_position[0]; delete [] initial_position;
	}
	else{
	
		double(*lp)(double *param, mcmc_data_interface *interface, void *parameters);
		if(generation_method.find("IMRPhenomD") != std::string::npos && (dimension-total_mods) == 4){
			if(total_mods == 0){
				lp = &standard_log_prior_D_intrinsic;
			}
			else{
				lp = &standard_log_prior_D_intrinsic_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD_NRT") != std::string::npos && (dimension-total_mods) == 6){
			
			if(total_mods == 0){
				lp = &standard_log_prior_D_intrinsic_NRT;
			}
			else{
				lp = &standard_log_prior_D_intrinsic_NRT_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD") != std::string::npos && (dimension-total_mods) == 11){
			if(total_mods == 0){
				lp = &standard_log_prior_D;
			}
			else{
				lp = &standard_log_prior_D_mod;
			}
		}
		else if(generation_method.find("IMRPhenomD_NRT") != std::string::npos && (dimension-total_mods) == 13){
			if(total_mods == 0){
				lp = &standard_log_prior_D_NRT;
			}
			else{
				lp = &standard_log_prior_D_NRT_mod;
			}
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && (dimension-total_mods) == 8){
			if(total_mods == 0){
				lp = &standard_log_prior_Pv2_intrinsic;
			}
			else{
				lp = &standard_log_prior_Pv2_intrinsic_mod;
			}
		}
		else if(generation_method.find("IMRPhenomPv2") != std::string::npos && (dimension-total_mods) == 15){
			if(total_mods == 0){
				lp = &standard_log_prior_Pv2;
			}
			else{
				lp = &standard_log_prior_Pv2_mod;
			}
		}
		else{
			std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
			return 1;
		}
		//if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 11){
		//	lp = &standard_log_prior_D;
		//}
		//else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 15){
		//	lp = &standard_log_prior_Pv2;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomPv2") != std::string::npos  
		//	|| generation_method.find("gIMRPhenomPv2") !=std::string::npos||
		//	generation_method.find("PNSeries_ppE_IMRPhenomPv2") !=std::string::npos||
		//	generation_method.find("dCS_IMRPhenomPv2") != std::string::npos || 
		//	generation_method.find("EdGB_IMRPhenomPv2") != std::string::npos ||
		//	generation_method.find("EdGB_GHOv1_IMRPhenomPv2") != std::string::npos|| 
		//	generation_method.find("EdGB_GHOv2_IMRPhenomPv2") != std::string::npos|| 
		//	generation_method.find("EdGB_GHOv3_IMRPhenomPv2") != std::string::npos)
		//	&& dimension >= 16 ){
		//	lp = &standard_log_prior_Pv2_mod;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("gIMRPhenomD") != std::string::npos
		//	|| generation_method.find("PNSeries_ppE_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("dCS_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("EdGB_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("EdGB_GHOv1_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("EdGB_GHOv2_IMRPhenomD") != std::string::npos
		//	|| generation_method.find("EdGB_GHOv3_IMRPhenomD") != std::string::npos)
		//	&& dimension >= 11){
		//	lp = &standard_log_prior_D_mod;
		//}
		//else if(generation_method.find("IMRPhenomPv2") != std::string::npos && dimension == 8){
		//	lp = &standard_log_prior_Pv2_intrinsic;
		//}
		//else if(generation_method.find("IMRPhenomD") != std::string::npos && dimension == 4){
		//	lp = &standard_log_prior_D_intrinsic;
		//}
		//else if((generation_method.find("ppE_IMRPhenomD") != std::string::npos 
		//	|| generation_method.find("gIMRPhenomD") != std::string::npos) 
		//	&& dimension >= 4){
		//	lp = &standard_log_prior_D_intrinsic_mod;
		//}
		//else if( (generation_method.find("ppE_IMRPhenomPv2") != std::string::npos  
		//	|| generation_method.find("gIMRPhenomPv2") !=std::string::npos)
		//	&& dimension >= 9 ){
		//	lp = &standard_log_prior_Pv2_intrinsic_mod;
		//}
		//else if((generation_method.find("IMRPhenomD_NRT") != std::string::npos) && (dimension == 6)){
		//	lp = &standard_log_prior_D_intrinsic_NRT;

		//}
		//else if((generation_method.find("IMRPhenomD_NRT") != std::string::npos) && (dimension == 13)){
		//	lp = &standard_log_prior_D_NRT;

		//}
		//else{
		//	std::cout<<"ERROR -- wrong detector/dimension combination for this tool -- Check mcmc_gw for general support"<<std::endl;
		//	return 1;
		//}

		if(continue_from_checkpoint){
			mcmc_sampler_output sampler_output(chain_N, dimension);
			continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(initial_checkpoint_file,&sampler_output,output, samples,  
					max_thermo_chain_N, chain_temps, 
					swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
					lp,threads, pool,show_progress,detector_N, 
					data, psd,freqs, data_lengths,gps_time, detectors,&mod_struct,
					generation_method,stat_file,output_file, "",check_file);	

		}
		else{
			double **initial_position=NULL;
			initial_position = new double*[1];
			initial_position[0] = new double[dimension];
			double *seeding_var = NULL;
			double **ensemble_initial_position = NULL;
			if(initial_position_file != ""){
				read_file(initial_position_file, initial_position,1,dimension);
			}
			if(initial_ensemble_position_file != ""){
				ensemble_initial_position = new double*[chain_N];
				for(int i = 0 ; i<chain_N; i++){
					ensemble_initial_position[i] = new double[dimension];
				}
				read_file(initial_ensemble_position_file, ensemble_initial_position,chain_N,dimension);
				//for(int i = 0 ; i<chain_N; i++){
				//	for(int j = 0 ; j<dimension; j++){
				//		std::cout<<ensemble_initial_position[i][j]<<" ";
				//	}	
				//	std::cout<<std::endl;
				//}
			}
			std::cout<<"Running uncorrelated sampler "<<std::endl;
			mcmc_sampler_output sampler_output(chain_N, dimension);
			PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dimension, samples, chain_N, 
					max_thermo_chain_N, initial_position[0],seeding_var,ensemble_initial_position,chain_temps, 
					swap_freq, t0, nu,max_chunk_size,allocation_scheme, 
					lp,threads, pool,show_progress,detector_N, 
					//data, psd,freqs, data_lengths,gps_time, detectors,Nmod, bppe,
					data, psd,freqs, data_lengths,gps_time, detectors,&mod_struct,
					generation_method,stat_file,output_file, "",check_file);	
			sampler_output.create_data_dump(true, false, output_file);
			delete [] initial_position[0]; delete [] initial_position;
			if(initial_ensemble_position_file != ""){
				for(int i = 0 ; i<chain_N; i++){
					delete [] ensemble_initial_position[i]; 
				}
				delete [] ensemble_initial_position;
			}
		}

	}



	deallocate_2D_array(output, samples, dimension);
	free(data_lengths);
	deallocate_2D_array(psd,detector_N, psd_length);
	deallocate_2D_array(freqs,detector_N, psd_length);
	for(int i =0; i<detector_N; i++)
		free(data[i]);
	free(data);

	if(generation_method.find("ppE") != std::string::npos){
		delete [] bppe;	
		for(int i = 0 ; i<Nmod ; i++){
			delete [] mod_priors[i] ;
		}
		delete [] mod_priors;
	}
	
	return 0;

}
double standard_log_prior_D_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<-.95 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<-.95 || (pos[10])>.95){return a;}//chi2
	for(int i = 0 ; i<dim - 11; i++){
		if( pos[11+i] <mod_priors[i][0] || pos[11+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}
double standard_log_prior_D_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<dim - 13; i++){
		if( pos[13+i] <mod_priors[i][0] || pos[13+i] >mod_priors[i][1]){return a;}
	}
	return standard_log_prior_D_NRT(pos,interface, parameters);

}
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	if(pos[11]<tidal1_prior[0] || pos[11]>tidal1_prior[1]){return a;}
	if(pos[12]<tidal2_prior[0] || pos[12]>tidal2_prior[1]){return a;}
	return standard_log_prior_D(pos,interface, parameters);

}
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<-.95 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<-.95 || (pos[10])>.95){return a;}//chi2
	//return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}
double standard_log_prior_Pv2(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//####################
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//####################

	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<0 || (pos[9])>.95){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//phip
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];}
}
double standard_log_prior_Pv2_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//####################
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//####################

	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<DL_prior[0] || std::exp(pos[6])>DL_prior[1]){return a;}//DL
	if ((pos[9])<0 || (pos[9])>.95){return a;}//a1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//a2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//theta1
	if ((pos[12])<-1 || (pos[12])>1){return a;}//theta2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//phip
	if ((pos[14])<0 || (pos[14])>2*M_PI){return a;}//phip
	for(int i = 0 ; i<dim - 15; i++){
		if( pos[15+i] <mod_priors[i][0] || pos[15+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta))+3*pos[6];
}
double standard_log_prior_D_intrinsic_NRT(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	if(pos[4]<tidal1_prior[0] || pos[4]>tidal1_prior[1]){return a;}
	if(pos[5]<tidal2_prior[0] || pos[5]>tidal2_prior[1]){return a;}
	return standard_log_prior_D_intrinsic(pos,interface, parameters);

}
double standard_log_prior_D_intrinsic_NRT_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	for(int i = 0 ; i<dim - 6; i++){
		if( pos[6+i] <mod_priors[i][0] || pos[6+i] >mod_priors[i][1]){return a;}
	}
	return standard_log_prior_D_intrinsic_NRT(pos,interface,parameters) ;

}
double standard_log_prior_D_intrinsic(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<-.95 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<-.95 || (pos[3])>.95){return a;}//chi2
	//else {return log(chirpmass_eta_jac(chirp,eta)) ;}
	else {return log(aligned_spin_prior(pos[2]))+log(aligned_spin_prior(pos[3])) + log(chirpmass_eta_jac(chirp,eta));} 


}
double standard_log_prior_D_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<-.95 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<-.95 || (pos[3])>.95){return a;}//chi2
	for(int i = 0 ; i<dim - 4; i++){
		if( pos[4+i] <mod_priors[i][0] || pos[4+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta)) ;


}
double standard_log_prior_Pv2_intrinsic(double *pos, mcmc_data_interface *interface ,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<0 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>.95){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	if ((pos[7])<0 || (pos[7])>2*M_PI){return a;}//chi2
	return log(chirpmass_eta_jac(chirp,eta)) ;

}
double standard_log_prior_Pv2_intrinsic_mod(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	//###########
	double chirp = exp(pos[0]);
	double eta = pos[1];
	if (eta<.0 || eta>.25){return a;}//eta
	double m1 = calculate_mass1(chirp,eta );
	double m2 = calculate_mass2(chirp,eta );
	if(m1<mass1_prior[0] || m1>mass1_prior[1]){return a;}
	if(m2<mass2_prior[0] || m2>mass2_prior[1]){return a;}
	//###########
	if ((pos[2])<0 || (pos[2])>.95){return a;}//chi1 
	if ((pos[3])<0 || (pos[3])>.95){return a;}//chi2
	if ((pos[4])<-1 || (pos[4])>1){return a;}//chi1 
	if ((pos[5])<-1 || (pos[5])>1){return a;}//chi2
	if ((pos[6])<0 || (pos[6])>2*M_PI){return a;}//chi2
	if ((pos[7])<0 || (pos[7])>2*M_PI){return a;}//chi2
	for(int i = 0 ; i<dim - 8; i++){
		if( pos[8+i] <mod_priors[i][0] || pos[8+i] >mod_priors[i][1]){return a;}
	}
	return log(chirpmass_eta_jac(chirp,eta)) ;

}

double standard_log_prior_skysearch(double *pos, mcmc_data_interface *interface, void *parameters){

	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){return a;}//RA
	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//PhiRef
	if( pos[5] < (T_merger - .1) || pos[5] > (T_merger + .1)) { return a; }
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	else { return 3*pos[6];}
}

//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}

//Uniform in m1 and m2, transformed to lnM and q
double chirpmass_q_jac(double chirpmass, double q){
	return chirpmass*chirpmass/(pow(q/pow_int(q+1,2),1./5.) * q);
}

//For use with spin aligned system 
//-- mimics uniform magnitude and uniform cos \tilt, but mapped to chi
//-- Fit in mathematica 
double aligned_spin_prior(double chi){
	double a=0.0039132 , b= 3.95381;
	return a * exp(-b * abs(chi));	
}
