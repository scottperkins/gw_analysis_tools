#include <iostream>
#include <gwat/waveform_util.h>
//#include <gwat/waveform_generator.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_gw.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/EA_IMRPhenomD_NRT.h>

double T_merger;
double mass1_prior[2];
double mass2_prior[2];
double tidal1_prior[2];
double tidal2_prior[2];
double spin1_prior[2];
double spin2_prior[2];
double tidal_s_prior[2];
double EA_prior[6];
double DL_prior[2];
bool tidal_love;
bool tidal_love_error;
bool alpha_param;
bool EA_region1;

bool tidal_love_boundary_violation(double q,double lambda_s);
double EA_current_constraints(double *pos, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_D_NRT_EA(double *pos, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters);
double aligned_spin_prior(double chi);
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters);
double chirpmass_eta_jac(double chirpmass, double eta);
int main(int argc, char *argv[])
{
	std::cout<<"Injection EA"<<std::endl;
	

	//######################################
	//######################################
	//######################################
	//Injection Parameters
	//######################################
	//######################################
	//######################################

	gen_params injection;

	//Masses in M_SOL, plus useful other parameters
	injection.mass1 = 1.44;
	injection.mass2 = 1.29399; //These two together give a chirpmass of 1.188 (like in GW170817)
	//injection.mass1 = 1.81;
	//injection.mass2 = 1.514;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	double q = injection.mass2/injection.mass1;

	//Spins
	injection.spin1[2] = .003;
	injection.spin2[2] = -.002;
	injection.spin1[1] = .0;
	injection.spin2[1] = .0;
	injection.spin1[0] = .0;
	injection.spin2[0] = .0;

	//Extrinsic
	//injection.Luminosity_Distance = 159;
	injection.Luminosity_Distance = 63;
	injection.psi = .2;
	//injection.RA = .275;
        injection.RA = 3.42;
	injection.DEC = -.37;
	//injection.incl_angle = 0.785398;
	//injection.incl_angle = .51;
	injection.incl_angle = 2.532207345558998;
	double gps = 1187008882.4;
	injection.gmst = gps_to_GMST_radian(gps);

	//Tidal parameters -- either use tidal love or not
	//injection.tidal1 = 100;	
	//injection.tidal2 = 100;	
	//injection.tidal_s = .5*(injection.tidal1 + injection.tidal2);
	injection.tidal_love = true;
	injection.tidal_s = 242;
	
	//Trivial constants
	injection.phiRef = 2.;
	injection.f_ref = 20.;

	//Specific flags -- probably don't modify
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;

	//EA specific injection parameters set to small values of coupling constants
	//Determine which EA parameterization to use
	alpha_param = false;
	EA_region1 = false;

	if(alpha_param){
	  //Initial point when using alpha parameterization
	  //injection.alpha1_EA = -5.53584498257556007e-07;
	  injection.alpha2_EA = -6.5861481401721529e-08;
	  //injection.alpha1_EA = -0.087895; //distinct from GR injection
	  //injection.alpha2_EA = -0.00767322; //distinct from GR injection
	  //injection.alpha1_EA =  -0.185022041097348;
	  //injection.alpha2_EA = -0.0070102640180704;
	  injection.alpha1_EA = -0.245;
	  injection.cbarw_EA = 0.956938; 
	  //injection.cbarw_EA = 0.999640898572741; 
	  //injection.cbarw_EA = 0.163453;	  
	  injection.csigma_EA = 0;
	}
	else{
	  if(EA_region1){ //region of parameter space with ctheta = 3ca. 
	    injection.ca_EA = 1.76044705650715010e-06; 
	    injection.delta_ctheta_EA = .001; 
	    injection.cw_EA = 3.52089e-06;
	    injection.csigma_EA = 0;
	  }
	  else{
	    //This is a value of coupling constants that I know satisfies all the constraints in the prior
	    /*injection.ca_EA = 5.01510850153863e-07; 
	      injection.ctheta_EA = 1.50544346457309e-06; 
	      injection.cw_EA = 5.11795863509178; */
	    injection.csigma_EA = 0; //Note that because the code is now set to run with only 15 dimensions in EA, changing this value from zero shouldn't change anything. 
	    //injection.csigma_EA =  1.43323203548789e-16;

	    //Test of point with (c_V-1) ~ 1e-15 and (c_S -1) ~ 1e-15
	    injection.ca_EA = 1.76044705650715010e-06; 
	    injection.ctheta_EA = 5.28136e-06; 
	    injection.cw_EA = 3.52089e-06; 
	    
	    //Test of point at the edge of the priors for injection with EA
	    /*injection.ca_EA = 1.6e-06; 
	      injection.ctheta_EA = 5.0e-06; 
	      injection.cw_EA = 5.11795863509178; 
	      injection.csigma_EA = 0;*/
	  }
	}

	//ppE coefficients -- leave off for now
	//injection.Nmod = 1;
	//injection.bppe = new double[1];
	//injection.betappe = new double[1];
	//injection.bppe[0] = -1;
	//injection.betappe[0] = 0;
	
	//PhenomD_NRT injection
	std::string injection_method = "IMRPhenomD_NRT";
	//std::string injection_method = "EA_IMRPhenomD_NRT";
	
	//Detector properties
	int detect_number = 3;
	//std::string detectors[4] = {"Hanford","Livingston","Virgo","Kagra"};
        //std::string SN[4] = {"AdLIGODesign","AdLIGODesign","AdLIGODesign","KAGRA_pess"};
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//Frequency array properties
	double fmin = 10;
	double fmax =2048;
	//Length of total signal -- leave at 4 seconds for now
	double Tsignal = 4;
	//double Tsignal = 128; 
	double deltaF = 1./Tsignal;
	//Merger time -- 3/4 of total signal length. If using >4 seconds, change to Tsignal - 2
	T_merger= Tsignal*3./4.;
	//T_merger = Tsignal - 2; 
	int length = (int)((fmax-fmin)/deltaF);
	
	//Input PSD, freq, and data arrays -- these will actually be passed into the MCMC
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];

	//Set merger time of injection
	injection.tc=Tsignal-T_merger;
	
	//Populate freq, psd, and allocate array for data
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i] = new double[length];
			
		//Populate freq array
		for(int j =0 ; j<length; j++){
			freq[i][j] = fmin + j*deltaF;
		}
		//Populate PSD array
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
	}

	//Set data_length array 
	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}

	//Create coherent gravitational wave detector between all detectors
	create_coherent_GW_detection(detectors,detect_number, freq, data_lengths,  true, &injection, injection_method, data);
	
	//Calculate the SNR and print
	double total_snr=0;
	for(int i = 0 ; i<detect_number; i++){
		total_snr += pow_int( calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false), 2);
	}
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;
	//###################################################################
	//###################################################################
	//MCMC Settings
	//###################################################################
	//###################################################################

	int dim = 15;
	std::string recovery_method = "EA_IMRPhenomD_NRT";
	//std::string recovery_method = "IMRPhenomD_NRT";

	int ensemble = 20;
	//int ensemble = 60;
	//int chains = 180;
	//int ensemble = 3;
	//int ensemble = 15;
	int chains = 100;
        //int chains = 6;
	//int chains = 60;

	double temps[chains];
	double c = 1.2;
	temps[0]=1;
	for(int i = 1 ; i<chains; i++){
		temps[i]=temps[i-1]*c;
	}
	
	//Initial Guess  -- also the exact right answer.. Perks of injections!!
	//Needs to match dimension/model we're using -- e.g. if using IMRPhenomD_NRT with tidal love, 12 dimensions. If not, 13.
	double initial_position[dim]; 
	initial_position[0]= injection.RA; 
	initial_position[1]= sin(injection.DEC);
	initial_position[2]= injection.psi; 
	initial_position[3]= cos(injection.incl_angle); 
	initial_position[4]= injection.phiRef;
	initial_position[5]= T_merger;
	initial_position[6]= log(injection.Luminosity_Distance);
	initial_position[7]= log(chirpmass);
	initial_position[8]= eta;
	initial_position[9]= injection.spin1[2];
	initial_position[10]= injection.spin2[2];
	initial_position[11]= injection.tidal_s;
	/*
	double initial_position[dim]= {
	  injection.RA, 
	  sin(injection.DEC),
	  injection.psi, 
	  cos(injection.incl_angle), 
	  injection.phiRef, 
	  T_merger, 
	  log(injection.Luminosity_Distance),
	  log(chirpmass), 
	  eta, 
	  injection.spin1[2],
	  injection.spin2[2],
	  injection.tidal_s,
	  0.,
	  0.,
	  0.
	  //injection.csigma_EA
	  };*/
	
       	if(alpha_param){
	  initial_position[12]= injection.alpha1_EA;
	  initial_position[13]= injection.alpha2_EA;
	  initial_position[14]= injection.cbarw_EA;
	}
	else{
	  if(EA_region1){
	    initial_position[12]= injection.ca_EA;
	    initial_position[13]= injection.delta_ctheta_EA;
	    initial_position[14]= injection.cw_EA;
	  }
	  else{
	    initial_position[12]= injection.ca_EA;
	    initial_position[13]= injection.ctheta_EA;
	    initial_position[14]= injection.cw_EA;
	  }
	} 
	  
	//Write it out to file
	write_file("data/injections.csv",initial_position,dim);


	//MCMC parameters
	//Keep at NULL -- outdated parameter
	double *seeding = NULL;
	//Frequency of swaps -- every N steps attempt a swap
	int swap_freq = 5;
	//Number of threads to use
	//int threads = 64;
	int threads = 8;
	//Whether to pool or not, leave at true
	bool pool = true;
	//Whether to show progress or not, leave at true
	bool show_progress = true;
	//Input files
	//std::string checkpoint_file_start = "data/injection_checkpoint_start.csv"; 
	//Output files -- leave ll_file/ac_file as blank -- outdated
	std::string stat_file = "data/injection_stat.txt";
	std::string output_file = "data/injection_output.hdf5";
	std::string ll_file = "";
	std::string checkpoint_file = "data/injection_checkpoint.csv";
	//std::string ac_file = "data/injection_ac.csv";
	std::string ac_file = "";
	
	//Additional modifications
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = 0;
	//double b = -1;
	//mod_struct.bppe =&b ;

	//Set parameterization to use
	mod_struct.alpha_param = alpha_param;
	mod_struct.EA_region1 = EA_region1;
	
	//Use tidal love or not
	tidal_love = true;
	mod_struct.tidal_love = tidal_love;
        tidal_love_error = false;
	mod_struct.tidal_love_error = tidal_love_error; 
	
	//Samples, burn in parameters 	
	int samples = 2000;
	double t0 = 100;
	double nu = 100;


	//M1 prior range in M_SOL
	mass1_prior[0] = .5;
	mass1_prior[1] = 3;

	//M2 prior range in M_SOL
	mass2_prior[0] = .5;
	mass2_prior[1] = 3;

	//spin1 prior range 
	spin1_prior[0] = -0.01;
	spin1_prior[1] = .01;

	//spin2 prior range 
	spin2_prior[0] = -0.01;
	spin2_prior[1] = .01;

	//Luminosity prior range in MPC
	DL_prior[0] = 1;
	DL_prior[1] = 500;
	
	//If using tidal_love, only tidal_s is used. Otherwise, tidal 1 and tidal 2 are used.
	//tidal 1 prior range 
	tidal1_prior[0] = 1;
	tidal1_prior[1] = 5000;

	//tidal 2 prior range 
	tidal2_prior[0] = 1;
	tidal2_prior[1] = 5000;

	//tidal s prior range 
	tidal_s_prior[0] = 1;
	tidal_s_prior[1] = 5000;

	//EA coupling constants prior range
	if(alpha_param){
	  EA_prior[0]= -.25; //alpha1 min
	  //EA_prior[0]= -8.; //alpha1 min
	  EA_prior[1]= .25; //alpha1 max
	  //EA_prior[1]= 8.; //alpha1 max
	  EA_prior[2]= -.025; //alpha2 
	  //EA_prior[2]= -8.; //alpha2 min
	  EA_prior[3]= .025; //alpha2 max
	  //EA_prior[3]= 8.; //alpha2 max
	  EA_prior[4]= -1.;//+pow(10,-10.); //cbarw min 
	  EA_prior[5]=1.;//-pow(10,-10.); //cbarw max
	}
	else if(EA_region1){
	  EA_prior[0]= -3.; //ca min
	  EA_prior[1]= 3.; //ca max
	  EA_prior[2]= -.003; //delta_ctheta min 
	  EA_prior[3]=.003; //delta_ctheta max
	  EA_prior[4]= -3; //cw min 
	  EA_prior[5]=3; //cw max
	}
	else{
	  EA_prior[0]= -3.; //ca min
	  EA_prior[1]= 3.; //ca max
	  EA_prior[2]= -3.; //ctheta min
	  EA_prior[3]= 3.; //ctheta max
	  EA_prior[4]= -3.; //cw min 
	  EA_prior[5]= 3.; //cw max
	}

	/* Test initial point for prior compliance */
	mcmc_data_interface interface ;
	interface.max_dim = dim;
	
	double prior = standard_log_prior_D_NRT_EA(initial_position,&interface, (void *)nullptr);
	std::cout<<"Prior for initial point: "<<prior<<std::endl;
	exit(1);
	
	//Output structures
	double **output  = allocate_2D_array( samples, dim);
	std::string chain_distribution="double";
	int max_chunk_size = 1e6;
	mcmc_sampler_output sampler_output(chains,dim);
	
	//Run the sampler! -- If need to, change the prior as needed.
	//Run with EA recovery
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dim , samples, chains, ensemble,initial_position, seeding, (double**)NULL,temps, swap_freq, t0,nu,max_chunk_size,chain_distribution,standard_log_prior_D_NRT_EA, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//Run with EA recovery from checkpoint
	//continue_PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(checkpoint_file_start, &sampler_output,output, samples, ensemble,temps, swap_freq, t0,nu,max_chunk_size,chain_distribution,standard_log_prior_D_NRT_EA, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//Run with NRT recovery
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(&sampler_output,output, dim , samples, chains, ensemble,initial_position, seeding, (double**)NULL,temps, swap_freq, t0,nu,max_chunk_size,chain_distribution,standard_log_prior_D_NRT, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ll_file, checkpoint_file);



	//Cleanup
	deallocate_2D_array(output,  samples, dim);

	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] injection.bppe;
	delete [] injection.betappe;
	delete [] freq;
	delete [] psd;
	delete [] data;

	return 0;
}
double EA_current_constraints(double *pos, mcmc_data_interface *interface, void *parameters)
{
  //std::cout<<"Checking EA_constraints"<<std::endl; 

  int dim =  interface->max_dim;
  double a = -std::numeric_limits<double>::infinity();

  source_parameters<double> sp;
  double lnChirpmass = pos[7];//ln(M_sol)
  double eta = pos[8];
  
  sp.mass1 = calculate_mass1(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.mass2 = calculate_mass2(std::exp(lnChirpmass)*MSOL_SEC,eta);
  sp.M = sp.mass1 + sp.mass2;
  sp.alpha_param = alpha_param;
  sp.EA_region1 = EA_region1;
  
  //Setting tidal1 and tidal2 from tidal_s
  if(tidal_love)
    {
      IMRPhenomD_NRT<double> modelNRT;
      modelNRT.binary_love_relation(pos[11], tidal_love_error, &sp);
      //std::cout<<"tidal1= "<<sp.tidal1<<"tidal2= "<<sp.tidal2<<std::endl; 
    }
  sp.csigma_EA = 0;
  
  if(tidal_love){
    if(alpha_param){
      sp.alpha1_EA = pos[12]; //alpha1
      sp.alpha2_EA = pos[13]; //alpha2
      sp.cbarw_EA = pos[14]; //cbarw
    }
    else{
      if(EA_region1){
	sp.ca_EA = pos[12]; //ca
	sp.delta_ctheta_EA = pos[13];
	sp.cw_EA = pos[14]; //cw
      }
      else{
	sp.ca_EA = pos[12]; //ca
	sp.ctheta_EA = pos[13]; //ctheta
	sp.cw_EA = pos[14]; //cw
	//sp.csigma_EA = pos[15]; //csigma
      }
    }
  }
  else{
    if(alpha_param){
      sp.alpha1_EA = pos[13]; //alpha1
      sp.alpha2_EA = pos[14]; //alpha2
      sp.cbarw_EA = pos[15]; //cbarw
    }
    else{
      if(EA_region1){
	sp.ca_EA = pos[13]; //ca
	sp.delta_ctheta_EA = pos[14];
	sp.cw_EA = pos[15]; //cw
      }
      else{
	sp.ca_EA = pos[13]; //ca
	sp.ctheta_EA = pos[14]; //ctheta
	sp.cw_EA = pos[15]; //cw
	//sp.csigma_EA = pos[16]; //csigma
      }
    }
  }
  sp.EA_nan_error_message = false;
  EA_IMRPhenomD_NRT<double> model;
  model.pre_calculate_EA_factors(&sp);

  if(sp.ca_EA < 0 || sp.ca_EA > 2.){return a;}
  /* Throws out points with ca < 0 or ca >2 because these violate 
   * the positive energy condition for the spin-0 mode (scalar mode).   
   * See equation 40 of arXiv:gr-qc/0507059v3.
   */
  if(sp.cw_EA < (-sp.csigma_EA/(1. - sp.csigma_EA))){return a;}
  /* Throws out points with cw < -csigma/(1 - csigma) because these violate 
   * the positive energy condition for the spin-1 mode (vector mode).
   * Note that the positive energy condition for the spin-2 modes is always 
   * satisfied. See equation 40 of arXiv:gr-qc/0507059v3.
   */
  //std::cout<<"EA constraints test 1"<<std::endl; 
  if(sp.cTsq_EA < 0 || sp.cVsq_EA < 0 || sp.cSsq_EA < 0){return a;}
  /* Throws out points with speeds not greater than or equal to zero (these 
   * would produce gradient instabilities or ghosts)
   * arXiv:gr-qc/0402005 and arXiv:1108.1835
   */
  //std::cout<<"EA constraints test 2"<<std::endl; 
  
  if(isnan(sp.kappa3_EA))
    {
      //std::cout<<"kappa3:"<<sp.kappa3_EA<<std::endl; 
      return a; 
    }
  
  
  //std::cout<<"EA constraints test 3"<<std::endl; 

  
  //if(sp.cT_EA - 1. < -3*pow(10, -15.) || sp.cT_EA -1. > 7*pow(10, -16.)){return a;}
  /* Throws out points that don't obey the cT constraint from GW170817 and GRB170817A
   * arXiv:1710.05834 just in case a user feeds in a bad csigma. 
   */
  //std::cout<<"EA constraints test 4"<<std::endl; 

  if(sp.ctheta_EA < 0 || (sp.ctheta_EA + (8.*sp.ca_EA/7)) > (2./7.)){return a;}
  /*
   * Throws out points that violate Big Bang Nucleosynthesis constraints.
   *  arXiv:hep-th/0407149v3 
   */

  if(sp.s1_EA >= 1 || sp.s2_EA >= 1){return a;}
  
  //if(fabs(sp.alpha1_EA) > pow(10, -4.) || fabs(sp.alpha2_EA) > 4.*pow(10, -7.)){return a;}
  /* Throws out points that do not obey observational solar system constraints on 
   * alpha1 and alpha2
   * arXiv:1403.7377 and arXiv:gr-qc/0509114
   */
  //std::cout<<"EA constraints test 5"<<std::endl; 

  bool Cherenkov = true;
  if(Cherenkov)
    {
      bool violate = false;
      if(abs(sp.c13_EA * sp.c13_EA *(sp.c1_EA * sp.c1_EA + 2*sp.c1_EA*sp.c3_EA + sp.c3_EA * sp.c3_EA - 2*sp.c4_EA)/(2*sp.c1_EA*sp.c1_EA)) >= 7*pow(10, -32.))
	    //enforcing constraint from Eq. 4.7 of arXiv:hep-ph/0505211
	    {
	      violate = true;
	    }
	  
      if(violate){return a;}
      //std::cout<<"EA constraints test 6"<<std::endl; 
      
      if(sp.cS_EA < 1)
	{
	  if(fabs((sp.c2_EA + sp.c3_EA - sp.c4_EA)/sp.c1_EA) > pow(10, -22.))
	    {
	      if((sp.c3_EA - sp.c4_EA)*(sp.c3_EA - sp.c4_EA)/fabs(sp.c14_EA) >= pow(10, -30.)){return a;}
	      //enforcing constraint from Eq.4.15 of arXiv:hep-ph/0505211
	    }
	  //std::cout<<"EA constraints test 7"<<std::endl; 
	  if(fabs((sp.c4_EA - sp.c2_EA - sp.c3_EA)/sp.c1_EA) >= 3*pow(10,-19.)){return a;}
	  //enforcing constraint from Eq.5.14 of arXiv:hep-ph/0505211
	  
	}
    }
  

  //std::cout<<"Made it to end of EA constraints"<<std::endl;
 /* 
  double sigma = 1.021*pow(10, -5.); 
  double mu = -0.563*pow(10, -5.);
  double prob;
  
  prob = exp(-(1./2.)*((sp.alpha1_EA - mu)*(sp.alpha1_EA - mu))/(sigma*sigma));
  */

  sp.EA_nan_error_message = false;
  model.EA_check_nan(&sp);
/*
  if(sp.s1_EA >1)
    {
      std::cout<<"s1 = "<<sp.s1_EA<<std::endl; 
    }
  if(sp.s2_EA >1)
    {
      std::cout<<"s2 = "<<sp.s2_EA<<std::endl; 
    }
  */
  //return log(prob);
  return 0; 
}
//Prior for EA with NRT
double standard_log_prior_D_NRT_EA(double *pos, mcmc_data_interface *interface, void *parameters)
{
  //std::cout<<"Starting prior"<<std::endl; 
  //std::cout<<pos[12]<<" "<<pos[13]<<" "<<pos[14]<<" "<<pos[15]<<std::endl;
  int dim =  interface->max_dim;
  double a = -std::numeric_limits<double>::infinity();
  if(tidal_love){
    if( pos[12] <EA_prior[0] || pos[12] >EA_prior[1]){return a;} //ca or alpha1
    if( pos[13] <EA_prior[2] || pos[13] >EA_prior[3]){return a;} //ctheta or alpha2
    if( pos[14] <EA_prior[4] || pos[14] >EA_prior[5]){return a;} //cw or cbarw
      //if( pos[15] <EA_prior[6] || pos[15] >EA_prior[7]){return a;} //csigma 
  }
  else{
    if( pos[13] <EA_prior[0] || pos[13] >EA_prior[1]){return a;} //ca or alpha1
    if( pos[14] <EA_prior[2] || pos[14] >EA_prior[3]){return a;} //ctheta or alpha2
    if( pos[15] <EA_prior[4] || pos[15] >EA_prior[5]){return a;} //cw or cbarw
    //if( pos[16] <EA_prior[6] || pos[16] >EA_prior[7]){return a;} //csigma
  }
  
  double NS = standard_log_prior_D_NRT(pos,interface, parameters);
  if(NS == a){return a;}
  double EA_constraints =  EA_current_constraints(pos, interface, parameters);
  if (EA_constraints ==a){return a;}
  //double EA_constraints = 0; 
  //std::cout<<"Made it to the end of the prior"<<std::endl; 
  
  return EA_constraints + NS;
}

//Prior for NRT, with option for tidal_love
double standard_log_prior_D_NRT(double *pos, mcmc_data_interface *interface,void *parameters)
{
	int dim =  interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	double chirp = exp(pos[7]);
	double m1 = calculate_mass1(chirp,pos[8]);
	double m2 = calculate_mass2(chirp,pos[8]);
	double q = m2/m1;//<1
	if(tidal_love){
		if(pos[11]<tidal_s_prior[0] || pos[11]>tidal_s_prior[1]){return a;}
		if(tidal_love_boundary_violation(q,pos[11])){return a;}
		
	}
	else{
		if(pos[11]<tidal1_prior[0] || pos[11]>tidal1_prior[1]){return a;}
		if(pos[12]<tidal2_prior[0] || pos[12]>tidal2_prior[1]){return a;}
	}
	return standard_log_prior_D(pos,interface, parameters);

}
//Prior for PhenomD -- i.e. BBH
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
	if ((pos[9])<spin1_prior[0] || (pos[9])>spin1_prior[1]){return a;}//chi1 
	if ((pos[10])<spin2_prior[0] || (pos[10])>spin2_prior[1]){return a;}//chi1
	//std::cout<<"Made it to end of prior"<<std::endl; 
	return log(aligned_spin_prior(pos[9]))+log(aligned_spin_prior(pos[10])) + log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;

}
//Jacobian factor for \mathcal{M}/\eta to m1,m2
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}

//Try to model uniform spin prior for aligned spins
double aligned_spin_prior(double chi){
	double a=0.0039132 , b= 3.95381;
	return a * exp(-b * abs(chi));	
}

//Tidal love validity boundary in lambda_s / q space
bool tidal_love_boundary_violation(double q,double lambda_s)
{
	//Relation from arXiv:1903.03909v7 fit as rough threshhold for the 
	//validity of the tidal_s-tidal_a-q relation
	//Fit in log(lambda_s) - q space
	if(  q< 1.2321 - .124616*log(lambda_s)){return true;}
	
	return false;
}
