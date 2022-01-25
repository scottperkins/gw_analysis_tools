#include "mcmc_gw_extended.h"


#include "mcmc_gw.h"
#include "waveform_generator.h"
#include "util.h"
#include "io_util.h"
#include "detector_util.h"
#include "ppE_utilities.h"
#include "waveform_util.h"
#include "ortho_basis.h"
#include "fisher.h"

#include <ptrjmcmc/PtrjmcmcSampler.h>

/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
void PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(mcmc_sampler_output *sampler_output,
	double **output,
	int dimension,
	int N_steps,
	int chain_N,
	int max_chain_N_thermo_ensemble,
	double *initial_pos,
	double *seeding_var,
	double **ensemble_initial_pos,
	double *chain_temps,
	int swp_freq,
	int t0,
	int nu,
	int max_chunk_size,
	std::string chain_distribution_scheme,
	double(*log_prior)(double *param, mcmc_data_interface *interface,void *parameters),
	int numThreads,
	bool pool,
	bool show_prog,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string statistics_filename,
	std::string chain_filename,
	std::string likelihood_log_filename,
	std::string checkpoint_filename
	)
{
	//Create fftw plan for each detector (length of data stream may be different)
	fftw_outline *plans= (fftw_outline *)malloc(sizeof(fftw_outline)*num_detectors);
	for (int i =0;i<num_detectors;i++)
	{	
		allocate_FFTW_mem_forward(&plans[i] , data_length[i]);
	}

	std::mutex fisher_mutex;

//##########################################################
//##########################################################
	mcmcVariables mcmcVar ;
	mcmcVar.mcmc_noise = noise_psd;
	mcmcVar.mcmc_init_pos = initial_pos;
	mcmcVar.mcmc_frequencies = frequencies;
	mcmcVar.mcmc_data = data;
	mcmcVar.mcmc_data_length = data_length;
	mcmcVar.mcmc_detectors = detectors;
	mcmcVar.mcmc_generation_method = generation_method;
	mcmcVar.mcmc_fftw_plans = plans;
	mcmcVar.mcmc_num_detectors = num_detectors;
	mcmcVar.mcmc_gps_time = gps_time;
	mcmcVar.mcmc_gmst = gps_to_GMST_radian(gps_time);
	mcmcVar.mcmc_mod_struct = mod_struct;
	mcmcVar.mcmc_save_waveform = true;



//##########################################################
//##########################################################




	mcmc_noise = noise_psd;	
	mcmc_init_pos = initial_pos;
	mcmc_frequencies = frequencies;
	mcmc_data = data;
	mcmc_data_length = data_length;
	mcmc_detectors = detectors;
	mcmc_generation_method = generation_method;
	mcmc_fftw_plans = plans;
	mcmc_num_detectors = num_detectors;
	mcmc_gps_time = gps_time;
	mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	//mcmc_Nmod = mod_struct->ppE_Nmod;
	//mcmc_bppe = mod_struct->bppe;
	mcmc_mod_struct = mod_struct;
	mcmc_log_beta = false;
	mcmc_intrinsic = false;

	//To save time, intrinsic waveforms can be saved between detectors, if the 
	//frequencies are all the same
	mcmc_save_waveform = true;
	for(int i =1 ;i<mcmcVar.mcmc_num_detectors; i++){
		if( mcmcVar.mcmc_data_length[i] != mcmcVar.mcmc_data_length[0] ||
			mcmcVar.mcmc_frequencies[i][0]!= mcmcVar.mcmc_frequencies[0][0] ||
			mcmcVar.mcmc_frequencies[i][mcmcVar.mcmc_data_length[i]-1] 
				!= mcmcVar.mcmc_frequencies[0][mcmcVar.mcmc_data_length[0]-1]){
			mcmcVar.mcmc_save_waveform= false;
		}
			
	}

	bool local_seeding ;
	if(!seeding_var){
		local_seeding = true;
	}
	else{
		local_seeding = false;
	}

	double **seeding_var_ptr = &seeding_var;
	PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);
	
	//######################################################
	int T = (int)(1./(mcmcVar.mcmc_frequencies[0][1]-mcmcVar.mcmc_frequencies[0][0]));
	debugger_print(__FILE__,__LINE__,T);
	int burn_factor = T/4; //Take all sources to 4 seconds
	debugger_print(__FILE__,__LINE__,burn_factor);
	std::complex<double> **burn_data = new std::complex<double>*[mcmcVar.mcmc_num_detectors];
	double **burn_freqs = new double*[mcmcVar.mcmc_num_detectors];
	double **burn_noise = new double*[mcmcVar.mcmc_num_detectors];
	int *burn_lengths = new int[mcmcVar.mcmc_num_detectors];
	fftw_outline *burn_plans= new fftw_outline[mcmcVar.mcmc_num_detectors];
	for(int j = 0; j<mcmcVar.mcmc_num_detectors; j++){
		burn_lengths[j] = mcmcVar.mcmc_data_length[j]/burn_factor;
		burn_data[j]= new std::complex<double>[burn_lengths[j]];
		burn_freqs[j]= new double[burn_lengths[j]];
		burn_noise[j]= new double[burn_lengths[j]];
		allocate_FFTW_mem_forward(&burn_plans[j], burn_lengths[j]);
		int ct = 0;
		for( int k = 0 ; k<mcmcVar.mcmc_data_length[j]; k++){
			if(k%burn_factor==0 && ct<burn_lengths[j]){
				burn_data[j][ct] = mcmcVar.mcmc_data[j][k];
				burn_freqs[j][ct] = mcmcVar.mcmc_frequencies[j][k];
				burn_noise[j][ct] = mcmcVar.mcmc_noise[j][k];
				ct++;
			}
		}
	}
	

	//######################################################
	//######################################################
	//Fishers sometimes need AD, but that's slow and single 
	//threaded -- use GAUSS quad with precomputed weights 
	//and abscissa 
	//double *fish_freqs = NULL;
	//double *fish_weights = NULL;
	//double **fish_psd = NULL; //Needs to be interpolated from data given
	//int fish_length = 100;
	//double flow = mcmc_frequencies[0][0];
	//double fhigh = mcmc_frequencies[0][mcmc_data_length[0]-1];
	//if(mod_struct->GAUSS_QUAD){
	//	fish_freqs = new double[fish_length];	
	//	fish_weights = new double[fish_length];	
	//	fish_psd = new double*[mcmc_num_detectors];	
	//	
	//	gauleg(log10(flow), log10(fhigh), fish_freqs, fish_weights, fish_length);
	//	for(int i = 0 ; i<fish_length; i++){
	//		fish_freqs[i]=pow(10,fish_freqs[i]);
	//	}
	//	for(int i = 0 ; i<mcmc_num_detectors; i++){
	//		fish_psd[i] = new double[fish_length];
	//		gsl_interp_accel *accel = gsl_interp_accel_alloc();
	//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, mcmc_data_length[i]);
	//		gsl_spline_init(spline, mcmc_frequencies[i], mcmc_noise[i],mcmc_data_length[i]);
	//		for(int j = 0 ; j<fish_length; j++){
	//			fish_psd[i][j] = gsl_spline_eval(spline, fish_freqs[j],accel);
	//		}
	//		gsl_spline_free(spline);
	//	}
	//}
	//######################################################
	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chain_N];
	for(int i = 0 ;i<chain_N; i++){
		user_parameters[i] = new MCMC_user_param;
		
		user_parameters[i]->burn_data = burn_data;
		user_parameters[i]->burn_freqs = burn_freqs;
		user_parameters[i]->burn_noise = burn_noise;
		user_parameters[i]->burn_lengths = burn_lengths;
		user_parameters[i]->burn_plans=burn_plans;

		user_parameters[i]->mFish= &fisher_mutex;
		user_parameters[i]->GAUSS_QUAD= mod_struct->GAUSS_QUAD;
		user_parameters[i]->log10F = mod_struct->log10F;

		if(mod_struct->weights){
			user_parameters[i]->weights = mod_struct->weights;			
		}
		else{
			user_parameters[i]->weights = new double*[num_detectors];			
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->weights[j]=NULL;
			}
		}

		user_parameters[i]->fisher_GAUSS_QUAD = mod_struct->fisher_GAUSS_QUAD;
		user_parameters[i]->fisher_log10F = mod_struct->fisher_log10F;
		user_parameters[i]->fisher_freq= mod_struct->fisher_freq;
		if(mod_struct->fisher_weights){
			user_parameters[i]->fisher_weights= mod_struct->fisher_weights;
		}
		else{
			user_parameters[i]->fisher_weights = new double*[num_detectors];
			for(int j = 0 ; j<num_detectors; j++){
				user_parameters[i]->fisher_weights[j]=NULL;
			}
		}	
		user_parameters[i]->fisher_PSD= mod_struct->fisher_PSD;
		user_parameters[i]->fisher_length= mod_struct->fisher_length;


		user_parameters[i]->mod_struct = mod_struct;

		
		//user_parameters[i]->burn_freqs = mcmc_frequencies;
		//user_parameters[i]->burn_data = mcmc_data;
		//user_parameters[i]->burn_noise = mcmc_noise;
		//user_parameters[i]->burn_lengths = mcmc_data_length;
	}
	//######################################################
	//###########################################################
	//double trial[11] = {1,.9,.2,.9,.2,3,7,4,.24,0,0};
	//mcmc_data_interface i;
	//i.min_dim = 11; 
	//i.max_dim = 11; 
	//i.chain_id  = 0; 
	//i.chain_number  = 11; 
	//double ll = MCMC_likelihood_wrapper(trial, &i ,user_parameters[0]);
	//debugger_print(__FILE__,__LINE__,ll);
	//###########################################################


	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
		max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
		swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
		log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void**)user_parameters,numThreads, pool, 
		//log_prior,MCMC_likelihood_wrapper, NULL,(void**)user_parameters,numThreads, pool, 
		show_prog,statistics_filename,
		chain_filename, likelihood_log_filename,checkpoint_filename);
	if(!mod_struct->fisher_weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mod_struct->weights){
		for(int i = 0 ;i<chain_N;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	//if(mod_struct->GAUSS_QUAD){
	//	delete [] fish_freqs;
	//	delete [] fish_weights;
	//	for(int i = 0 ;i<mcmc_num_detectors; i++){
	//		delete [] fish_psd[i];
	//	}
	//	delete [] fish_psd;
	//}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chain_N; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	//#################################################
	free(plans);
	if(local_seeding){ delete [] seeding_var;}
}



double MCMC_likelihood_wrapper_v2(double *param, mcmc_data_interface *interface ,void *parameters)
{
	//return 2;
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;

	int dimension = interface->max_dim;
	double ll = 0;
	double *temp_params = new double[dimension];
	//#########################################################################
	gen_params_base<double> gen_params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&gen_params, dimension, mcmc_generation_method,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################

	//repack_non_parameters(temp_params, &gen_params, 
		//"MCMC_"+mcmc_generation_method, dimension, NULL);
	repack_parameters(temp_params, &gen_params, 
		"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//return 1;
	std::complex<double> **local_data = mcmc_data;
	double **local_freqs = mcmc_frequencies;
	double **local_noise = mcmc_noise;
	double **local_weights = user_param->weights;
	int *local_lengths = mcmc_data_length;
	fftw_outline *local_plans = mcmc_fftw_plans;
	std::string local_integration_method="SIMPSONS";
	//if(interface->burn_phase && user_param->burn_data){
	if(false){
		local_data = user_param->burn_data;
		local_freqs = user_param->burn_freqs;
		local_noise = user_param->burn_noise;
		local_lengths = user_param->burn_lengths;
		local_plans = user_param->burn_plans;
	}
	if(user_param->GAUSS_QUAD){
		local_integration_method="GAUSSLEG";
	}
	if(mcmc_intrinsic){
		if(mcmc_generation_method.find("IMRPhenomD") != std::string::npos){
			if(!mcmc_save_waveform){
				for(int i=0; i < mcmc_num_detectors; i++){
					gen_params.theta=0;	
					gen_params.phi=0;	
					gen_params.psi=0;	
					gen_params.phiRef = 1;
					gen_params.f_ref = 10;
					gen_params.incl_angle=0;	
					gen_params.tc =1;
					std::complex<double> *response =
						(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
					fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmc_detectors[i], local_gen, &gen_params);
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
					//ll += maximized_Log_Likelihood(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		(size_t) mcmc_data_length[i],
					//		&gen_params,
					//		mcmc_detectors[i],
					//		local_gen,
					//		&mcmc_fftw_plans[i]
					//		);
					free(response);
				}
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 10;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				std::complex<double> *response =
					(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[0]);
				fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmc_detectors[0], local_gen, &gen_params);
				//std::complex<double> *hc =
				//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
				//std::complex<double> *hp =
				//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
				//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_aligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							response,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
					//ll += maximized_Log_Likelihood(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		(size_t) mcmc_data_length[i],
					//		&gen_params,
					//		mcmc_detectors[i],
					//		local_gen,
					//		&mcmc_fftw_plans[i]
					//		);
					//ll += maximized_Log_Likelihood_unaligned_spin_internal(mcmc_data[i], 
					//		mcmc_noise[i],
					//		mcmc_frequencies[i],
					//		hp,
					//		hc,
					//		(size_t) mcmc_data_length[i],
					//		&mcmc_fftw_plans[i]
					//		);
				}
				free(response);
				//free(hp);
				//free(hc);

			}

		}
		else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
			//if(!mcmc_save_waveform){
			if(false){
			}
			else{
				gen_params.theta=0;	
				gen_params.phi=0;	
				gen_params.psi=0;	
				gen_params.phiRef = 1;
				gen_params.f_ref = 20;
				gen_params.incl_angle=0;	
				gen_params.tc =1;
				waveform_polarizations<double> wp;
				assign_polarizations(mcmc_generation_method,&wp);
				wp.allocate_memory(local_lengths[0]);
				fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
				for(int i=0; i < mcmc_num_detectors; i++){
					ll += maximized_Log_Likelihood_unaligned_spin_internal(local_data[i], 
							local_noise[i],
							local_freqs[i],
							wp.hplus,
							wp.hcross,
							(size_t) local_lengths[i],
							&local_plans[i]
							);
				}
				wp.deallocate_memory();
			}

		}
	}
	else{
		double RA = gen_params.RA;
		double DEC = gen_params.DEC;
		double PSI = gen_params.psi;
		//if(mcmc_generation_method.find("IMRPhenomD") != std::string:npos){
		
		ll =  MCMC_likelihood_extrinsic(mcmc_save_waveform, 
			&gen_params,local_gen, local_lengths, 
			local_freqs, local_data, local_noise, local_weights, local_integration_method, user_param->log10F,mcmc_detectors, 
			 mcmc_num_detectors);
		//ll=2;

		//ll = Log_Likelihood(mcmc_data[0], 
		//		mcmc_noise[0],
		//		mcmc_frequencies[0],
		//		mcmc_data_length[0],
		//		&gen_params,
		//		mcmc_detectors[0],
		//		local_gen,
		//		&mcmc_fftw_plans[0]
		//		);

		//}
		//else if(mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){

		//}
	}
	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		//if( local_gen.find("ppE") != std::string::npos ||
		//	local_gen.find("dCS") != std::string::npos ||
		//	local_gen.find("EdGB") != std::string::npos){
		//	delete [] gen_params.betappe;
		//}
		if( local_gen.find("ppE") != std::string::npos ||
			check_theory_support(local_gen)){
			delete [] gen_params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] gen_params.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] gen_params.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] gen_params.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] gen_params.delta_alpha;
			}

		}
	}
	return ll;

}

void MCMC_fisher_wrapper_v2(double *param,  double **output, mcmc_data_interface *interface,void *parameters)
{
	MCMC_user_param *user_param = (MCMC_user_param *)parameters;
	int dimension = interface->max_dim;
	double *temp_params = new double[dimension];
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&params, dimension, mcmc_generation_method,mcmc_mod_struct);
	//#########################################################################
	//#########################################################################
	repack_parameters(temp_params, &params, 
		"MCMC_"+mcmc_generation_method, dimension, NULL);
	//std::cout<<temp_params[11]<<" "<<temp_params[12]<<std::endl;
	//repack_parameters(mcmc_init_pos, &params, 
	//	"MCMC_"+mcmc_generation_method, dimension, NULL);
	//#########################################################################
	//#########################################################################
	//std::cout<<"INCL angle fisher: "<<params.incl_angle<<std::endl;
	for(int j =0; j<dimension; j++){
		for(int k =0; k<dimension; k++)
		{
			output[j][k] =0;
		}
	} 
	double **local_freq=mcmc_frequencies;
	double **local_noise=mcmc_noise;
	double **local_weights=user_param->weights;
	int *local_lengths= mcmc_data_length;
	std::string local_integration_method = "SIMPSONS";
	bool local_log10F = user_param->fisher_log10F;
	if(user_param->fisher_freq){
		local_freq = user_param->fisher_freq;
	}
	if(user_param->fisher_PSD){
		local_noise = user_param->fisher_PSD;
	}
	if(user_param->fisher_length){
		local_lengths = user_param->fisher_length;
	}
	if(user_param->fisher_weights){
		local_weights = user_param->fisher_weights;
	}
	if(user_param->fisher_GAUSS_QUAD){
		local_integration_method = "GAUSSLEG";
	}
	double **temp_out = allocate_2D_array(dimension,dimension);
	for(int i =0 ; i <mcmc_num_detectors; i++){
		
		//Use AD 
		if(user_param->fisher_AD)
		{	
			std::unique_lock<std::mutex> lock{*(user_param->mFish)};
			//fisher_autodiff(mcmc_frequencies[i], mcmc_data_length[i],
			//	"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmc_noise[i]);
			fisher_autodiff(local_freq[i], local_lengths[i],
				"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
				(gen_params *)(&params),  local_integration_method,local_weights[i],true,local_noise[i]);
		}
		else{
			fisher_numerical(local_freq[i], local_lengths[i],
				"MCMC_"+mcmc_generation_method, mcmc_detectors[i],mcmc_detectors[0],temp_out,dimension, 
				&params, mcmc_deriv_order, NULL, NULL, local_noise[i]);

		}
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++)
			{
				output[j][k] +=temp_out[j][k];
				//if(std::isnan(output[j][k]))
				//{
				//      std::cout<<j<<" "<<k<<" "<<temp_out[j][k]<<std::endl;
				//}
			}
		} 
	}
	//Add prior information to fisher
	//if(mcmc_generation_method.find("Pv2") && !mcmc_intrinsic){
	
	MCMC_fisher_transformations(temp_params, output,dimension,local_gen,mcmc_intrinsic,
		interface,mcmc_mod_struct, parameters);
	deallocate_2D_array(temp_out, dimension,dimension);
	//////////////////////////////////////////////
	//if(!interface->burn_phase)
	//{
	//	debugger_print(__FILE__,__LINE__,"Fisher MCMC");
	//	double **cov = allocate_2D_array( dimension,dimension);
	//	gsl_cholesky_matrix_invert(output, cov, dimension);
	//	for(int i = 0 ; i<dimension; i++){
	//		std::cout<<sqrt(cov[i][i])<<" ";	
	//		//for(int j = 0 ; j<dimension; j++){
	//		//	std::cout<<cov[i][j]<<" ";	
	//		//}
	//		//std::cout<<std::endl;	
	//		
	//	}
	//	std::cout<<std::endl;	
	//	deallocate_2D_array(cov, dimension,dimension);
	//}
	//////////////////////////////////////////////

	//Cleanup
	delete [] temp_params;
	if(check_mod(local_gen)){
		//if(local_gen.find("ppE") != std::string::npos ||
		//	local_gen.find("dCS")!=std::string::npos||
		//	local_gen.find("EdGB")!=std::string::npos){
		//	delete [] params.betappe;
		//}
		if(local_gen.find("ppE") != std::string::npos ||
			check_theory_support(local_gen)){
			delete [] params.betappe;
		}
		else if( local_gen.find("gIMR") != std::string::npos){
			if(mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}
