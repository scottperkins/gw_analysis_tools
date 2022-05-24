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

#include <bayesship/bayesshipSampler.h>
#include <bayesship/dataUtilities.h>
#include <bayesship/proposalFunctions.h>
#include <bayesship/utilities.h>

void MCMC_fisher_wrapper_v2(bayesship::positionInfo *pos,   double **output, void *userParameters);



class MCMC_likelihood_wrapper_v2: public bayesship::probabilityFn
{
public:
	bayesship::bayesshipSampler *sampler;
	mcmcVariables *mcmcVar;
	virtual double eval(bayesship::positionInfo *pos, int chainID)
	{
		//return 2;
		//mcmcVariables *mcmcVar = (mcmcVariables *)userParameters;
		//MCMC_user_param *user_param = (MCMC_user_param *)userParameters;
	
		int dimension = sampler->maxDim;
		double ll = 0;
		double *temp_params = new double[dimension];
		//#########################################################################
		gen_params_base<double> gen_params;
		std::string local_gen = MCMC_prep_params_v2(pos->parameters, 
			temp_params,&gen_params, dimension, mcmcVar->mcmc_generation_method,mcmcVar->mcmc_mod_struct, mcmcVar->mcmc_intrinsic,mcmcVar->mcmc_gmst );
		//#########################################################################
		//#########################################################################
	
		//repack_non_parameters(temp_params, &gen_params, 
			//"MCMC_"+mcmc_generation_method, dimension, NULL);
		repack_parameters(temp_params, &gen_params, 
			"MCMC_"+mcmcVar->mcmc_generation_method, dimension, NULL);
		//#########################################################################
		//#########################################################################
		//return 1;
		std::complex<double> **local_data = mcmcVar->mcmc_data;
		double **local_freqs = mcmcVar->mcmc_frequencies;
		double **local_noise = mcmcVar->mcmc_noise;
		double **local_weights = (mcmcVar->user_parameters)->weights;
		int *local_lengths = mcmcVar->mcmc_data_length;
		fftw_outline *local_plans = mcmcVar->mcmc_fftw_plans;
		std::string local_integration_method="SIMPSONS";
		//if(interface->burn_phase && user_param->burn_data){
		if(false){
			local_data = mcmcVar->user_parameters->burn_data;
			local_freqs = mcmcVar->user_parameters->burn_freqs;
			local_noise = mcmcVar->user_parameters->burn_noise;
			local_lengths = mcmcVar->user_parameters->burn_lengths;
			local_plans = mcmcVar->user_parameters->burn_plans;
		}
		if(mcmcVar->user_parameters->GAUSS_QUAD){
			local_integration_method="GAUSSLEG";
		}
		if(mcmcVar->mcmc_intrinsic){
			if(mcmcVar->mcmc_generation_method.find("IMRPhenomD") != std::string::npos){
				if(!mcmcVar->mcmc_save_waveform){
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
						gen_params.theta=0;	
						gen_params.phi=0;	
						gen_params.psi=0;	
						gen_params.phiRef = 1;
						gen_params.f_ref = 10;
						gen_params.incl_angle=0;	
						gen_params.tc =1;
						std::complex<double> *response =
							(std::complex<double> *) malloc(sizeof(std::complex<double>) * local_lengths[i]);
						fourier_detector_response_horizon(local_freqs[i], local_lengths[i], response, mcmcVar->mcmc_detectors[i], local_gen, &gen_params);
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
					fourier_detector_response_horizon(local_freqs[0], local_lengths[0], response, mcmcVar->mcmc_detectors[0], local_gen, &gen_params);
					//std::complex<double> *hc =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//std::complex<double> *hp =
					//	(std::complex<double> *) malloc(sizeof(std::complex<double>) * mcmc_data_length[0]);
					//fourier_waveform(mcmc_frequencies[0], mcmc_data_length[0], hp,hc, local_gen, &gen_params);
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
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
			else if(mcmcVar->mcmc_generation_method.find("IMRPhenomP")!=std::string::npos){
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
					assign_polarizations(mcmcVar->mcmc_generation_method,&wp);
					wp.allocate_memory(local_lengths[0]);
					fourier_waveform(local_freqs[0],local_lengths[0], &wp, local_gen, &gen_params);
					for(int i=0; i < mcmcVar->mcmc_num_detectors; i++){
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
			
			ll =  MCMC_likelihood_extrinsic(mcmcVar->mcmc_save_waveform, 
				&gen_params,local_gen, local_lengths, 
				local_freqs, local_data, local_noise, local_weights, local_integration_method, mcmcVar->user_parameters->log10F,mcmcVar->mcmc_detectors, 
				 mcmcVar->mcmc_num_detectors);
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
};



/*! \brief Takes in an MCMC checkpoint file and continues the chain
 *
 * Obviously, the user must be sure to correctly match the dimension, number of chains, the generation_method, 
 * the prior function, the data, psds, freqs, and the detectors (number and name), and the gps_time to the 
 * previous run, otherwise the behavior of the sampler is undefined.
 *
 * numThreads and pool do not necessarily have to be the same
 */
bayesship::bayesshipSampler *  PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW_v2(
	int dimension,
	int independentSamples,
	int ensembleSize,
	int ensembleN,
	bayesship::positionInfo *initialPosition,
	bayesship::positionInfo **initialEnsemble,
	double swapProb,
	int burnIterations,
	int burnPriorIterations,
	int priorIterations,
	bool writePriorData,
	int max_chunk_size,
	double **priorRanges,
	bayesship::probabilityFn *log_prior,
	int numThreads,
	bool pool,
	int num_detectors,
	std::complex<double> **data,
	double **noise_psd,
	double **frequencies,
	int *data_length,
	double gps_time,
	std::string *detectors,
	MCMC_modification_struct *mod_struct,
	std::string generation_method,
	std::string outputDir,
	std::string outputFileMoniker
	)
{
	int chainN = ensembleSize*ensembleN;
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
	//mcmcVar.mcmc_init_pos = initial_pos;
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
	mcmcVar.maxDim = dimension;



//##########################################################
//##########################################################




	//mcmc_noise = noise_psd;	
	//mcmc_init_pos = initial_pos;
	//mcmc_frequencies = frequencies;
	//mcmc_data = data;
	//mcmc_data_length = data_length;
	//mcmc_detectors = detectors;
	//mcmc_generation_method = generation_method;
	//mcmc_fftw_plans = plans;
	//mcmc_num_detectors = num_detectors;
	//mcmc_gps_time = gps_time;
	//mcmc_gmst = gps_to_GMST_radian(mcmc_gps_time);
	////mcmc_Nmod = mod_struct->ppE_Nmod;
	////mcmc_bppe = mod_struct->bppe;
	//mcmc_mod_struct = mod_struct;
	//mcmc_log_beta = false;
	//mcmc_intrinsic = false;

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
	local_seeding = false;
	//if(!seeding_var){
	//	local_seeding = true;
	//}
	//else{
	//	local_seeding = false;
	//}

	double **seeding_var_ptr = nullptr;
	//PTMCMC_method_specific_prep(generation_method, dimension, seeding_var_ptr, local_seeding);
	PTMCMC_method_specific_prep_v2(generation_method, dimension, seeding_var_ptr, local_seeding,&(mcmcVar.mcmc_intrinsic), mcmcVar.mcmc_mod_struct);
	
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
	

	MCMC_user_param **user_parameters=NULL;
	user_parameters = new MCMC_user_param*[chainN];
	for(int i = 0 ;i<chainN; i++){
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
	//mcmcVar.user_parameters = user_parameters;
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
	
	MCMC_likelihood_wrapper_v2 *ll = new MCMC_likelihood_wrapper_v2();
	ll->mcmcVar = &mcmcVar;
	bayesship::bayesshipSampler *sampler = new bayesship::bayesshipSampler(ll,log_prior );
	ll->sampler = sampler;
	sampler->independentSamples = independentSamples;
	sampler->outputDir = outputDir;
	sampler->outputFileMoniker = outputFileMoniker;
	sampler->burnIterations = burnIterations;
	sampler->burnPriorIterations = burnPriorIterations;
	sampler->priorIterations = priorIterations;
	sampler->writePriorData = writePriorData;
	sampler->threads = numThreads;
	sampler->threadPool = pool;
	sampler->maxDim = dimension;
	sampler->ensembleSize = ensembleSize;
	sampler->ensembleN = ensembleN;
	sampler->priorRanges = priorRanges;
	sampler->initialPosition = initialPosition;

	//Testing
	//sampler->coldOnlyStorage = false;
	sampler->coldOnlyStorage = true;

	mcmcVariables **mcmcVarVec  = new mcmcVariables*[chainN];
	for(int i = 0 ; i<chainN; i++){
		mcmcVarVec[i] = &mcmcVar;
		mcmcVarVec[i]->user_parameters = user_parameters[i];
	}
	sampler->userParameters = (void **) mcmcVarVec;

	//##########################################################
	
	int proposalFnN = 4;
	bayesship::proposal **propArray = new bayesship::proposal*[proposalFnN];
	propArray[0] = new bayesship::gaussianProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler);
	propArray[1] = new bayesship::differentialEvolutionProposal(sampler);
	propArray[2] = new bayesship::KDEProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, sampler, false );
	propArray[3] = new bayesship::fisherProposal(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, MCMC_fisher_wrapper_v2,   sampler->userParameters,  100,sampler);

	//Rough estimate of the temperatures
	double betaTemp[sampler->ensembleSize];
	betaTemp[0] = 1;
	betaTemp[sampler->ensembleSize-1] = 0;
	double deltaBeta = pow((1e2),1./sampler->ensembleSize);
	for(int i = 1 ; i<sampler->ensembleSize-1; i++){
		betaTemp[i] = betaTemp[i-1]/deltaBeta;
	}

	double **propProb = new double*[chainN];	
	for(int i = 0 ; i<chainN; i++){
		int ensemble = i / sampler->ensembleN;
		propProb[i] = new double[proposalFnN];	
		propProb[i][2] = 0.0;

		//propProb[i][0] = 0.05;
		//propProb[i][1] = 0.25;
		//propProb[i][3] = 0.7;
		propProb[i][1] = .7 - 0.45*( betaTemp[ensemble]); //.25 to .7
		propProb[i][3] = 0.25 + .45*( betaTemp[ensemble]); //.7 to .25

		propProb[i][0] = 1. - propProb[i][3] - propProb[i][1]- propProb[i][2];
		//std::cout<<propProb[i][0]<<" "<<propProb[i][1]<<" "<<propProb[i][3]<<std::endl;

	}
	
	//propProb[0] = 0.05;
	//propProb[1] = 0.3;
	//propProb[2] = 0.05;
	//propProb[3] = 0.6;
	//propProb[0] = 0.05;
	//propProb[1] = 0.25;
	//propProb[2] = 0.0;
	//propProb[3] = 0.7;

	//propProb[0] = 0.1;
	//propProb[1] = 0.0;
	//propProb[2] = 0.0;
	//propProb[3] = 0.9;



	bayesship::proposalData *propData = new bayesship::proposalData(chainN,proposalFnN, propArray, (double *)nullptr, propProb);
	
	
	
	
	//int proposalFnN = 4;
	//bayesship::proposalFn propArray[proposalFnN];
	//void *proposalFnVariables[proposalFnN];	
	//propArray[0] = bayesship::gaussianProposal;
	//propArray[1] = bayesship::differentialEvolutionProposal;
	//propArray[2] = bayesship::KDEProposal;
	//propArray[3] = bayesship::FisherProposal;

	////Rough estimate of the temperatures
	//double betaTemp[sampler->ensembleSize];
	//betaTemp[0] = 1;
	//betaTemp[sampler->ensembleSize-1] = 0;
	//double deltaBeta = pow((1e2),1./sampler->ensembleSize);
	//for(int i = 1 ; i<sampler->ensembleSize-1; i++){
	//	betaTemp[i] = betaTemp[i-1]/deltaBeta;
	//}

	//float **propProb = new float*[chainN];	
	//for(int i = 0 ; i<chainN; i++){
	//	int ensemble = i / sampler->ensembleN;
	//	//std::cout<<ensemble<<std::endl;
	//	//std::cout<<betaTemp[ensemble]<<std::endl;
	//	propProb[i] = new float[proposalFnN];	
	//	propProb[i][2] = 0.0;

	//	//propProb[i][0] = 0.05;
	//	//propProb[i][1] = 0.25;
	//	//propProb[i][3] = 0.7;
	//	propProb[i][1] = .7 - 0.45*( betaTemp[ensemble]); //.25 to .7
	//	propProb[i][3] = 0.25 + .45*( betaTemp[ensemble]); //.7 to .25

	//	propProb[i][0] = 1. - propProb[i][3] - propProb[i][1]- propProb[i][2];
	//	//std::cout<<propProb[i][0]<<" "<<propProb[i][1]<<" "<<propProb[i][3]<<std::endl;

	//}
	//
	////propProb[0] = 0.05;
	////propProb[1] = 0.3;
	////propProb[2] = 0.05;
	////propProb[3] = 0.6;
	////propProb[0] = 0.05;
	////propProb[1] = 0.25;
	////propProb[2] = 0.0;
	////propProb[3] = 0.7;

	////propProb[0] = 0.1;
	////propProb[1] = 0.0;
	////propProb[2] = 0.0;
	////propProb[3] = 0.9;


	//bayesship::gaussianProposalVariables *gpv = new bayesship::gaussianProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);

	////bayesship::KDEProposalVariables *kdepv = new bayesship::KDEProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);
	//bayesship::KDEProposalVariables kdepv(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim);

	//bayesship::FisherProposalVariables *fpv = new bayesship::FisherProposalVariables(sampler->ensembleN*sampler->ensembleSize, sampler->maxDim, MCMC_fisher_wrapper_v2,   sampler->userParameters,  100);

	//proposalFnVariables[0] = (void *)gpv;
	//proposalFnVariables[1] = (void *)nullptr;
	//proposalFnVariables[2] = (void *)&kdepv;
	////proposalFnVariables[2] = (void *)nullptr;
	//proposalFnVariables[3] = (void *)fpv;

	//bayesship::proposalFnWriteCheckpoint *writeCheckpointFns = new bayesship::proposalFnWriteCheckpoint[proposalFnN];
	//writeCheckpointFns[0] = bayesship::gaussianProposalWriteCheckpoint;
	//writeCheckpointFns[1] = nullptr;
	//writeCheckpointFns[2] = nullptr;
	//writeCheckpointFns[3] = nullptr;

	//bayesship::proposalFnLoadCheckpoint *loadCheckpointFns = new bayesship::proposalFnLoadCheckpoint[proposalFnN];
	//loadCheckpointFns[0] = bayesship::gaussianProposalLoadCheckpoint;
	//loadCheckpointFns[1] = nullptr;
	//loadCheckpointFns[2] = nullptr;
	//loadCheckpointFns[3] = nullptr;



	//bayesship::proposalFnData *propData = new bayesship::proposalFnData(chainN,proposalFnN, propArray,proposalFnVariables, (float *)nullptr, propProb,writeCheckpointFns,loadCheckpointFns);

	//#######################################################################################
	sampler->proposalFns = propData;
	//##########################################################



	sampler->sample();

	for(int i =0 ; i<chainN ; i++){
		delete [] propProb[i];
	}
	delete [] propProb;


	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated(sampler_output,output, dimension, N_steps, chain_N, 
	//	max_chain_N_thermo_ensemble,initial_pos,seeding_var,ensemble_initial_pos, chain_temps, 
	//	swp_freq, t0, nu,max_chunk_size,chain_distribution_scheme,
	//	log_prior,MCMC_likelihood_wrapper, MCMC_fisher_wrapper,(void**)user_parameters,numThreads, pool, 
	//	//log_prior,MCMC_likelihood_wrapper, NULL,(void**)user_parameters,numThreads, pool, 
	//	show_prog,statistics_filename,
	//	chain_filename, likelihood_log_filename,checkpoint_filename);
	
	

	if(!mcmcVar.mcmc_mod_struct->fisher_weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->fisher_weights;
		}
	}
	if(!mcmcVar.mcmc_mod_struct->weights){
		for(int i = 0 ;i<chainN;i++){
			delete[] user_parameters[i]->weights;
		}
	}
	
	//Deallocate fftw plans
	for (int i =0;i<num_detectors;i++)
		deallocate_FFTW_mem(&plans[i]);
	//#################################################
	delete propData->proposals[0];
	delete propData->proposals[1];
	delete propData->proposals[2];
	delete propData->proposals[3];
	delete propData;
	delete [] propArray;
	for(int i = 0 ; i<num_detectors; i++){
		delete [] burn_data[i];
		delete [] burn_freqs[i];
		delete [] burn_noise[i];
		deallocate_FFTW_mem(&burn_plans[i]);
	}
	delete [] burn_data;
	delete [] burn_lengths;
	delete [] burn_noise;
	delete [] burn_freqs;
	delete [] burn_plans;
	for(int i = 0 ; i<chainN; i++){
		delete user_parameters[i];
	}
	delete [] user_parameters;
	delete [] mcmcVarVec;
	//#################################################
	free(plans);
	//if(local_seeding){ delete [] seeding_var;}
	return sampler;
}
/*! \brief Unpacks MCMC parameters for method specific initiation 
 *
 * Populates seeding vector if non supplied, populates mcmc_Nmod, populates mcmc_log_beta, populates mcmc_intrinsic
 */
void PTMCMC_method_specific_prep_v2(std::string generation_method, int dimension,double **seeding_var, bool local_seeding, bool *intrinsic,MCMC_modification_struct *mod_struct)
{
	int totalmod = (mod_struct->gIMR_Nmod_phi + mod_struct->gIMR_Nmod_sigma + mod_struct->gIMR_Nmod_beta + mod_struct->gIMR_Nmod_alpha  + mod_struct->ppE_Nmod);
	debugger_print(__FILE__,__LINE__,totalmod);
	if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 4)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+4] = 1;
			}
		}
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 6)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, chi1, chi2, tidal1,  tidal2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=10;
			(*seeding_var)[5]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+6] = 1;
			}
		}
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomPv2") != std::string::npos && (dimension - totalmod) == 8)
	{
		std::cout<<"Sampling in parameters: ln chirpmass, eta, a1, a2, tilt1, tilt2, phi1, phi2";
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.5;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.1;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+8] = 1;
			}
		}
		*intrinsic=true;
	} 
	else if(generation_method.find("PhenomD") != std::string::npos && (dimension - totalmod) == 11)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.1;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.5;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			(*seeding_var)[8]=.1;
			(*seeding_var)[9]=.1;
			(*seeding_var)[10]=.5;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+11] = 1;
			}
		}
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomD_NRT") != std::string::npos && (dimension - totalmod) == 13)
	{
		std::cout<<"Sampling in parameters: RA, sin  DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, chi1, chi2, tidal1, tidal2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=.1;
			(*seeding_var)[1]=.1;
			(*seeding_var)[2]=.1;
			(*seeding_var)[3]=.1;
			(*seeding_var)[4]=.5;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=.1;
			(*seeding_var)[7]=.1;
			(*seeding_var)[8]=.1;
			(*seeding_var)[9]=.1;
			(*seeding_var)[10]=.5;
			(*seeding_var)[11]=10;
			(*seeding_var)[12]=10;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+13] = 1;
			}
		}
		*intrinsic=false;
	} 
	else if(generation_method.find("PhenomPv2") != std::string::npos && (dimension - totalmod) == 15)
	{
		std::cout<<"Sampling in parameters: RA, sin DEC, psi, cos iota,phi_ref, tc,  ln DL, ln chirpmass, eta, a1, a2,cos tilt1, cos tilt2, phi1, phi2"<<std::endl;
		for(int i =0; i<totalmod; i++){
			std::cout<<", mod_"<<i;
		}
		std::cout<<std::endl;
		if(local_seeding){
			(*seeding_var) = new double[dimension];
			(*seeding_var)[0]=1.;
			(*seeding_var)[1]=.5;
			(*seeding_var)[2]=1;
			(*seeding_var)[3]=1;
			(*seeding_var)[4]=1;
			(*seeding_var)[5]=.1;
			(*seeding_var)[6]=1;
			(*seeding_var)[7]=1.;
			(*seeding_var)[8]=.3;
			(*seeding_var)[9]=.5;
			(*seeding_var)[10]=.5;
			(*seeding_var)[11]=1.;
			(*seeding_var)[12]=1.;
			(*seeding_var)[13]=1.;
			(*seeding_var)[14]=1.;
			for(int i = 0  ;i <totalmod; i++){
				(*seeding_var)[i+15] = 1;
			}
		}
		*intrinsic=false;
	} 
	else{
		std::cout<<
			"Input parameters not valid, please check that input is compatible with the supported methods - dimension combinations"<<std::endl;
		exit(1);
	}
}



/*! \brief utility to do MCMC specific transformations on the input param vector before passing to the repacking utillity
 *
 * Returns the local generation method to be used in the LL functions
 */
std::string MCMC_prep_params_v2(double *param, double *temp_params, gen_params_base<double> *gen_params, int dimension, std::string generation_method, MCMC_modification_struct *mod_struct, bool intrinsic, double gmst)
{
	if(intrinsic) gen_params->sky_average = true;
	else gen_params->sky_average = false;
	gen_params->f_ref = 20;
	gen_params->shift_time = true;
	gen_params->shift_phase = true;
	//gen_params->shift_time = false;
	//gen_params->shift_phase = false;
	gen_params->gmst = gmst;
	gen_params->equatorial_orientation=false;
	gen_params->horizon_coord=false;

	gen_params->NSflag1 = mod_struct->NSflag1;
	gen_params->NSflag2 = mod_struct->NSflag2;
	//gen_params->NSflag1 = false;
	//gen_params->NSflag2 = false;
	for(int i = 0 ; i <dimension; i++){
		temp_params[i]=param[i];
	}
	int base = dimension;
	if(check_mod(generation_method)){
		//if(generation_method.find("ppE") != std::string::npos ||
		//	generation_method.find("dCS") !=std::string::npos||
		//	generation_method.find("EdGB") != std::string::npos){
		//	gen_params->bppe=mcmc_mod_struct->bppe;
		//	gen_params->Nmod=mcmc_mod_struct->ppE_Nmod;
		//	gen_params->betappe=new double[gen_params->Nmod];
		//}
		if(generation_method.find("ppE") != std::string::npos ||
			check_theory_support(generation_method)){
			gen_params->bppe=mod_struct->bppe;
			gen_params->Nmod=mod_struct->ppE_Nmod;
			gen_params->betappe=new double[gen_params->Nmod];
			base = dimension - mod_struct->ppE_Nmod;
		}
		else if(generation_method.find("gIMR") != std::string::npos){
			gen_params->Nmod_phi=mod_struct->gIMR_Nmod_phi;
			gen_params->phii=mod_struct->gIMR_phii;
			if(gen_params->Nmod_phi !=0){
				gen_params->delta_phi=new double[gen_params->Nmod_phi];
			}
			gen_params->Nmod_sigma=mod_struct->gIMR_Nmod_sigma;
			gen_params->sigmai=mod_struct->gIMR_sigmai;
			if(gen_params->Nmod_sigma !=0){
				gen_params->delta_sigma=new double[gen_params->Nmod_sigma];
			}
			gen_params->Nmod_beta=mod_struct->gIMR_Nmod_beta;
			gen_params->betai=mod_struct->gIMR_betai;
			if(gen_params->Nmod_beta !=0){
				gen_params->delta_beta=new double[gen_params->Nmod_beta];
			}
			gen_params->Nmod_alpha=mod_struct->gIMR_Nmod_alpha;
			gen_params->alphai=mod_struct->gIMR_alphai;
			if(gen_params->Nmod_alpha !=0){
				gen_params->delta_alpha=new double[gen_params->Nmod_alpha];
			}
			base = dimension 
				- mod_struct->gIMR_Nmod_phi 
				- mod_struct->gIMR_Nmod_sigma 
				- mod_struct->gIMR_Nmod_beta 
				- mod_struct->gIMR_Nmod_alpha; 
		}
		if((generation_method.find("dCS")!= std::string::npos ||
			generation_method.find("EdGB")!=std::string::npos)){
			//temp_params[base] = pow(temp_params[base],.25)/(c*1000);
			temp_params[base] = 
				pow_int(temp_params[base]/(c/1000.) , 4);
		}
	}
	return generation_method;
}




void MCMC_fisher_transformations_v2(
	double *param, 
	double **fisher, 
	int dimension,
	std::string generation_method,
	bool intrinsic,
	MCMC_modification_struct *mod_struct
	)
{
	if(!intrinsic){
		fisher[0][0] += 1./(4*M_PI*M_PI);//RA
		fisher[1][1] += 1./4;//sin DEC
		fisher[2][2] += 1./(4*M_PI*M_PI);//psi
		fisher[3][3] += 1./(4);//cos iota
		fisher[4][4] += 1./(4*M_PI*M_PI);//phiref
		fisher[5][5] += 1./(.01);//tc
		fisher[8][8] += 1./.25;//eta
		fisher[9][9] += 1./4;//spin1
		fisher[10][10] += 1./4;//spin2
		if(generation_method.find("PhenomPv2") != std::string::npos){
			fisher[11][11] += 1./4;//cos theta1
			fisher[12][12] += 1./4;//cos theta2
			fisher[13][13] += 1./(4*M_PI*M_PI);//phi1
			fisher[14][14] += 1./(4*M_PI*M_PI);//phi2
		}
	}
	else{
		if(generation_method.find("PhenomPv2") != std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4);//spin1
			fisher[3][3] =1./(4);//spin2
			fisher[4][4] =1./(4);//cos theta1
			fisher[5][5] =1./(4);//cos theta2
			fisher[6][6] =1./(4*M_PI*M_PI) ;//phi1
			fisher[7][7] =1./(4*M_PI*M_PI) ;//phi2
		}
		else if (generation_method.find("PhenomD")!=std::string::npos){
			fisher[1][1] =1./(.25) ;//eta
			fisher[2][2] =1./(4) ;//spin1
			fisher[3][3] =1./(4) ;//spin2
	
		}
	}

	if(generation_method.find("dCS") != std::string::npos ||
		generation_method.find("EdGB") != std::string::npos){
		int base = dimension - mod_struct->ppE_Nmod;
		for(int i = 0 ; i<dimension; i++){
			//Transform to root alpha from alpha^2
			double factor = 4* pow(param[base], 3./4.);
			//Transform to km from sec
			factor *= 1000 / c ;
			fisher[base][i] *= factor ;
			fisher[i][base] *= factor;
		}
	}
	return;

}



void MCMC_fisher_wrapper_v2(bayesship::positionInfo *pos,   double **output, void *userParameters)
{
	mcmcVariables *mcmcVar= (mcmcVariables *)userParameters;

	int dimension = mcmcVar->maxDim;
	double *temp_params = new double[dimension];
	double param[dimension];
	for(int i = 0 ; i<dimension; i++){	
		param[i] = pos->parameters[i];
	}
	//#########################################################################
	gen_params_base<double> params;
	std::string local_gen = MCMC_prep_params(param, 
		temp_params,&params, dimension, mcmcVar->mcmc_generation_method,mcmcVar->mcmc_mod_struct);
	//#########################################################################
	//#########################################################################
	repack_parameters(temp_params, &params, 
		"MCMC_"+mcmcVar->mcmc_generation_method, dimension, NULL);
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
	double **local_freq=mcmcVar->mcmc_frequencies;
	double **local_noise=mcmcVar->mcmc_noise;
	double **local_weights=mcmcVar->user_parameters->weights;
	int *local_lengths= mcmcVar->mcmc_data_length;
	std::string local_integration_method = "SIMPSONS";
	bool local_log10F = mcmcVar->user_parameters->fisher_log10F;
	if(mcmcVar->user_parameters->fisher_freq){
		local_freq = mcmcVar->user_parameters->fisher_freq;
	}
	if(mcmcVar->user_parameters->fisher_PSD){
		local_noise = mcmcVar->user_parameters->fisher_PSD;
	}
	if(mcmcVar->user_parameters->fisher_length){
		local_lengths = mcmcVar->user_parameters->fisher_length;
	}
	if(mcmcVar->user_parameters->fisher_weights){
		local_weights = mcmcVar->user_parameters->fisher_weights;
	}
	if(mcmcVar->user_parameters->fisher_GAUSS_QUAD){
		local_integration_method = "GAUSSLEG";
	}
	double **temp_out = allocate_2D_array(dimension,dimension);
	for(int i =0 ; i <mcmcVar->mcmc_num_detectors; i++){
		
		//Use AD 
		if(mcmcVar->user_parameters->fisher_AD)
		{	
			std::unique_lock<std::mutex> lock{*(mcmcVar->user_parameters->mFish)};
			//fisher_autodiff(mcmcVar->mcmc_frequencies[i], mcmcVar->mcmc_data_length[i],
			//	"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
			//	(gen_params *)(&params),  "SIMPSONS",(double *)NULL,false,mcmcVar->mcmc_noise[i]);
			fisher_autodiff(local_freq[i], local_lengths[i],
				"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
				(gen_params *)(&params),  local_integration_method,local_weights[i],true,local_noise[i]);
		}
		else{
			fisher_numerical(local_freq[i], local_lengths[i],
				"MCMC_"+mcmcVar->mcmc_generation_method, mcmcVar->mcmc_detectors[i],mcmcVar->mcmc_detectors[0],temp_out,dimension, 
				&params, 4, NULL, NULL, local_noise[i]);

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
	//if(mcmcVar->mcmc_generation_method.find("Pv2") && !mcmcVar->mcmc_intrinsic){
	
	MCMC_fisher_transformations_v2(temp_params, output,dimension,local_gen,mcmcVar->mcmc_intrinsic,
		mcmcVar->mcmc_mod_struct);
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
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_phi !=0){
				delete [] params.delta_phi;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_sigma !=0){
				delete [] params.delta_sigma;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_beta !=0){
				delete [] params.delta_beta;
			}
			if(mcmcVar->mcmc_mod_struct ->gIMR_Nmod_alpha !=0){
				delete [] params.delta_alpha;
			}

		}
	}


}
