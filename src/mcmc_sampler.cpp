#include "mcmc_sampler.h"
#include "autocorrelation.h"
#include "util.h"
#include "mcmc_sampler_internals.h"
#include "threadPool.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include <time.h>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <unistd.h>

#ifndef _OPENMP
#define omp ignore
#endif


/*!\file 
 * Source file for the sampler foundation
 *
 * Source file for generic MCMC sampler. Sub routines that are 
 * application agnostic are housed in mcmc_sampler_internals
 */

//Random number variables
const gsl_rng_type *T;
gsl_rng * r;
sampler *samplerptr;

/*! \brief Class to facilitate the comparing of chains for priority
 *
 * 3 levels of priority: 0 (high) 1 (default) 2 (low)
 */
class Comparator
{
public:
	bool operator()(int i, int j)
	{
		return samplerptr->priority[i]>samplerptr->priority[j];	
	}
};
class Comparatorswap
{
public:
	bool operator()(int i, int j)
	{
		return false;	
	}
};
class ThreadPool

{
public:
	//using Task = std::function<void()>;
	
	
	explicit ThreadPool(std::size_t numThreads)
	{
		start(numThreads);
	}

	~ThreadPool()
	{
		stop();
	}


	void enqueue(int i)
	{
		{
			std::unique_lock<std::mutex> lock{mEventMutex};
			mTasks.emplace(std::move(i));
		}
		mEventVar.notify_one();
	}

	void enqueue_swap(int i)
	{
		{
			std::unique_lock<std::mutex> lock{mEventMutexSWP};
			mSwaps.emplace(std::move(i));
		}
		mEventVarSWP.notify_one();
	}
	
	void public_stop()
	{
		stop();
	}
private:
	std::vector<std::thread> mThreads;
	
	std::condition_variable mEventVar;

	std::mutex mEventMutex;

	bool mStopping = false;
	
	int numSwpThreads= 1;
		
	std::condition_variable mEventVarSWP;

	std::mutex mEventMutexSWP;

	//std::queue<Task> mTasks;
	//std::queue<int> mTasks;
	//std::queue<int> mSwaps;
	std::priority_queue<int,std::vector<int>,Comparator> mTasks;
	//std::priority_queue<int,std::vector<int>,Comparator> mSwaps;
	std::priority_queue<int,std::vector<int>,Comparatorswap> mSwaps;

	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads-numSwpThreads; i++)
		{
			mThreads.emplace_back([=]{
				while(true)
				{
					int j;
					{
						std::unique_lock<std::mutex> lock{mEventMutex};

						//mEventVar.wait(lock,[=]{return mStopping || !mTasks.empty(); });
						mEventVar.wait(lock,[=]{return mStopping || !mTasks.empty(); });
						
						if (mStopping && mTasks.empty())
							break;	
						//j = std::move(mTasks.front());
						j = std::move(mTasks.top());
						mTasks.pop();
						//std::cout<<mTasks.empty();
					}
					mcmc_step_threaded(j);
					
				}
			});
		}

		//Swapping thread
		for(auto i =0u; i<numSwpThreads; i++)
		{
			mThreads.emplace_back([=]{
				while(true)
				{
					int j, k;
					{
						std::unique_lock<std::mutex> lock{mEventMutexSWP};

						//mEventVarSWP.wait(lock,[=]{return mStopping || !mSwaps.empty(); });
						mEventVarSWP.wait(lock,[=]{return mStopping || !(mSwaps.size()<2); });
						
						if (mStopping && mSwaps.size()<2)
							break;	
						//j = std::move(mSwaps.front());
						j = std::move(mSwaps.top());
						mSwaps.pop();
						//k = std::move(mSwaps.front());
						k = std::move(mSwaps.top());
						mSwaps.pop();
					}
					mcmc_swap_threaded(j,k);
					
				}
			});
		}
	}
	void stop() noexcept
	{
		std::cout<<std::endl;
		std::cout<<"Stop initiated -- waiting for threads to finish"<<std::endl;
		{
			std::unique_lock<std::mutex> lock{mEventMutex};
			//std::unique_lock<std::mutex> lock{mEventMutexSWP};
			mStopping = true;
		}
		
		mEventVar.notify_all();
		mEventVarSWP.notify_all();
		
		for(auto &thread: mThreads)
			thread.join();
	}
};
ThreadPool *poolptr;
/*!\brief Generic reversable jump sampler, where the likelihood, prior, and reversable jump proposal are parameters supplied by the user
 *
 * Currently, no dynamic PT option, as it would be too many free parameters for the sampler to converge to a reasonable temperature distribution in a reasonable amount of time. Best use case, use the PTMCMC_MH_dyanmic_PT for the ``base'' dimension space, and use that temperature ladder.
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
 *
 * 3 modes to use - 
 *
 * single threaded (numThreads = 1) runs single threaded
 *
 * multi-threaded ``deterministic'' (numThreads>1 ; pool = false) progresses each chain in parallel for swp_freq steps, then waits for all threads to complete before swapping temperatures in sequenctial order (j, j+1) then (j+1, j+2) etc (sequenctially)
 *
 * multi-threaded ``stochastic'' (numThreads>2 ; pool = true) progresses each chain in parallel by queueing each temperature and evaluating them in the order they were submitted. Once finished, the threads are queued to swap, where they swapped in the order they are submitted. This means the chains are swapped randomly, and the chains do NOT finish at the same time. The sampler runs until the the 0th chain reaches the step number
 *
 * Note on limits: In the prior function, if a set of parameters should be disallowed, return -std::numeric_limits<double>::infinity()  -- (this is in the <limits> file in std)
 *
 * The parameter array uses the dimensions [0,min_dim] always, and [min_dim, max_dim] in RJPTMCMC fashion
 *
 * Format for the auto_corr file (compatable with csv, dat, txt extensions): each row is a dimension of the cold chain, with the first row being the lengths used for the auto-corr calculation:
 *
 * lengths: length1 , length2 ...
 *
 * dim1: length1 , length2 ...
 *
 * .
 *
 * .
 *
 * .
 *
 *
 * Format for the chain file (compatable with csv, dat, txt extensions): each row is a step, each column a dimension:
 *
 * Step1: dim1 , dim2 , ..., max_dim, param_status1, param_status2, ...
 *
 * Step2: dim1 , dim2 , ..., max_dim, param_status1, param_status2, ...
 *
 * .
 *
 * .
 *
 * .
 *
 * Statistics_filename : should be txt extension
 *
 * checkpoint_file : This file saves the final position of all the chains, as well as other metadata, and can be loaded by the function <FUNCTION> to continue the chain from the point it left off. Not meant to be read by humans, the data order is custom to this software library. An empty string ("") means no checkpoint will be saved. For developers, the contents are:
 *
 * dimension, # of chains
 *
 * temps of chains
 *
 * Stepping widths of all chains
 *
 * Final position of all chains
 */
void RJPTMCMC_MH_internal(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int max_dimension, 	/**< dimension of the parameter space being explored*/
	int min_dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	int initial_dim, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[max_dimension] -- initial seeding of zero corresponds to the dimension turned off initially*/
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	std::function<double(double*,int,int)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int,int)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int,double**,int)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	std::function<double(double*,int,int)> RJ_proposal,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	//To Do:
	//	Retrofit the sampler_internal functions to accept dim-status array -- for regular PTMCMC, just submit with all 1's
	//		Includes all the steps
	//		Includes statistics file, percentage of time spent on/off, etc
	//	retrofit PTMCMC/PTMCMC_dynamic_PT, thread pool functions
	//	Maybe get rid of sampler->output in favor of sampler->output with dimension [2xmax_dim][N_steps], for parameter status always
		
}

/*! \brief Dyanmically tunes an MCMC for optimal spacing. step width, and chain number
 *
 * Based on arXiv:1501.05823v3
 *
 * Currently, Chain number is fixed
 *
 * max_chain_N_thermo_ensemble sets the maximium number of chains to use to in successively hotter chains to cover the likelihood surface while targeting an optimal swap acceptance target_swp_acc. 
 *
 * max_chain_N determines the total number of chains to run once thermodynamic equilibrium has been reached. This results in chains being added after the initial PT dynamics have finished according to chain_distribution_scheme.
 *
 * If no preference, set max_chain_N_thermo_ensemble = max_chain_N = numThreads = (number of cores (number of threads if hyperthreaded))-- this will most likely be the most optimal configuration. If the number of cores on the system is low, you may want to use n*numThreads for some integer n instead, depending on the system.
 *
 * chain_distribution_scheme:
 *
 * "cold": All chains are added at T=1 (untempered)
 *
 * "refine": Chains are added between the optimal temps geometrically -- this may be a good option as it will be a good approximation of the ideal distribution of chains, while keeping the initial dynamical time low 
 *
 * "double": Chains are added in order of rising temperature that mimic the distribution achieved by the earier PT dynamics
 */
void PTMCMC_MH_dynamic_PT_alloc_internal(double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	std::function<double(double*,int,int)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int,int)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int,double**,int)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler samplerobj;
	samplerptr = &samplerobj;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = swp_freq;
	//For PT dynamics
	samplerptr->N_steps = N_steps;

	samplerptr->dimension = dimension;
	samplerptr->num_threads = numThreads;
	samplerptr->output =output;

	//Start out with geometrically spaced chain
	//samplerptr->chain_temps =new double [max_chain_N_thermo_ensemble];
	samplerptr->chain_temps = chain_temps;
	samplerptr->chain_temps[0] = 1.;
	samplerptr->chain_N = max_chain_N_thermo_ensemble;
	double c = 1.5;
	for(int i = 1; i<samplerptr->chain_N-1; i++){
		samplerptr->chain_temps[i] = c * samplerptr->chain_temps[i-1];
	}
	//Set last chain essentially at infinity
	samplerptr->chain_temps[samplerptr->chain_N -1] = 1.e14;

	//NOTE: currently, update the history every step until length is 
	//reached, then the history is updated every 20th step, always only
	//keeping history of length history_length (overwrites the list as 
	//it walks forward when it reaches the end)
	//samplerptr->history_length = 500;
	samplerptr->history_length = 500;
	samplerptr->history_update = 5;
	//Number of steps to take with the fisher before updating the fisher 
	//to a new value 
	//NOTE: if this is too low, detailed balance isn't maintained without 
	//accounting for the changing fisher (doesn't cancel in MH ratio)
	//but if the number is high enough, detailed balance is approximately 
	//kept without calculating second fisher
	samplerptr->fisher_update_number = 200;

	//samplerptr->output = output;
	//During chain allocation, pooling isn't used
	samplerptr->pool = false;
	samplerptr->numThreads = numThreads;
	samplerptr->A = new int[samplerptr->chain_N];
	samplerptr->PT_alloc = true;

	samplerptr->chain_N = chain_N;//For allocation purposes, this needs to be the maximium number of chains
	allocate_sampler_mem(samplerptr);
	samplerptr->chain_N = max_chain_N_thermo_ensemble;


	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	
	int  k=0;
	//Assign initial position to start chains
	assign_initial_pos(samplerptr, initial_pos,seeding_var);	
	
	//NOTE: instead of dynamics, use variance over accept ratios over \nu steps
	//Average percent change in temperature 
	double ave_dynamics= 1.;
	//tolerance in percent change of temperature to determine equilibrium
	double tolerance= .01;
	int stability_ct = 0;
	int stability_tol = 5;
	//Frequency to check for equilibrium
	int equilibrium_check_freq=10*nu;
	
	double *old_temps = new double[samplerptr->chain_N];
	for(int i =0; i<samplerptr->chain_N; i++){
		std::cout<<samplerptr->chain_temps[i]<<std::endl;
		old_temps[i]=samplerptr->chain_temps[i];
	}

	//For each loop, we walk forward till one more swap has happened, then we update temps
	//samplerptr->N_steps = swp_freq;
	std::cout<<"Dynamical PT allocation (measured by average percent change in temperature): "<<std::endl;
	
	int t = 0;
	samplerptr->show_progress = false;
	while( t <= (N_steps-equilibrium_check_freq))
	//while(true)
	{
		if(ave_dynamics<tolerance ){
			if(stability_ct > stability_tol)
				break;
			else
				stability_ct++;
		}
		else stability_ct = 0;
		//step equilibrium_check_freq
		for(int i =0; i<equilibrium_check_freq/swp_freq; i++){
			//steps swp_freq
			//samplerptr->N_steps += swp_freq;
			PTMCMC_MH_step_incremental(samplerptr, samplerptr->swp_freq);	
			t+= samplerptr->swp_freq;
			//std::cout<<"TIME: "<<t<<std::endl;
			//Move temperatures
			update_temperatures(samplerptr, t0, nu, t);
		}
		//Calculate average percent change in temperature
		double sum = 0;
		for (int j =0; j<samplerptr->chain_N; j++){
			//std::cout<<samplerptr->chain_temps[j]<<std::endl;
			sum += std::abs((samplerptr->chain_temps[j] - old_temps[j])/old_temps[j]);
			//std::cout<<"SUM: "<<sum<<std::endl;
		}
		ave_dynamics = sum / samplerptr->chain_N;
		if(show_prog){
			printProgress(1. -  (ave_dynamics - tolerance)/(tolerance+ave_dynamics));
			std::cout<<"DYNAMICS: "<<ave_dynamics<<std::endl;
			//std::cout<<"TIME: "<<t<<std::endl;
		}
	}
	int acc, rej;
	for (int j =0; j<samplerptr->chain_N; j++){
		std::cout<<"TEMP "<<j<<": "<<samplerptr->chain_temps[j]<<std::endl;
		acc = samplerptr->swap_accept_ct[j];	
		rej = samplerptr->swap_reject_ct[j];	
		std::cout<<"Accept ratio "<<j<<": "<<(double)acc/(acc+rej)<<std::endl;
		std::cout<<"Swap attempts "<<j<<": "<<(acc+rej)<<std::endl;
		
	}

	//#################################################################
	//
	//Copy sampler to new sampler and run loop, skip checkpoint nonsense
	//
	//#################################################################
	
	//if(chain_distribution_scheme =="cold"){
	//	for(int i =max_chain_N_thermo_ensemble;i<chain_N; i++){
	//		samplerptr->chain_temps[i] = 1;
	//	}
	//}
	//write_file(chain_filename, samplerptr->output[0], N_steps,samplerptr->dimension);
	
	sampler static_sampler;
	initiate_full_sampler(&static_sampler, samplerptr, max_chain_N_thermo_ensemble, chain_N, chain_distribution_scheme);

	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_file(chain_filename, samplerptr->output[0], samplerptr->N_steps,samplerptr->dimension);
	delete [] old_temps;
	delete [] samplerptr->A;
	deallocate_sampler_mem(samplerptr);

	static_sampler.show_progress=show_prog;
	static_sampler.pool=pool;

	samplerptr = &static_sampler;	
	write_checkpoint_file(samplerptr, checkpoint_file);
	//PTMCMC_MH_loop(&static_sampler);

	end =clock();
	wend =omp_get_wtime();

	static_sampler.time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	static_sampler.time_elapsed_wall = (double)(wend-wstart);

	std::cout<<std::endl;
	
	acend =clock();
	wacend =omp_get_wtime();
	static_sampler.time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	static_sampler.time_elapsed_wall_ac = (double)(wacend - wend);

	std::cout<<std::endl;
	//double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	//double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	//std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	//accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	//rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	//std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	//double nansum=0;
	//for (int i =0; i< chain_N; i++)
	//	nansum+= samplerptr->nan_counter[i];
	//std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	

	if(checkpoint_file !=""){
		write_checkpoint_file(&static_sampler, checkpoint_file);
	}

	//free(step_accepted);
	//free(step_rejected);
	deallocate_sampler_mem(&static_sampler);
	//##################################################################
	//##################################################################
	//##################################################################
	//std::string internal_checkpoint;	
	//if(checkpoint_file !=""){
	//	internal_checkpoint = checkpoint_file;
	//	write_checkpoint_file(samplerptr, checkpoint_file);
	//}
	//else {
	//	internal_checkpoint = "temp_checkpoint_file.csv";;
	//	write_checkpoint_file(samplerptr, internal_checkpoint);
	//}

	////write_file(chain_filename, samplerptr->output[0], N_steps,samplerptr->dimension);
	//deallocate_sampler_mem(samplerptr);
	//continue_MCMC_MH_internal(internal_checkpoint, output, N_steps, swp_freq, log_prior,
	//	log_likelihood, fisher, numThreads, pool, show_prog, statistics_filename,
	//	 chain_filename,auto_corr_filename, checkpoint_file);
	//if(checkpoint_file =="")
	//	remove("temp_checkpoint_file.csv");
	//delete [] samplerptr->chain_temps;
}

/*!\brief Generic sampler, where the likelihood, prior are parameters supplied by the user
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
 *
 * 3 modes to use - 
 *
 * single threaded (numThreads = 1) runs single threaded
 *
 * multi-threaded ``deterministic'' (numThreads>1 ; pool = false) progresses each chain in parallel for swp_freq steps, then waits for all threads to complete before swapping temperatures in sequenctial order (j, j+1) then (j+1, j+2) etc (sequenctially)
 *
 * multi-threaded ``stochastic'' (numThreads>2 ; pool = true) progresses each chain in parallel by queueing each temperature and evaluating them in the order they were submitted. Once finished, the threads are queued to swap, where they swapped in the order they are submitted. This means the chains are swapped randomly, and the chains do NOT finish at the same time. The sampler runs until the the 0th chain reaches the step number
 *
 * Note on limits: In the prior function, if a set of parameters should be disallowed, return -std::numeric_limits<double>::infinity()  -- (this is in the <limits> file in std)
 *
 * Format for the auto_corr file (compatable with csv, dat, txt extensions): each row is a dimension of the cold chain, with the first row being the lengths used for the auto-corr calculation:
 *
 * lengths: length1 , length2 ...
 *
 * dim1: length1 , length2 ...
 *
 * .
 *
 * .
 *
 * .
 *
 *
 * Format for the chain file (compatable with csv, dat, txt extensions): each row is a step, each column a dimension:
 *
 * Step1: dim1 , dim2 , ...
 *
 * Step2: dim1 , dim2 , ...
 *
 * .
 *
 * .
 *
 * .
 *
 * Statistics_filename : should be txt extension
 *
 * checkpoint_file : This file saves the final position of all the chains, as well as other metadata, and can be loaded by the function <FUNCTION> to continue the chain from the point it left off. Not meant to be read by humans, the data order is custom to this software library. An empty string ("") means no checkpoint will be saved. For developers, the contents are:
 *
 * dimension, # of chains
 *
 * temps of chains
 *
 * Stepping widths of all chains
 *
 * Final position of all chains
 */
void PTMCMC_MH_internal(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	std::function<double(double*,int,int)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int,int)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int,double**,int)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler sampler;
	samplerptr = &sampler;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = swp_freq;
	samplerptr->chain_temps = chain_temps;
	samplerptr->chain_N = chain_N;
	samplerptr->N_steps = N_steps;
	samplerptr->dimension = dimension;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;

	//NOTE: currently, update the history every step until length is 
	//reached, then the history is updated every 20th step, always only
	//keeping history of length history_length (overwrites the list as 
	//it walks forward when it reaches the end)
	samplerptr->history_length = 500;
	samplerptr->history_update = 5;
	//Number of steps to take with the fisher before updating the fisher 
	//to a new value 
	//NOTE: if this is too low, detailed balance isn't maintained without 
	//accounting for the changing fisher (doesn't cancel in MH ratio)
	//but if the number is high enough, detailed balance is approximately 
	//kept without calculating second fisher
	samplerptr->fisher_update_number = 200;

	samplerptr->output = output;
	samplerptr->pool = pool;
	samplerptr->numThreads = numThreads;

	allocate_sampler_mem(samplerptr);

	//########################################################
	//########################################################
	//Set chains with temp 1 to highest priority
	for(int i =0 ;i<samplerptr->chain_N; i++){
		if(samplerptr->chain_temps[i]==1)
			samplerptr->priority[i] = 0;
	}
	//########################################################
	//########################################################

	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	

	int  k=0;
	//int swp_accepted=0, swp_rejected=0;
	int *step_accepted = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	int *step_rejected = (int *)malloc(sizeof(int)*samplerptr->chain_N);
	for(int j=0; j<samplerptr->chain_N; j++){
		step_accepted[j]=0;
		step_rejected[j]=0;
	}
	
	//samplerptr = &sampler;
	//Assign initial position to start chains
	assign_initial_pos(samplerptr, initial_pos, seeding_var);	
	//if(!seeding_var){ 
	//	for (int j=0;j<samplerptr->chain_N;j++){
	//		samplerptr->de_primed[j]=false;
	//		for (int i = 0; i<samplerptr->dimension; i++)
	//		{
	//			//Only doing this last loop because there is sometimes ~5 elements 
	//			//not initialized on the end of the output, which screw up plotting
	//			for(int l =0; l<samplerptr->N_steps; l++)
	//				samplerptr->output[j][l][i] = initial_pos[i];
	//			
	//		}
	//		samplerptr->current_likelihoods[j] =
	//			 samplerptr->ll(samplerptr->output[j][0],samplerptr->dimension, j)/samplerptr->chain_temps[j];
	//		step_accepted[j]=0;
	//		step_rejected[j]=0;
	//	}
	//}
	////Seed non-zero chains normally with variance as specified
	//else{
	//	int attempts = 0;
	//	int max_attempts = 10;
	//	double temp_pos[samplerptr->dimension];
	//	for (int j=0;j<samplerptr->chain_N;j++){
	//		samplerptr->de_primed[j]=false;
	//		if(j == 0){
	//			for (int i = 0; i<samplerptr->dimension; i++)
	//			{
	//			//Only doing this last loop because there is sometimes ~5 elements 
	//			//not initialized on the end of the output, which screw up plotting
	//				for(int l =0; l<samplerptr->N_steps; l++)
	//					samplerptr->output[j][l][i] = initial_pos[i];
	//			}
	//		}
	//		else{
	//			do{
	//				for(int i =0; i<samplerptr->dimension; i++){
	//					temp_pos[i] = gsl_ran_gaussian(samplerptr->rvec[j],seeding_var[i]) + initial_pos[i];
	//				}
	//				attempts+=1;
	//			}while(samplerptr->lp(temp_pos, samplerptr->dimension,j) == limit_inf && attempts<max_attempts);
	//			attempts =0;
	//			if(samplerptr->lp(temp_pos, samplerptr->dimension,j) != limit_inf ){
	//				for(int i =0; i<samplerptr->dimension;i++){
	//					for(int l =0; l<samplerptr->N_steps; l++)
	//						samplerptr->output[j][l][i] = temp_pos[i];
	//				}
	//			}
	//			else{
	//				for(int i =0; i<samplerptr->dimension;i++){
	//					for(int l =0; l<samplerptr->N_steps; l++)
	//						samplerptr->output[j][l][i] = initial_pos[i];
	//				}
	//		
	//			}
	//		}
	//		
	//		samplerptr->current_likelihoods[j] =
	//			 samplerptr->ll(samplerptr->output[j][0],samplerptr->dimension, j)/samplerptr->chain_temps[j];
	//		step_accepted[j]=0;
	//		step_rejected[j]=0;
	//	}
	//}

		
	PTMCMC_MH_loop(samplerptr);	
	
	//############################################################
	//Write ll lp to file
	//write_file("testing/data/mcmc_ll_lp.csv",samplerptr->ll_lp_output[0],samplerptr->N_steps,2);
	//write_file("testing/data/mcmc_ll_lp_hot.csv",samplerptr->ll_lp_output[samplerptr->chain_N-1],samplerptr->N_steps,2);
	//############################################################
	//##############################################################
	int swp_accepted=0, swp_rejected=0;
	for (int i =0;i<samplerptr->chain_N; i++)
	{
		swp_accepted+=samplerptr->swap_accept_ct[i];
		swp_rejected+=samplerptr->swap_reject_ct[i];
		step_accepted[i]+=samplerptr->step_accept_ct[i];
		step_rejected[i]+=samplerptr->step_reject_ct[i];
	}
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	std::cout<<std::endl;
	
	//###########################################################
	//Auto-correlation
	if(auto_corr_filename != ""){
		std::cout<<"Calculating Autocorrelation: "<<std::endl;
		int segments = 50;
		double target_corr = .01;
		write_auto_corr_file_from_data(auto_corr_filename, samplerptr->output[0],samplerptr->N_steps,samplerptr->dimension,segments, target_corr, samplerptr->num_threads);
	}
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	std::cout<<std::endl;
	double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	double nansum=0;
	for (int i =0; i< chain_N; i++)
		nansum+= samplerptr->nan_counter[i];
	std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_file(chain_filename, samplerptr->output[0], samplerptr->N_steps,samplerptr->dimension);

	if(checkpoint_file !=""){
		write_checkpoint_file(samplerptr, checkpoint_file);
	}

	free(step_accepted);
	free(step_rejected);
	deallocate_sampler_mem(samplerptr);
}

/*! \brief Routine to take a checkpoint file and begin a new chain at said checkpoint
 *
 * See MCMC_MH_internal for more details of parameters (pretty much all the same)
 */
void continue_PTMCMC_MH_internal(std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	std::function<double(double*,int,int)> log_prior,/**< std::function for the log_prior function -- takes double *position, int dimension, int chain_id*/
	std::function<double(double*,int,int)> log_likelihood,/**< std::function for the log_likelihood function -- takes double *position, int dimension, int chain_id*/
	std::function<void(double*,int,double**,int)>fisher,/**< std::function for the fisher function -- takes double *position, int dimension, double **output_fisher, int chain_id*/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();

	sampler sampler;
	samplerptr = &sampler;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		samplerptr->fisher_exist = false;
	}
	else 
		samplerptr->fisher_exist = true;
	
	
	//Construct sampler structure
	samplerptr->lp = log_prior;
	samplerptr->ll = log_likelihood;
	samplerptr->fish = fisher;
	samplerptr->swp_freq = swp_freq;
	samplerptr->N_steps = N_steps;
	samplerptr->show_progress = show_prog;
	samplerptr->num_threads = numThreads;

	//NOTE: currently, update the history every step until length is 
	//reached, then the history is updated every 20th step, always only
	//keeping history of length history_length (overwrites the list as 
	//it walks forward when it reaches the end)
	samplerptr->history_length = 500;
	samplerptr->history_update = 5;
	//Number of steps to take with the fisher before updating the fisher 
	//to a new value 
	//NOTE: if this is too low, detailed balance isn't maintained without 
	//accounting for the changing fisher (doesn't cancel in MH ratio)
	//but if the number is high enough, detailed balance is approximately 
	//kept without calculating second fisher
	samplerptr->fisher_update_number = 200;

	samplerptr->output = output;
	samplerptr->pool = pool;

	//Unpack checkpoint file -- allocates memory internally -- separate call unneccessary
	load_checkpoint_file(start_checkpoint_file, samplerptr);


	//allocate other parameters
	for (int chain_index=0; chain_index<samplerptr->chain_N; chain_index++)
		assign_probabilities(samplerptr, chain_index);
	for (int j=0;j<samplerptr->chain_N;j++){
		samplerptr->current_likelihoods[j] =
			 samplerptr->ll(samplerptr->output[j][0],samplerptr->dimension, j)/samplerptr->chain_temps[j];
		//std::cout<<samplerptr->current_likelihoods[j]<<std::endl;
		//step_accepted[j]=0;
		//step_rejected[j]=0;
	}
	
	//########################################################
	//########################################################
	//Set chains with temp 1 to highest priority
	for(int i =0 ;i<samplerptr->chain_N; i++){
		if(samplerptr->chain_temps[i]==1)
			samplerptr->priority[i] = 0;
	}
	//########################################################
	//########################################################
	
	//########################################################
	PTMCMC_MH_loop(samplerptr);	
	//##############################################################
	
	//int swp_accepted=0, swp_rejected=0;
	//for (int i =0;i<samplerptr->chain_N; i++)
	//{
	//	swp_accepted+=samplerptr->swap_accept_ct[i];
	//	swp_rejected+=samplerptr->swap_reject_ct[i];
	//	step_accepted[i]+=samplerptr->step_accept_ct[i];
	//	step_rejected[i]+=samplerptr->step_reject_ct[i];
	//}
	end =clock();
	wend =omp_get_wtime();

	samplerptr->time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall = (double)(wend-wstart);

	std::cout<<std::endl;
	
	//###########################################################
	//Auto-correlation
	if(auto_corr_filename != ""){
		std::cout<<"Calculating Autocorrelation: "<<std::endl;
		int segments = 50;
		double target_corr = .01;
		write_auto_corr_file_from_data(auto_corr_filename, samplerptr->output[0],samplerptr->N_steps,samplerptr->dimension,segments, target_corr, samplerptr->num_threads);
	}
	//###########################################################
	acend =clock();
	wacend =omp_get_wtime();
	samplerptr->time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	samplerptr->time_elapsed_wall_ac = (double)(wacend - wend);

	//std::cout<<std::endl;
	//double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	//double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	//std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	//accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	//rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	//std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	//std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	//double nansum=0;
	//for (int i =0; i< chain_N; i++)
	//	nansum+= samplerptr->nan_counter[i];
	//std::cout<<"NANS in Fisher Calculations (all chains): "<<nansum<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(samplerptr, statistics_filename);
	
	if(chain_filename != "")
		write_file(chain_filename, samplerptr->output[0], samplerptr->N_steps,samplerptr->dimension);

	if(end_checkpoint_file !=""){
		write_checkpoint_file(samplerptr, end_checkpoint_file);
	}

	//free(step_accepted);
	//free(step_rejected);
	//temps usually allocated by user, but for continued chains, this is done internally
	free(samplerptr->chain_temps);
	deallocate_sampler_mem(samplerptr);
}
					
/*!\brief Internal function that runs the actual loop for the sampler -- increment version
 *
 * The regular loop function runs for the entire range, this increment version will only step ``increment'' steps -- asynchronous: steps are measured by the 0th chain
 *
 */
void PTMCMC_MH_step_incremental(sampler *sampler, int increment)
{
	//Make sure we're not going out of memory
	if (sampler->progress + increment > sampler->N_steps)
		increment = sampler->N_steps - sampler->progress;
	//Sampler Loop - ``Deterministic'' swapping between chains
	if (!sampler->pool)
	{
		int k =0;
		int cutoff ;
		int step_log;
		omp_set_num_threads(sampler->num_threads);
		#pragma omp parallel 
		{
		while (k<(increment-1) ){
			#pragma omp for
			for (int j=0; j<sampler->chain_N; j++)
			{
				if( sampler->N_steps-sampler->chain_pos[j] <= sampler->swp_freq) 
					cutoff = sampler->N_steps-sampler->chain_pos[j]-1;	
				else cutoff = sampler->swp_freq;	
				if(j==0)
					step_log = cutoff;
				for (int i = 0 ; i< cutoff;i++)
				{
					int success;
					if(!sampler->RJMCMC){
						success = mcmc_step(sampler, sampler->output[j][sampler->chain_pos[j]], sampler->output[j][sampler->chain_pos[j]+1],sampler->param_status[j][0],sampler->param_status[j][0],j);	
					}
					else{
						success = mcmc_step(sampler, sampler->output[j][sampler->chain_pos[j]], sampler->output[j][sampler->chain_pos[j]+1],sampler->param_status[j][sampler->chain_pos[j]],sampler->param_status[j][sampler->chain_pos[j]+1],j);	
					}
					//success = mcmc_step(sampler, sampler->output[j][sampler->chain_pos[j]], sampler->output[j][sampler->chain_pos[j]+1],j);	
					sampler->chain_pos[j]+=1;
					//if(success==1){step_accepted[j]+=1;}
					//else{step_rejected[j]+=1;}
					if(success==1){sampler->step_accept_ct[j]+=1;}
					else{sampler->step_reject_ct[j]+=1;}
					if(!sampler->de_primed[j])
						update_history(sampler,sampler->output[j][sampler->chain_pos[j]], j);
					else if(sampler->chain_pos[j]%sampler->history_update==0)
						update_history(sampler,sampler->output[j][sampler->chain_pos[j]], j);
					//Log LogLikelihood and LogPrior	
					if(sampler->log_ll){
						samplerptr->ll_lp_output[j][sampler->chain_pos[j]][0]= 
							samplerptr->current_likelihoods[j];
					}
					if(sampler->log_lp){
						samplerptr->ll_lp_output[j][sampler->chain_pos[j]][1]= 
							samplerptr->lp(
							samplerptr->output[j][sampler->chain_pos[j]],
							samplerptr->dimension, j);
					}
					
				}
				if(!sampler->de_primed[j]) 
				{
					if ((sampler->chain_pos[j])>sampler->history_length)
					{
						sampler->de_primed[j]=true;
						assign_probabilities(sampler,j);	
					}
				}
			}
			#pragma omp single
			{
				k+= step_log;
				sampler->progress+=step_log;
				int swp_accepted=0, swp_rejected=0;
				///NEEDS TO CHANGE
				//chain_swap(sampler, sampler->output, sampler->chain_pos[0], &swp_accepted, &swp_rejected);
				for(int i =0 ; i<sampler->chain_N-1; i++){
					int k = sampler->chain_pos[i];
					int l = sampler->chain_pos[i+1];
					int success;
					success = single_chain_swap(sampler, sampler->output[i][k], sampler->output[i+1][l],i,i+1);
					//success = -1;
					if(success ==1){
						sampler->swap_accept_ct[i]+=1;	
						sampler->swap_accept_ct[i+1]+=1;	
						if(sampler->PT_alloc)
							sampler->A[i+1] = 1;
					}
					else{
						sampler->swap_reject_ct[i]+=1;	
						sampler->swap_reject_ct[+1]+=1;	
						if(sampler->PT_alloc)
							sampler->A[i+1] = 0;
					}
				}
				//sampler->swap_accept_ct+=swp_accepted;
				//sampler->swap_reject_ct+=swp_rejected;
				if(sampler->show_progress)
					printProgress((double)sampler->progress/sampler->N_steps);	
			}
		}
		}
	}

	//POOLING  -- ``Stochastic'' swapping between chains
	else
	{
		ThreadPool pool(sampler->num_threads);
		poolptr = &pool;
		while(sampler->progress<increment-1)
		{
			for(int i =0; i<sampler->chain_N; i++)
			{
				if(sampler->waiting[i]){
					if(sampler->chain_pos[i]<(sampler->N_steps-1))
					{
						sampler->waiting[i]=false;
						//if(i==0) samplerptr->progress+=samplerptr->swp_freq;
						poolptr->enqueue(i);
					}
					//If a chain finishes before chain 0, it's wrapped around 
					//and allowed to keep stepping at low priority-- 
					//not sure if this is the best
					//method for keeping the 0th chain from finishing last or not
					else if(i !=0){

						sampler->waiting[i]=false;
						std::cout<<"Chain "<<i<<" finished-- being reset"<<std::endl;
						sampler->priority[i] = 2;
						int pos = sampler->chain_pos[i];
						for (int k =0; k<sampler->dimension; k++){
							sampler->output[i][0][k] = 
								sampler->output[i][pos][k] ;
						}
						sampler->chain_pos[i] = 0;

						poolptr->enqueue(i);
					}
				}
				
			}
			if(sampler->show_progress)
				printProgress((double)sampler->progress/sampler->N_steps);
			//usleep(300);
		}
	}
}

/*!\brief Internal function that runs the actual loop for the sampler
 *
 */
void PTMCMC_MH_loop(sampler *sampler)
{
	int k =0;
	int cutoff ;
	//Sampler Loop - ``Deterministic'' swapping between chains
	if (!sampler->pool)
	{
		omp_set_num_threads(sampler->num_threads);
		#pragma omp parallel 
		{
		while (k<sampler->N_steps-1){
			#pragma omp single
			{
				if( sampler->N_steps-k <= sampler->swp_freq) 
					cutoff = sampler->N_steps-k-1;	
				else cutoff = sampler->swp_freq;	
			}
			#pragma omp for
			for (int j=0; j<sampler->chain_N; j++)
			{
				for (int i = 0 ; i< cutoff;i++)
				{
					int success;
					if(!sampler->RJMCMC){
						success = mcmc_step(sampler, sampler->output[j][k+i], sampler->output[j][k+i+1],sampler->param_status[j][0],sampler->param_status[j][0],j);	
					}
					else{
						success = mcmc_step(sampler, sampler->output[j][k+i], sampler->output[j][k+i+1],sampler->param_status[j][k+i],sampler->param_status[j][k+i+1],j);	
					}
					sampler->chain_pos[j]+=1;
					//if(success==1){step_accepted[j]+=1;}
					//else{step_rejected[j]+=1;}
					if(success==1){sampler->step_accept_ct[j]+=1;}
					else{sampler->step_reject_ct[j]+=1;}
					if(!sampler->de_primed[j])
						update_history(sampler,sampler->output[j][k+i+1], j);
					else if(sampler->chain_pos[j]%sampler->history_update==0)
						update_history(sampler,sampler->output[j][k+i+1], j);
					//Log LogLikelihood and LogPrior	
					if(sampler->log_ll){
						samplerptr->ll_lp_output[j][k+i+1][0]= 
							samplerptr->current_likelihoods[j];
					}
					if(sampler->log_lp){
						samplerptr->ll_lp_output[j][k+i+1][1]= 
							samplerptr->lp(
							samplerptr->output[j][k+i+1],
							samplerptr->dimension, j);
					}
					
				}
				if(!sampler->de_primed[j]) 
				{
					if ((k+cutoff)>sampler->history_length)
					{
						sampler->de_primed[j]=true;
						assign_probabilities(sampler,j);	
					}
				}
			}
			#pragma omp single
			{
				k+= cutoff;
				int swp_accepted=0, swp_rejected=0;
				chain_swap(sampler, sampler->output, k, &swp_accepted, &swp_rejected);
				//sampler->swap_accept_ct+=swp_accepted;
				//sampler->swap_reject_ct+=swp_rejected;
				if(sampler->show_progress)
					printProgress((double)k/sampler->N_steps);	
			}
		}
		}
	}

	//POOLING  -- ``Stochastic'' swapping between chains
	else
	{
		ThreadPool pool(sampler->num_threads);
		poolptr = &pool;
		while(sampler->progress<sampler->N_steps-1)
		{
			for(int i =0; i<sampler->chain_N; i++)
			{
				if(sampler->waiting[i]){
					if(sampler->chain_pos[i]<(sampler->N_steps-1))
					{
						sampler->waiting[i]=false;
						//if(i==0) samplerptr->progress+=samplerptr->swp_freq;
						poolptr->enqueue(i);
					}
					//If a chain finishes before chain 0, it's wrapped around 
					//and allowed to keep stepping at low priority-- 
					//not sure if this is the best
					//method for keeping the 0th chain from finishing last or not
					else if(sampler->chain_temps[i] !=1){

						sampler->waiting[i]=false;
						std::cout<<"Chain "<<i<<" finished-- being reset"<<std::endl;
						sampler->priority[i] = 2;
						int pos = sampler->chain_pos[i];
						for (int k =0; k<sampler->dimension; k++){
							sampler->output[i][0][k] = 
								sampler->output[i][pos][k] ;
						}
						sampler->chain_pos[i] = 0;

						poolptr->enqueue(i);
					}
					else{
						//If 0 T chain, just wait till everything else is done
						sampler->waiting[i] = false;	
					}
				}
				
			}
			if(sampler->show_progress)
				printProgress((double)sampler->progress/sampler->N_steps);
			//usleep(300);
		}
	}
}


void mcmc_step_threaded(int j)
{
	int k = samplerptr->chain_pos[j];
	int cutoff;
	if( samplerptr->N_steps-k <= samplerptr->swp_freq) cutoff = samplerptr->N_steps-k-1;	
	else cutoff = samplerptr->swp_freq;	
	for (int i = 0 ; i< cutoff;i++)
	{
		int success;
		//success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],j);	
		if(!samplerptr->RJMCMC){
			success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],samplerptr->param_status[j][0],samplerptr->param_status[j][0],j);	
		}
		else{
			success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],samplerptr->param_status[j][k+i],samplerptr->param_status[j][k+i+1],j);	
		}
	
		if(success==1){samplerptr->step_accept_ct[j]+=1;}
		else{samplerptr->step_reject_ct[j]+=1;}
		if(!samplerptr->de_primed[j])
			update_history(samplerptr,samplerptr->output[j][k+i+1], j);
		else if(samplerptr->chain_pos[j]%samplerptr->history_update==0)
			update_history(samplerptr,samplerptr->output[j][k+i+1], j);
		//##############################################################
		//Trouble shooting -- save lp and ll
		if(samplerptr->log_ll){
			samplerptr->ll_lp_output[j][k+i+1][0]= 
				samplerptr->current_likelihoods[j];
		}
		if(samplerptr->log_lp){
			samplerptr->ll_lp_output[j][k+i+1][1]= 
				samplerptr->lp(samplerptr->output[j][k+i+1],
				samplerptr->dimension, j);
		}
		//##############################################################
	}
	if(!samplerptr->de_primed[j]) 
	{
		if ((k+cutoff)>samplerptr->history_length)
		{
			samplerptr->de_primed[j]=true;
			assign_probabilities(samplerptr,j);	
		}
	}
	samplerptr->chain_pos[j]+=cutoff;
	if(j==0) samplerptr->progress+=cutoff;

	//update stepsize to maximize step efficiency
	//increases in stepsizes of 10%
	double frac, acc, rej;
	if(samplerptr->chain_pos[j]%samplerptr->check_stepsize_freq[j] == 0){
		//Gaussian
		if(samplerptr->step_prob[j][0]!= 0){
			acc = samplerptr->gauss_accept_ct[j] - samplerptr->gauss_last_accept_ct[j];	
			rej = samplerptr->gauss_reject_ct[j] - samplerptr->gauss_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][0] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][0] *=1.1;	
			}
			samplerptr->gauss_last_accept_ct[j]=samplerptr->gauss_accept_ct[j];
			samplerptr->gauss_last_reject_ct[j]=samplerptr->gauss_reject_ct[j];
		}	
		//de
		if(samplerptr->step_prob[j][1]!= 0){
			acc = samplerptr->de_accept_ct[j] - samplerptr->de_last_accept_ct[j];	
			rej = samplerptr->de_reject_ct[j] - samplerptr->de_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][1] *=1.1;	
			}
			samplerptr->de_last_accept_ct[j]=samplerptr->de_accept_ct[j];
			samplerptr->de_last_reject_ct[j]=samplerptr->de_reject_ct[j];
		}	
		//fisher
		if(samplerptr->step_prob[j][3]!= 0){
			acc = samplerptr->fish_accept_ct[j] - samplerptr->fish_last_accept_ct[j];	
			rej = samplerptr->fish_reject_ct[j] - samplerptr->fish_last_reject_ct[j];	
			frac = acc / (acc + rej);
			if(frac<samplerptr->min_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3] *=.9;	
			}
			else if(frac>samplerptr->max_target_accept_ratio[j]){
				samplerptr->randgauss_width[j][3] *=1.1;	
			}
			samplerptr->fish_last_accept_ct[j]=samplerptr->fish_accept_ct[j];
			samplerptr->fish_last_reject_ct[j]=samplerptr->fish_reject_ct[j];
		}	
		
	}
	poolptr->enqueue_swap(j);
}
void mcmc_swap_threaded(int i, int j)
{
	int k = samplerptr->chain_pos[i];
	int l = samplerptr->chain_pos[j];
	int success;
	success = single_chain_swap(samplerptr, samplerptr->output[i][k], samplerptr->output[j][l],i,j);
	//success = -1;
	if(success ==1){
		samplerptr->swap_accept_ct[i]+=1;	
		samplerptr->swap_accept_ct[j]+=1;	
	}
	else{
		samplerptr->swap_reject_ct[i]+=1;	
		samplerptr->swap_reject_ct[j]+=1;	
	}
	samplerptr->waiting[i]=true;
	samplerptr->waiting[j]=true;
}
void PTMCMC_MH(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	double (*log_prior)(double *param, int dimension),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	auto ll = [&log_likelihood](double *param, int dim, int chain_id){
		return log_likelihood(param, dim);};

	auto lp = [&log_prior](double *param, int dim, int chain_id){
		return log_prior(param, dim);};
	std::function<void(double*,int,double**,int)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int dim, double **fisherm, int chain_id){
			fisher(param, dim, fisherm);};
	}
	
	PTMCMC_MH_internal(output, dimension, N_steps, chain_N, initial_pos, seeding_var,chain_temps, swp_freq, 
			lp, ll, f, numThreads, pool, show_prog, 
			statistics_filename, chain_filename, auto_corr_filename, checkpoint_file);
}
void PTMCMC_MH(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain*/
	int chain_N,	/**< Number of chains*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps,	/**<Double array of temperatures for the chains*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	double (*log_prior)(double *param, int dimension, int chain_id),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension, int chain_id),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	std::function<double(double*,int,int)> lp = log_prior;
	std::function<double(double*,int,int)> ll = log_likelihood;
	std::function<void(double*,int,double**,int)>f = fisher;
	PTMCMC_MH_internal(output, dimension, N_steps, chain_N, initial_pos, seeding_var,chain_temps, swp_freq, 
			lp, ll, f, numThreads, pool, show_prog, 
			statistics_filename, chain_filename, auto_corr_filename, checkpoint_file);
}
void continue_PTMCMC_MH(std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	double (*log_prior)(double *param, int dimension, int chain_id),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension, int chain_id),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{

	std::function<double(double*,int,int)> lp = log_prior;
	std::function<double(double*,int,int)> ll = log_likelihood;
	std::function<void(double*,int,double**,int)>f = fisher;
	continue_PTMCMC_MH_internal(start_checkpoint_file,
			output,
			N_steps,
			swp_freq,
			lp,
			ll,
			f,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			auto_corr_filename,
			end_checkpoint_file);
}
void continue_PTMCMC_MH(std::string start_checkpoint_file,/**< File for starting checkpoint*/
	double ***output,/**< [out] output array, dimensions: output[chain_N][N_steps][dimension]*/
	int N_steps,/**< Number of new steps to take*/
	int swp_freq,/**< frequency of swap attempts between temperatures*/
	double (*log_prior)(double *param, int dimension),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads,/**<Number of threads to use*/
	bool pool,/**<Boolean for whether to use ``deterministic'' vs ``stochastic'' sampling*/
	bool show_prog,/**< Boolean for whether to show progress or not (turn off for cluster runs*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string auto_corr_filename,/**< Filename to output auto correlation in some interval, if empty string, not output*/
	std::string end_checkpoint_file/**< Filename to output data for checkpoint at the end of the continued run, if empty string, not saved*/
	)
{

	auto ll = [&log_likelihood](double *param, int dim, int chain_id){
		return log_likelihood(param, dim);};

	auto lp = [&log_prior](double *param, int dim, int chain_id){
		return log_prior(param, dim);};
	std::function<void(double*,int,double**,int)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int dim, double **fisherm, int chain_id){
			fisher(param, dim, fisherm);};
	}
	continue_PTMCMC_MH_internal(start_checkpoint_file,
			output,
			N_steps,
			swp_freq,
			lp,
			ll,
			f,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			auto_corr_filename,
			end_checkpoint_file);
}
void PTMCMC_MH_dynamic_PT_alloc(double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, int dimension),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	auto ll = [&log_likelihood](double *param, int dim, int chain_id){
		return log_likelihood(param, dim);};

	auto lp = [&log_prior](double *param, int dim, int chain_id){
		return log_prior(param, dim);};
	std::function<void(double*,int,double**,int)> f =NULL;
	if(fisher){
		f = [&fisher](double *param, int dim, double **fisherm, int chain_id){
			fisher(param, dim, fisherm);};
	}
	PTMCMC_MH_dynamic_PT_alloc_internal(output,
			dimension,
			N_steps,
			chain_N,
			max_chain_N_thermo_ensemble,
			initial_pos,
			seeding_var,
			chain_temps,
			swp_freq,
			t0,
			nu,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			checkpoint_file);

}
void PTMCMC_MH_dynamic_PT_alloc(double ***output, /**< [out] Output chains, shape is double[max_chain_N, N_steps,dimension]*/
	int dimension, 	/**< dimension of the parameter space being explored*/
	int N_steps,	/**< Number of total steps to be taken, per chain AFTER chain allocation*/
	int chain_N,/**< Maximum number of chains to use */
	int max_chain_N_thermo_ensemble,/**< Maximum number of chains to use in the thermodynamic ensemble (may use less)*/
	double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
	double *seeding_var, 	/**<Variance of the normal distribution used to seed each chain higher than 0 - shape double[dimension]*/
	double *chain_temps, /**< Final chain temperatures used -- should be shape double[chain_N]*/
	int swp_freq,	/**< the frequency with which chains are swapped*/
	int t0,/**< Time constant of the decay of the chain dynamics  (~1000)*/
	int nu,/**< Initial amplitude of the dynamics (~100)*/
	std::string chain_distribution_scheme, /*How to allocate the remaining chains once equilibrium is reached*/
	double (*log_prior)(double *param, int dimension, int chain_id),	/**<Funcion pointer for the log_prior*/
	double (*log_likelihood)(double *param, int dimension, int chain_id),	/**<Function pointer for the log_likelihood*/
	void (*fisher)(double *param, int dimension, double **fisher, int chain_id),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
	int numThreads, /**< Number of threads to use (=1 is single threaded)*/
	bool pool, /**< boolean to use stochastic chain swapping (MUST have >2 threads)*/
	bool show_prog, /**< boolean whether to print out progress (for example, should be set to ``false'' if submitting to a cluster)*/
	std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
	std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
	std::string checkpoint_file/**< Filename to output data for checkpoint, if empty string, not saved*/
	)
{
	std::function<double(double*,int,int)> lp = log_prior;
	std::function<double(double*,int,int)> ll = log_likelihood;
	std::function<void(double*,int,double**,int)>f = fisher;
	PTMCMC_MH_dynamic_PT_alloc_internal(output,
			dimension,
			N_steps,
			chain_N,
			max_chain_N_thermo_ensemble,
			initial_pos,
			seeding_var,
			chain_temps,
			swp_freq,
			t0,
			nu,
			chain_distribution_scheme,
			lp,
			ll,
			f,
			numThreads,
			pool,
			show_prog,
			statistics_filename,
			chain_filename,
			checkpoint_file);

}
