#include "mcmc_sampler.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "util.h"
#include "mcmc_sampler_internals.h"
#include <omp.h>
#include <time.h>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#include <queue>

#ifndef _OPENMP
#define omp ignore
#endif

#include <unistd.h>

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


	//void enqueue(Task task)
	//{
	//	{
	//		std::unique_lock<std::mutex> lock{mEventMutex};
	//		mTasks.emplace(std::move(task));
	//	}
	//	mEventVar.notify_one();
	//}
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
	
	int numSwpThreads= 2;
		
	std::condition_variable mEventVarSWP;

	std::mutex mEventMutexSWP;

	//std::queue<Task> mTasks;
	std::queue<int> mTasks;
	std::queue<int> mSwaps;

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
						j = std::move(mTasks.front());
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
						
						//if (mStopping && mSwaps.empty())
						if (mStopping && mSwaps.size()<2)
							break;	
						j = std::move(mSwaps.front());
						mSwaps.pop();
						k = std::move(mSwaps.front());
						mSwaps.pop();
						//std::cout<<mTasks.empty();
					}
					mcmc_swap_threaded(j,k);
					
				}
			});
		}
	}
	void stop() noexcept
	{
		std::cout<<"STOPPING"<<std::endl;
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

/*!\brief Generic sampler, where the likelihood, prior are parameters supplied by the user
 *
 * Base of the sampler, generic, with user supplied quantities for most of the samplers
 * properties
 * 	
 * Uses the Metropolis-Hastings method, with the option for Fisher/MALA steps if the Fisher
 * routine is supplied.
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
 */
void MCMC_MH(	double ***output, /**< [out] Output chains, shape is double[chain_N, N_steps,dimension]*/
		int dimension, 	/**< dimension of the parameter space being explored*/
		int N_steps,	/**< Number of total steps to be taken, per chain*/
		int chain_N,	/**< Number of chains*/
		double *initial_pos, 	/**<Initial position in parameter space - shape double[dimension]*/
		double *chain_temps,	/**<Double array of temperatures for the chains*/
		int swp_freq,	/**< the frequency with which chains are swapped*/
		double (*log_prior)(double *param, int dimension),	/**<Funcion pointer for the log_prior*/
		double (*log_likelihood)(double *param, int dimension),	/**<Function pointer for the log_likelihood*/
		void (*fisher)(double *param, int dimension, double **fisher),	/**<Function pointer for the fisher - if NULL, fisher steps are not used*/
		std::string statistics_filename,/**< Filename to output sampling statistics, if empty string, not output*/
		std::string chain_filename,/**< Filename to output data (chain 0 only), if empty string, not output*/
		std::string auto_corr_filename/**< Filename to output auto correlation in some interval, if empty string, not output*/
		)
{
	clock_t start, end, acend;
	double wstart, wend, wacend;
	start = clock();
	wstart = omp_get_wtime();
	int numThread = 20;
	omp_set_num_threads(numThread);

	//random number generator initialization
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r = gsl_rng_alloc(T);

	//Array holding the probability of each 
	//type of step - Gaussian, differential evolution, MMALA, Fisher
	sampler sampler;

	//if Fisher is not provided, Fisher and MALA steps
	//aren't used
	if(fisher ==NULL){
		sampler.fisher_exist = false;
	}
	else 
		sampler.fisher_exist = true;
	
	
	//Construct sampler structure
	sampler.lp = log_prior;
	sampler.ll = log_likelihood;
	sampler.fish = fisher;
	sampler.swp_freq = swp_freq;
	sampler.chain_temps = chain_temps;
	sampler.chain_N = chain_N;
	sampler.N_steps = N_steps;
	sampler.dimension = dimension;
	sampler.r = r;
	sampler.history_length = 1000;
	sampler.fisher_update_number = 1000;
	sampler.output = output;
	//########################################################
	//########################################################
	//POOLING -- TESTING
	sampler.pool = false;
	//########################################################
	//########################################################

	allocate_sampler_mem(&sampler);
	for (int chain_index; chain_index<sampler.chain_N; chain_index++)
		assign_probabilities(&sampler, chain_index);
	

	//int i,j, k=0;
	//int j, k=0;
	int  k=0;
	int swp_accepted=0, swp_rejected=0;
	//int step_accepted=0, step_rejected=0;
	int *step_accepted = (int *)malloc(sizeof(int)*sampler.chain_N);
	int *step_rejected = (int *)malloc(sizeof(int)*sampler.chain_N);
	
	samplerptr = &sampler;
	//Assign initial position to start chains
	//Currently, just set all chains to same initial position
	for (int j=0;j<sampler.chain_N;j++){
		for (int i = 0; i<sampler.dimension; i++)
		{
			sampler.de_primed[j]=false;
			for(int l =0; l<N_steps; l++)
				output[j][l][i] = initial_pos[i];
		}
		step_accepted[j]=0;
		step_rejected[j]=0;
	}

		
	
	int cutoff ;
	//Sampler Loop
	if (!samplerptr->pool)
	{
		#pragma omp parallel //num_threads(4)
		{
		while (k<N_steps-1){
			#pragma omp single
			{
			if( N_steps-k <= sampler.swp_freq) cutoff = N_steps-k-1;	
			else cutoff = sampler.swp_freq;	
			}
			//#pragma omp parallel for reduction(+:step_accepted) reduction(+:step_rejected)
			//#pragma omp parallel for 
			#pragma omp for
			for (int j=0; j<chain_N; j++)
			{
				for (int i = 0 ; i< cutoff;i++)
				{
					int success;
					success = mcmc_step(&sampler, output[j][k+i], output[j][k+i+1],j);	
					if(success==1){step_accepted[j]+=1;}
					else{step_rejected[j]+=1;}
					update_history(&sampler,output[j][k+i+1], j);
				}
				if(!sampler.de_primed[j]) 
				{
					if ((k+cutoff)>sampler.history_length)
					{
						sampler.de_primed[j]=true;
						assign_probabilities(&sampler,j);	
					}
				}
			}
			#pragma omp single
			{
			k+= cutoff;
			chain_swap(&sampler, output, k, &swp_accepted, &swp_rejected);
			printProgress((double)k/N_steps);	
			}
		}
		}
	}


	else
	{
		ThreadPool pool(numThread);
		poolptr = &pool;
		while(samplerptr->progress<N_steps-samplerptr->swp_freq-2)
		{
			for(int i =0; i<samplerptr->chain_N; i++)
			{
				if(samplerptr->waiting[i] && samplerptr->chain_pos[i]<(N_steps-samplerptr->swp_freq -1))
				{
					samplerptr->waiting[i]=false;
					if(i==0) samplerptr->progress+=samplerptr->swp_freq;
					poolptr->enqueue(i);
				}
			}
			printProgress((double)samplerptr->progress/N_steps);
		}
		//##############################################################
		for (int i =0;i<samplerptr->chain_N; i++)
		{
			swp_accepted+=samplerptr->swap_accept_ct[i];
			swp_rejected+=samplerptr->swap_reject_ct[i];
			step_accepted[i]=samplerptr->step_accept_ct[i];
			step_rejected[i]=samplerptr->step_reject_ct[i];
		}
	}
	end =clock();
	wend =omp_get_wtime();

	//for (int j=0;j<sampler.chain_N;j++){
	//	std::cout<<std::endl;
	//	//for (int i = 0; i<sampler.dimension; i++)
	//	{
	//		std::cout<<exp(samplerptr->output[j][N_steps-100][0])/MSOL_SEC<<std::endl;
	//		std::cout<<samplerptr->output[j][N_steps-100][1]<<std::endl;
	//		std::cout<<samplerptr->output[j][N_steps-100][2]<<std::endl;
	//		std::cout<<samplerptr->output[j][N_steps-100][3]<<std::endl;
	//	}
	//}

	//###########################################################
	//Auto-correlation
	if(auto_corr_filename != ""){
		//Transpose the data for auto-correlation
		double **temp = (double **) malloc(sizeof(double*)*N_steps);
		for (int i = 0 ; i< sampler.dimension; i++){
			temp[i] = (double *)malloc(sizeof(double)*N_steps);
			for(int j =0; j< N_steps; j++){
				temp[i][j] = output[0][j][i];
			}
		}
		int segments = 50;
		double **ac = allocate_2D_array(sampler.dimension+1, segments);
		//First row is the step size for the given auto-corr length
		double stepsize = (double)N_steps/segments;
		for (int i =0 ; i<segments; i++)
			ac[0][i] = (int)(stepsize*(1.+i));
		
		#pragma omp parallel for 
		for (int i = 0 ; i< sampler.dimension; i++)
		{
			auto_corr_intervals(temp[i],N_steps, ac[i+1], segments, 0.01);
		}	

		write_file(auto_corr_filename, ac, sampler.dimension+1, segments);

		for (int i = 0 ; i< sampler.dimension; i++){
			free(temp[i]);
		}
		free(temp);
		deallocate_2D_array(ac,sampler.dimension+1, segments);
	}
	//###########################################################
	
	acend =clock();
	wacend =omp_get_wtime();
	sampler.time_elapsed_cpu = (double)(end-start)/CLOCKS_PER_SEC;
	sampler.time_elapsed_wall = (double)(wend-wstart);
	sampler.time_elapsed_cpu_ac = (double)(acend-end)/CLOCKS_PER_SEC;
	sampler.time_elapsed_wall_ac = (double)(wacend - wend);

	std::cout<<std::endl;
	double accepted_percent = (double)(swp_accepted)/(swp_accepted+swp_rejected);
	double rejected_percent = (double)(swp_rejected)/(swp_accepted+swp_rejected);
	std::cout<<"Accepted percentage of chain swaps (all chains): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of chain swaps (all chains): "<<rejected_percent<<std::endl;
	accepted_percent = (double)(step_accepted[0])/(step_accepted[0]+step_rejected[0]);
	rejected_percent = (double)(step_rejected[0])/(step_accepted[0]+step_rejected[0]);
	std::cout<<"Accepted percentage of steps (cold chain): "<<accepted_percent<<std::endl;
	std::cout<<"Rejected percentage of steps (cold chain): "<<rejected_percent<<std::endl;
	std::cout<<"NANS (all chains): "<<sampler.nan_counter<<std::endl;
	
	if(statistics_filename != "")
		write_stat_file(&sampler, statistics_filename, step_accepted, step_rejected,
				swp_accepted, swp_rejected);	
	
	if(chain_filename != "")
		write_file(chain_filename, output[0], N_steps,sampler.dimension);

	free(step_accepted);
	free(step_rejected);
	deallocate_sampler_mem(&sampler);
}


void mcmc_step_threaded(int j)
{
	int k = samplerptr->chain_pos[j];
	int cutoff = samplerptr->swp_freq;
	for (int i = 0 ; i< cutoff;i++)
	{
		int success;
		success = mcmc_step(samplerptr, samplerptr->output[j][k+i], samplerptr->output[j][k+i+1],j);	
		if(success==1){samplerptr->step_accept_ct[j]+=1;}
		else{samplerptr->step_reject_ct[j]+=1;}
		update_history(samplerptr,samplerptr->output[j][k+i+1], j);
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
