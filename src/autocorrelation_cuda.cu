#include "autocorrelation_cuda.h"
#include "autocorrelation_cuda.hu"
#include "util.h"
#include <iostream>
#include <condition_variable>
#include <thread>
#include <queue>
#include <functional>
#include <mutex>
#include <unistd.h>
#include <threadPool.h>
#include <cufft.h>


/*! \file
 */

GPUplan *plans_global;

/*! \brief Internal function to calculate the autocorrelation for a given lag
 * Customized for the thread pool architecture, with extra arguments because of the way the memory is allocated
 */
__device__ __host__
void auto_corr_internal(double *arr, /**< Input array of data*/
			int length, /**< Length of input array*/
			int lag,  /**< Lag to be used to calculate the correlation*/
			double average,  /**< Average of the array arr*/
			double *corr,  /**< [out] output correlation*/
			int start_id /**< ID of location to start calculation -- input arrary arr is assumed to be contiguous for multiple dimensions*/
			)
{
	double sum = 0;
	for(int i =0; i< (length - lag); i++){
		sum+= (arr[i+lag+start_id] - average ) * ( arr[i+start_id] - average );
	}		
	*corr = sum / (length - lag);
}

/*! \brief Internal function to launch the CUDA kernel for a range of autocorrelations
 * 
 * Correlation function used:
 *
 * rho(lag) = 1 / (length - lag) \sum (arr[i+lag]-average) ( arr[i]- average)
 *
 * target_corr = rho(rho_index)/rho(0) = rho(rho_index)/var
 */
__global__
void auto_corr_internal_kernal(double *arr, /**< Input array of data*/
				int length,  /**< Length of data array*/
				double average, /**< Average of input data*/
				int *rho_index, /**< [out] Index of the lag that results ina correlation ratio target_corr*/
				double target_corr, /**< Target correlation ratio rho(lag)/rho(0) = target_corr*/
				double var, /**< Variance rho(0)*/
				int start_id/**< Starting index to use for the data array arr*/
				)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < *rho_index){
		double corr;
		auto_corr_internal(arr, length, id, average, &corr, start_id);
		if(corr/var<target_corr) atomicMin(rho_index, id);
	}

}

/*! \brief Write data file for autocorrelation lengths of the data given a data file name, as written by the mcmc_sampler
 */
void write_file_auto_corr_from_data_file_accel(std::string acfile, /**< Filename of the autocorrelation data*/
					std::string chains_file, /**<Filename of the data file for the chains*/
					int dimension, /**< Dimension of the data*/
					int N_steps, /**< Number of steps in the chain*/
					int num_segments,  /**< Number of segments to check the autocorrelation length for each dimension*/
					double target_corr/**< Target correlation ratio to use for the correlation length calculation*/
					)
{
	double **chains = allocate_2D_array(N_steps, dimension);
	read_file(chains_file, chains, N_steps, dimension);
	write_file_auto_corr_from_data_accel(acfile, chains, dimension, 
			N_steps, num_segments, target_corr);	
	deallocate_2D_array(chains,N_steps, dimension);
}

/*! \brief Write data file given output chains, as formatted by the mcmc_sampler
 */
void write_file_auto_corr_from_data_accel(std::string acfile, /**< Output autocorrelation filename */
					double **chains, /**< Chain data from MCMC_sampler*/
					int dimension, /**< Dimension of the data*/
					int N_steps, /**< Number of steps in the chain*/
					int num_segments,  /**< Number of segments to check the autocorrelation length for each dimension*/
					double target_corr/**< Target correlation ratio to use for the correlation length calculation*/
					)
{
	double **autocorr = allocate_2D_array(dimension, num_segments);
	double **autocorrout = allocate_2D_array(dimension+1, num_segments);
	auto_corr_from_data_accel(chains, dimension, N_steps, num_segments, target_corr, autocorr);
	int seg_step = N_steps/ num_segments;
	for(int i =0 ; i < num_segments; i ++)
	{
		autocorrout[0][i] = (i+1) * seg_step;
		for(int j =0; j<dimension; j++){
			autocorrout[j+1][i] = autocorr[j][i];
		}
	}
	write_file(acfile, autocorrout, dimension+1, num_segments);
	deallocate_2D_array(autocorr,dimension, num_segments);
	deallocate_2D_array(autocorrout,dimension+1, num_segments);
}

/*! \brief Find autocorrelation of data at different points in the chain length and output to autocorr
 */
void auto_corr_from_data_accel(double **output, /**< Chain data input*/
				int dimension, /**< Dimension of the data*/
				int N_steps, /**< Number of steps in the data*/
				int num_segments, /**< number of segments to calculate the autocorrelation length*/
				double target_corr, /**< Target correlation ratio*/
				double **autocorr /**<[out] Autocorrelation lengths for the different segments*/
				)
{
	int device_num;
	cudaGetDeviceCount(&device_num);
	chains_internal = output;
	target_corr_internal = target_corr;
	chain_length_internal = N_steps;
	num_segments_internal = num_segments;
	dimension_internal = dimension;
	autocorr_internal = autocorr;

	GPUplan plans[device_num];
	
	for(int i = 0 ; i< device_num; i++)
	{
		plans[i].device_id = i;
		allocate_gpu_plan(&plans[i],chain_length_internal, dimension_internal, num_segments_internal);
		copy_data_to_device(&plans[i], chains_internal, chain_length_internal, dimension_internal, num_segments_internal);
	}
	plans_global = plans;

	//ThreadPoolKernelLaunch kernelpool;
	{
		threadPool<> kernelpool(device_num, ac_gpu_wrapper);
		for(int i =0; i< dimension_internal*num_segments_internal; i++){
			kernelpool.enqueue(i);
		}
	}
	//Wait for final jobs to finish before deallocating gpu memory
	for(int i = 0 ; i< device_num; i++){
		cudaSetDevice(i);
		cudaStreamSynchronize(plans_global[i].stream);
	}

	//Copy over data from Device to Host
	double **lags = allocate_2D_array(device_num, num_segments_internal*dimension_internal);
	int *lags_transfer;
	cudaMallocHost((void **)&lags_transfer, sizeof(int)* num_segments_internal*dimension_internal);
	for(int i =0; i<device_num; i++){
		cudaSetDevice(i);
		cudaMemcpyAsync(lags_transfer, plans_global[i].device_lags, 
			sizeof(int)*dimension_internal*num_segments_internal, 
			cudaMemcpyDeviceToHost,plans_global[i].stream );
		cudaStreamSynchronize(plans_global[i].stream);
		for(int j =0; j<num_segments_internal*dimension_internal; j++)
			lags[i][j] = lags_transfer[j];
	}
	cudaFreeHost(lags_transfer);
	for(int i =0 ; i<num_segments_internal*dimension_internal; i++)
	{
		for(int j = 0; j<device_num; j++){
			if(lags[j][i] != 2*chain_length_internal){
				int dim = i/num_segments_internal;
				int k = i - dim*num_segments_internal;
				autocorr[dim][k] = lags[j][i];
			}
		}
	}
	deallocate_2D_array(lags, device_num, num_segments_internal);
		
	std::cout<<"STREAMS SYNCED"<<std::endl;
	std::cout<<"DEALLOCATING"<<std::endl;
	for(int i = 0 ; i< device_num; i++)
	{
		deallocate_gpu_plan(&plans_global[i], chain_length_internal, dimension_internal, num_segments);	
	}
	
	for(int i =0 ; i< device_num; i++){
		cudaSetDevice(i);
		cudaDeviceReset();
	}
	std::cout<<std::endl;

	
}
/*! \brief Wrapper function for the thread pool
 */
void ac_gpu_wrapper(int thread, /**< Host thread*/
			int job_id/**< Job ID*/
			)
{
	launch_ac_gpu(thread, job_id, chains_internal, 
		chain_length_internal, dimension_internal, 
		target_corr_internal, num_segments_internal);
}

/*! \brief Launch the GPU kernel, formatted for the thread pool
 */
void launch_ac_gpu(int device, int element, double **data, int length, int dimension, double target_corr, int num_segments)
{
	//std::cout<<"Thread (device): "<<device<<std::endl;
	cudaSetDevice(device);
	int dim = element/num_segments;
	int k = element-dim*num_segments;
	int length_step = length / num_segments;
	int length_seg = (k+1) * length_step;
	int *host_seg;
	int start_id = dim * length;
	//plans_global[device].initial_lag = &length_seg;

	double sum = 0;
	for (int i =start_id ; i< start_id + length_seg; i++){
		sum+=plans_global[device].host_data[i];
	}
	double average = sum/length_seg;

	double var=0;
	auto_corr_internal( plans_global[device].host_data, 
				length_seg, 0, average, &var, start_id);

	//cudaMemcpy(plans_global[device].device_data, 
	//		chains_internal, 
	//		sizeof(double)*length, 
	//		cudaMemcpyHostToDevice);
	//cudaMemcpy(&plans_global[device].device_lags[element], 
	//		&length_seg,
	//		sizeof(int), cudaMemcpyHostToDevice
	//		);
	//cudaMemcpyAsync(plans_global[device].device_data, 
	//		plans_global[device].host_data, 
	//		sizeof(double)*length_seg, 
	//		cudaMemcpyHostToDevice, plans_global[device].stream);
	//cudaMemcpyAsync(&plans_global[device].device_lags[element], 
	//		plans_global[device].initial_lag,
	//		sizeof(int), cudaMemcpyHostToDevice, 
	//		plans_global[device].stream);
	
	int N = length_seg;
	auto_corr_internal_kernal
		<<<N/THREADS_PER_BLOCK,THREADS_PER_BLOCK, 
		0, plans_global[device].stream>>>
		(plans_global[device].device_data, length_seg, average, 
		&plans_global[device].device_lags[element], target_corr, var, start_id);

	//auto_corr_internal_kernal
	//	<<<N/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>
	//	(plans_global[device].device_data, length_seg, average, 
	//	&plans_global[device].device_lags[element], target_corr, var, start_id);
	//cudaMemcpy(plans_global[device].host_lag, 
	//		&plans_global[device].device_lags[element], sizeof(int), 
	//		cudaMemcpyDeviceToHost);
	cudaStreamSynchronize(plans_global[device].stream);
	//cudaDeviceSynchronize();
	//return *plans_global[device].host_lag - start_id;
	//std::cout<<element<<std::endl;
	//return 1;//*plans_global[device].host_lag ;
}

/*! \brief Allocates memory for autocorrelation--GPU structure
 */
void allocate_gpu_plan(GPUplan *plan, /**< Structure for GPU plan*/
		int data_length, /**< Length of data*/
		int dimension, /**< Dimension of the data*/
		int num_segments /**< Number of segments to calculate the autocorrelation length*/
		)
{
	cudaSetDevice(plan->device_id);
	
	cudaMalloc((void **)&plan->device_data, sizeof(double)*data_length*dimension);
	cudaMallocHost((void **)&plan->host_data, sizeof(double)*data_length*dimension);
	cudaMalloc((void **)&plan->device_lag, sizeof(int));
	cudaMalloc((void **)&plan->device_lags, sizeof(int)*dimension*num_segments);
	cudaMallocHost((void **)&plan->host_lag, sizeof(int));
	cudaMallocHost((void **)&plan->initial_lag, sizeof(int));
	cudaStreamCreate(&plan->stream);
}
/*! \brief Deallocates memory for the autocorrelation--GPU structure 
 */
void deallocate_gpu_plan(GPUplan *plan, /**< Structure for the GPU plan*/
		int data_length, /**< Length of data*/
		int dimension, /**< Dimension of the data*/
		int num_segments /**< Number of segments to calculate the autocorrelation length*/
		)
{	
	cudaSetDevice(plan->device_id);
	cudaFree(plan->device_data);
	cudaFree(plan->device_lag);
	cudaFree(plan->device_lags);
	cudaFreeHost(plan->host_lag);
	cudaFreeHost(plan->initial_lag);
	cudaFreeHost(plan->host_data);
	cudaStreamDestroy(plan->stream);
}
/*! \brief Copy data to device before starting kernels
 */
void copy_data_to_device(GPUplan *plan, /**< GPU plan*/
		double **input_data, /**<Input chain data*/
		int data_length, /**< Length of data*/
		int dimension, /**< Dimension of the data*/
		int num_segments /**< Number of segments to calculate the autocorrelation length*/
		)
{
	cudaSetDevice(plan->device_id);
	for(int i =0; i< dimension; i++){
		for(int j =0; j<data_length; j++){
			plan->host_data[i*data_length + j] = input_data[j][i];
		}
	}
	cudaMemcpyAsync(plan->device_data, plan->host_data, 
		sizeof(double)*data_length*dimension, 
		cudaMemcpyHostToDevice,plan->stream );
	int * data_lengths;
	cudaMallocHost((void **)&data_lengths, sizeof(int)*num_segments*dimension);
	for(int i =0 ;i < num_segments*dimension; i++)
		data_lengths[i] = data_length*2;
	cudaMemcpyAsync(plan->device_lags, data_lengths, 
		sizeof(int)*dimension*num_segments, 
		cudaMemcpyHostToDevice,plan->stream );
	cudaFreeHost(data_lengths);
	
}



/*! \brief Faster approximation of the autocorrelation of a chain. Implements FFT/IFFT -- GPU accelerated
 *
 * Based on the Wiener-Khinchin Theorem.
 *
 * Algorithm used from https://lingpipe-blog.com/2012/06/08/autocorrelation-fft-kiss-eigen/
 *
 */
//void auto_correlation_spectral_accel(double *chain, int length, double *autocorr)
//{
//
//	//Normalize
//	double *x_cent = (double *)malloc(sizeof(double)*length);
//	//Calculate Average
//	double ave = 0;
//	for(int i =0; i<length; i++)
//		ave+= chain[i];
//	ave /= length;
//	
//	//Create normalized vector
//	for(int i = 0 ; i<length; i++){
//		x_cent[i] = chain[i]-ave;
//	}
//
//	//Padded length
//	int L = pow(2, std::ceil( std::log2(length) ) );	
//
//	//Padded Vector
//	double *x_pad = (double *)malloc(sizeof(double)*L);
//
//	//Copy centered vector
//	for(int i = 0 ; i < length; i++){
//		x_pad[i] = x_cent[i];
//	}
//
//	//Add padding
//	for(int i = length ; i < L; i++){
//		x_pad[i] = 0;
//	}
//
//	//Allocate FFTW3 memory
//	double *norm = (double *)malloc(sizeof(double)*L);
//	std::complex<double> *out_host = (std::complex<double> *)malloc(sizeof(std::complex<double>)*L);
//	//fftw_outline plan;
//	//cufftHandle plan; 
//	cufftReal *in;
//	cufftComplex *out;
//	int batch = 1;
//	
//	cudaMalloc((void **) &in, sizeof(cfftReal)*L*batch);
//	cudaMalloc((void **) &out, sizeof(cfftComplex)*L*batch);
//	cufftPlan1d(&plan, L, CUFFT D2Z, batch); 
//	cufftExecD2Z(plan, in, data, CUFFT FORWARD);
//	cudaMemcpy(out_host,out, sizeof(std::complex<double>)*L, cudaMemcpyDeviceToHost);
//	cufftDestroy(plan); 
//	cudaFree(in);
//	cudaFree(out);
//	//initiate_likelihood_function(&plan, L);
//	//
//	//fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
//	//fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
//	//for(int i =0 ; i<L; i++){
//	//	in[i][0] = x_pad[i];
//	//	in[i][1] = 0;
//	//}
//
//	////Execute Forward Transform
//	//fftw_execute_dft(plan.p, in, out);
//
//	//Take norm^2 of the output
//	for(int i =0 ; i<L; i++){
//		double re = std::real(out_host[i]);
//		double im = std::imag(out_host[i]);
//		norm[i] = re*re + im*im;
//	}
//	//Execute Reverse Transform
//	fftw_outline plan_inv;
//	allocate_FFTW3_mem_inverse(&plan_inv, L);
//	for(int i =0 ; i<L; i++){
//		in[i][0] = norm[i];
//		in[i][1] = 0;
//	}
//	fftw_execute_dft(plan_inv.p, in, out);
//	
//	//acov is the result
//	double *acov = (double *)malloc(sizeof(double)*length);
//	for(int i =0 ; i< length; i++){
//		acov[i] = out[i][0];	
//	}
//
//	//adjust the cov
//	double *mask = (double *)malloc(sizeof(double)*L);
//	//first length elements are 1
//	for(int i = 0 ; i < length; i++){
//		mask[i]=1;
//	}
//	// last L-length elements are 0
//	for(int i = length ; i < L; i++){
//		mask[i]=0;
//	}
//	for(int i =0 ; i<L; i++){
//		in[i][0] = mask[i];
//		in[i][1] = 0;
//	}
//
//	//execute fft
//	fftw_execute_dft(plan.p, in ,out);
//	
//	//output vector -- will be trimmed to length
//	double *normadj = (double *)malloc(sizeof(double)*length);
//	//trimmed output
//	for(int i =0 ; i< length; i++){
//		normadj[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
//	}
//	
//	//adjust the cov vector
//	//for(int i =0 ; i<length ; i++){
//	//	acov[i]/=normadj[i];
//	//}
//	
//	double var = acov[0];
//
//	for(int i = 0 ; i< length; i++)
//		autocorr[i] = acov[i]/var;
//
//	//Free memory
//	deactivate_likelihood_function(&plan);
//	deactivate_likelihood_function(&plan_inv);
//	free(norm);
//	free(out_host);
//	free(mask);
//	free(normadj);
//	fftw_free(in);
//	fftw_free(out);
//	free(x_cent);
//	free(x_pad);
//}
