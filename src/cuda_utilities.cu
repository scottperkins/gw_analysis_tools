#include "cuda_utilities.h"
#include "cuda_utilities.hu"
#include "util.h"
#include <iostream>
#include <condition_variable>
#include <thread>
#include <queue>
#include <functional>
#include <mutex>
#include <unistd.h>

/*! \file
 */
GPUplan *plans_global;
class ThreadPoolKernelLaunch
{
public:
	explicit ThreadPoolKernelLaunch()
	{
		int device_num;
		cudaGetDeviceCount(&device_num);
		start(device_num);
		//start(1);
	}

	~ThreadPoolKernelLaunch()
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

	void public_stop()
	{
		stop();
	}
	int get_queue_length()
	{
		return mTasks.size();
	}
	void end_pool()
	{
		while(true)
		{
			if(mTasks.empty()){
				stop();
			}
			usleep(100);	
		}
	}
private:
	
	std::vector<std::thread> mThreads;
	std::condition_variable mEventVar;
	std::mutex mEventMutex;
	bool mStopping = false;
	std::queue<int> mTasks;

	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads; i++)
		{
			mThreads.emplace_back([=]{
				while(true)
				{
					int j;
					{
						std::unique_lock<std::mutex> lock{mEventMutex};
						mEventVar.wait(lock,[=]{return mStopping || !mTasks.empty(); });
						
						if (mStopping && mTasks.empty())
							break;	
						j = std::move(mTasks.front());
						mTasks.pop();
					}
					//std::cout<<"DEVICE: "<<i<<std::endl;
					ac_gpu_wrapper(i, j);
					
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
			mStopping = true;
		}
		
		mEventVar.notify_all();
		
		for(auto &thread: mThreads)
			thread.join();
	}
	
};

__device__ __host__
void auto_corr_internal(double *arr, int length, int lag, double average, double *corr, int start_id)
{
	double sum = 0;
	for(int i =0; i< (length - lag); i++){
		sum+= (arr[i+lag+start_id] - average ) * ( arr[i+start_id] - average );
	}		
	*corr = sum / (length - lag);
}
__global__ 
void auto_corr_internal_prep(int *lag, int length)
{
	*lag = length;
}

__global__
void auto_corr_internal_kernal(double *arr, int length,  double average, int *rho_index, double target_corr, double var, int start_id)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < *rho_index){
		double corr;
		auto_corr_internal(arr, length, id, average, &corr, start_id);
		if(corr/var<target_corr) atomicMin(rho_index, id);
	}

}

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, int num_segments, double target_corr, double **autocorr)
{
	int device_num;
	//int current_dev;
	//int succ;
	cudaGetDeviceCount(&device_num);
	//std::cout<<"NUMBER OF DEVICES "<<device_num<<std::endl;
	//cudaGetDevice(&current_dev);
	//std::cout<<"CURRENT DEV: "<<current_dev<<std::endl;
	//succ = cudaSetDevice(0);
	//std::cout<<"SUCCESS: "<<succ<<std::endl;
	//cudaGetDevice(&current_dev);
	//std::cout<<"CURRENT DEV: "<<current_dev<<std::endl;

	if(device_num==1){
	//if(true){
		int *rho_index;
		//HERE
		cudaMallocManaged( (void**)&rho_index, sizeof(int) );
		int dim ;
		int length_step = N_steps / num_segments;
		int iterations = dimension * num_segments;
		for(dim=0; dim<dimension; dim ++){
			for(int k =0 ; k<num_segments; k++){
				//std::cout<<"DIM: "<<dim<<std::endl;
				//std::cout<<"k: "<<k<<std::endl;
				//std::cout<<"LOOP: "<<dim*num_segments + k<<std::endl;
				int length_seg = (k+1) * length_step;
				int laginit = length_seg;

				double *temp = (double*) malloc(sizeof(double)* length_seg);
				double *arr;
				//HERE
				cudaMallocManaged( (void**)&arr, length_seg*sizeof(double) );

				for(int	j = 0 ; j< length_seg; j++){
					temp[j] = output[j][dim];	
				}

				double sum = 0;
				for (int i =0 ; i< length_seg; i++){
					sum+=temp[i];
				}
				double average = sum/length_seg;

				double var=0;
				auto_corr_internal( temp, length_seg, 0, average, &var, 0);

				cudaMemcpy(arr, temp, sizeof(double)*length_seg, cudaMemcpyHostToDevice);
				cudaMemcpy(rho_index, &laginit, sizeof(int), cudaMemcpyHostToDevice);
				
				int N = length_seg;
				
				//std::cout<<"LAUNCHING KERNAL"<<std::endl;
				auto_corr_internal_kernal
					<<<N/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>
					(arr, length_seg, average, rho_index, target_corr, var, 0);

				int lag ;
				cudaMemcpy(&lag, rho_index, sizeof(int), cudaMemcpyDeviceToHost);
				//cudaDeviceSynchronize();
				//std::cout<<"COPYING RESULTS"<<std::endl;
				autocorr[k][dim] = lag;
				free(temp);
				cudaFree(arr);
				printProgress((double)(dim*num_segments + k)/iterations);
				
			}
		}
		cudaFree(rho_index);
	}
	else{
		chains_internal = output;
		target_corr_internal = target_corr;
		chain_length_internal = N_steps;
		num_segments_internal = num_segments;
		dimension_internal = dimension;
		autocorr_internal = autocorr;

		//GPUplan plans_local[gpu_count];
		//plans = plans_local;
		GPUplan plans[device_num];
		
		for(int i = 0 ; i< device_num; i++)
		{
			plans[i].device_id = i;
			allocate_gpu_plan(&plans[i],chain_length_internal, dimension_internal, num_segments_internal);
			copy_data_to_device(&plans[i], chains_internal, chain_length_internal, dimension_internal, num_segments_internal);
		}
		plans_global = plans;

		//ThreadPoolKernelLaunch kernelpool(output, dimension, N_steps, autocorr, target_corr, num_segments);
		ThreadPoolKernelLaunch kernelpool;
		for(int i =0; i< dimension_internal*num_segments_internal; i++){
			kernelpool.enqueue(i);
		}
		//Wait for queue to empty before deallocating gpu memory
		while(kernelpool.get_queue_length() != 0){
			usleep(50000);
		}	
		//Wait for final jobs to finish before deallocating gpu memory
		for(int i = 0 ; i< device_num; i++){
			cudaSetDevice(i);
			cudaStreamSynchronize(plans_global[i].stream);
		}
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
				//std::cout<<lags_transfer[j]<<std::endl;
				lags[i][j] = lags_transfer[j];
		}
		cudaFreeHost(lags_transfer);
		//Wait for final jobs to finish before deallocating gpu memory
		for(int i = 0 ; i< device_num; i++){
			cudaSetDevice(i);
			cudaStreamSynchronize(plans_global[i].stream);
		}
		for(int i =0 ; i<num_segments_internal*dimension_internal; i++)
		{
			for(int j = 0; j<device_num; j++){
				if(lags[j][i] != 2*chain_length_internal){
					int dim = i/num_segments_internal;
					int k = i - dim*num_segments_internal;
					autocorr[k][dim] = lags[j][i];
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
	}
	for(int i = 0; i<dimension; i++){
		for(int j =0; j<num_segments; j++){
			std::cout<<i<<" "<<j<<" "<<autocorr[j][i]<<std::endl;
		}
	}
	for(int i =0 ; i< device_num; i++){
		cudaSetDevice(i);
		cudaDeviceReset();
	}
	std::cout<<std::endl;

	
}
void ac_gpu_wrapper(int thread, int job_id)
{
	int dim = job_id/num_segments_internal;
	int k = job_id-dim*num_segments_internal;
	autocorr_internal[k][dim] = 
		launch_ac_gpu(thread, job_id, chains_internal, 
		chain_length_internal, dimension_internal, 
		target_corr_internal, num_segments_internal);
}

int launch_ac_gpu(int device, int element, double **data, int length, int dimension, double target_corr, int num_segments)
{
	cudaSetDevice(device);
	int dim = element/num_segments;
	int k = element-dim*num_segments;
	int length_step = length / num_segments;
	int length_seg = (k+1) * length_step;
	int *host_seg;
	int start_id = dim * length;
	plans_global[device].initial_lag = &length_seg;

	double sum = 0;
	for (int i =start_id ; i< start_id + length_seg; i++){
		sum+=plans_global[device].host_data[i];
	}
	double average = sum/length_seg;

	double var=0;
	auto_corr_internal( plans_global[device].host_data, 
				length_seg, 0, average, &var, start_id);

	//cudaMemcpyAsync(plans_global[device].device_data, 
	//		plans_global[device].host_data, 
	//		sizeof(double)*length_seg, 
	//		cudaMemcpyHostToDevice, 
	//		plans_global[device].stream);
	//cudaMemcpyAsync(plans_global[device].device_lag, 
	//		plans_global[device].initial_lag,
	//		sizeof(int), cudaMemcpyHostToDevice, 
	//		plans_global[device].stream);
	//printf("%d\n", *plans_global[device].device_lag);
	//std::cout<<*plans_global[device].initial_lag<<std::endl;
	
	int N = length_seg;
	auto_corr_internal_kernal
		<<<N/THREADS_PER_BLOCK,THREADS_PER_BLOCK, 
		0, plans_global[device].stream>>>
		(plans_global[device].device_data, length_seg, average, 
		&plans_global[device].device_lags[element], target_corr, var, start_id);

	//printf("%d\n", plans_global[device].device_lag);
	//cudaMemcpy(plans_global[device].host_lag, 
	//		&plans_global[device].device_lags[element], sizeof(int), 
	//		cudaMemcpyDeviceToHost);
	//cudaStreamSynchronize(plans_global[device].stream);
	//return *plans_global[device].host_lag - start_id;
	//std::cout<<element<<std::endl;
	return 1;//*plans_global[device].host_lag ;
}

void allocate_gpu_plan(GPUplan *plan, int data_length, int dimension, int num_segments)
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
void deallocate_gpu_plan(GPUplan *plan, int data_length, int dimension, int num_segments)
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
void copy_data_to_device(GPUplan *plan, double **input_data,int data_length, int dimension, int num_segments)
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
