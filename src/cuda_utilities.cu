#include "cuda_utilities.h"
#include "cuda_utilities.hu"
#include "util.h"
#include <iostream>

__device__ __host__
void auto_corr_internal(double *arr, int length, int lag, double average, double *corr)
{
	double sum = 0;
	for(int i =0; i< (length - lag); i++){
		sum+= (arr[i+lag] - average ) * ( arr[i] - average );
	}		
	*corr = sum / (length - lag);
}
__global__
void auto_corr_internal_kernal(double *arr, int length,  double average, int *rho_index, double target_corr, double var)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < *rho_index){
		double corr;
		auto_corr_internal(arr, length, id, average, &corr);
		if(corr/var<target_corr) atomicMin(rho_index, id);
	}

}

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, int num_segments, double target_corr, double **autocorr)
{

	int *rho_index;

	cudaMalloc( (void**)&rho_index, sizeof(int) );

	//double target_corr = .01;

	int dim ;
	int length_step = N_steps / num_segments;
	int threads_per_block = 512;
	int iterations = dimension * num_segments;

	for(dim=0; dim<dimension; dim ++){
		for(int k =0 ; k<num_segments; k++){
			int length_seg = (k+1) * length_step;
			int laginit = length_seg;

			double *temp = (double*) malloc(sizeof(double)* length_seg);
			double *arr;
			cudaMalloc( (void**)&arr, length_seg*sizeof(double) );

			for(int	j = 0 ; j< length_seg; j++){
				temp[j] = output[j][dim];	
			}

			double sum = 0;
			for (int i =0 ; i< length_seg; i++){
				sum+=temp[i];
			}
			double average = sum/length_seg;

			double var=0;
			auto_corr_internal( temp, length_seg, 0, average, &var);

			cudaMemcpy(arr, temp, sizeof(double)*length_seg, cudaMemcpyHostToDevice);
			cudaMemcpy(rho_index, &laginit, sizeof(int), cudaMemcpyHostToDevice);
			
			int N = length_seg;
			
			auto_corr_internal_kernal
				<<<N/threads_per_block,threads_per_block>>>
				(arr, length_seg, average, rho_index, target_corr, var);

			int lag ;
			cudaMemcpy(&lag, rho_index, sizeof(int), cudaMemcpyDeviceToHost);
			autocorr[k][dim] = lag;
			free(temp);
			cudaFree(arr);
			printProgress((double)(dim*num_segments + k)/iterations);
			
		}
	}
	std::cout<<std::endl;
	cudaFree(rho_index);

	
}
