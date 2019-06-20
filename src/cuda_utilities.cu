#include "cuda_utilities.h"
#include "cuda_utilities.hu"
#include <iostream>

__device__
void auto_corr_internal(double *arr, int length, int lag, double average, double *corr)
{
	double sum = 0;
	for(int i =0; i< (length - lag); i++){
		sum+= (arr[i+lag] - average ) * ( arr[i] - average );
	}		
	*corr = sum / (length - lag);
}
__global__
void auto_corr_internal_kernal(double *arr, int length, int *lag, double average, double *corr)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	//auto_corr_internal(arr, length, lag[id], average, &corr[id]);
	double sum = 0;
	for(int i =0; i< (length - lag[id]); i++){
		sum+= (arr[i+lag[id]] - average ) * ( arr[i] - average );
	}		
	double temp  = sum/ (length - lag[id]);
	//corr[id] = sum / (length - lag[id]);

}

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, double **autocorr)
{
	double *arr, *corr;
	int *lags;
	double average = 0;

	double *temp = (double*) malloc(sizeof(double)* N_steps);
	int *temp2 = (int*) malloc(sizeof(int)* N_steps);
	double *ac = (double*) malloc(sizeof(double)* N_steps);

	cudaMalloc( (void**)&arr, N_steps*sizeof(double) );
	cudaMalloc( (void**)&corr, N_steps*sizeof(double) );
	cudaMalloc( (void**)&lags, N_steps*sizeof(int) );

	int dim = 0;
	for(int	j = 0 ; j< N_steps; j++){
		//std::cout<<output[j][dim]<<std::endl;
		temp[j] = output[j][dim];	
		temp2[j] = j;
	}
	cudaMemcpy(arr, temp, sizeof(double)*N_steps, cudaMemcpyHostToDevice);
	cudaMemcpy(lags, temp2, sizeof(int)*N_steps, cudaMemcpyHostToDevice);
	
	int N = 750000;
	int threads_per_block = 512;
	
	auto_corr_internal_kernal<<<(int)((double)N/threads_per_block),threads_per_block>>>(arr, N_steps, lags, average, corr);
	//auto_corr_internal_kernal<<<N_steps,1>>>(arr, N_steps, lags, average, corr);

	//cudaMemcpy(ac, corr, sizeof(double)*N_steps, cudaMemcpyDeviceToHost);

	for(int i =0; i<100; i++)
		std::cout<<ac[i]<<std::endl;
	cudaFree(arr);
	cudaFree(corr);

	//int length = N_steps;
	//for(int j =0 ; j<N; j++){
	//	int id = j;
	//	double sum = 0;
	//	for(int i =0; i< (length - temp2[id]); i++){
	//		sum+= (temp[i+temp2[id]] - average ) * ( temp[i] - average );
	//	}		
	//	ac[id] = sum / (length - temp2[id]);
	//}
	free(temp);
	free(temp2);
	free(ac);
	
}
