#include "autocorrelation.h"
#include "util.h"
#include "threadPool.h"
#include <iostream>

/*! Max length of array to use serial calculation*/
#define MAX_SERIAL 20000

/*! \file 
 *
 * Turns out calculating the autocorrelation is more complicated if you want to do it fast, so it gets its own file now
 *
 * First row is the starting index of that segment
 *
 * Second row is the length of that segment
 *
 * If cumulative, the ac is calculated in the following format:
 *
 * |-------|
 * 
 * |--------------|
 * 
 * |-------------------------|
 *
 * ...
 *
 * Else, the ac is calculated as :
 *
 * |-------|
 * 
 *         |------|
 * 
 *                 |--------|
 *
 * ...
 */
void write_auto_corr_file_from_data_file(std::string autocorr_filename,
	std::string datafile,
	int length, /**< length of input data*/
	int dimension, /**< dimension of data*/
	int num_segments, /**< number of segements to compute the auto-corr length*/
	double target_corr, /**< Autocorrelation for which the autocorrelation length is defined (lag of autocorrelation for which it equals the target_corr)*/
	int num_threads, /**< Total number of threads to use*/
	bool cumulative /**< Boolean to calculate the autocorrelation cumulatively*/
	)
{
	double **chains = allocate_2D_array(length, dimension);
	read_file(datafile, chains, length, dimension);
	write_auto_corr_file_from_data(autocorr_filename, chains, length, 
			dimension, num_segments, target_corr, num_threads,cumulative);	
	deallocate_2D_array(chains,length, dimension);
}
/*! \brief Writes the autocorrelation file from a data array
 *
 * First row is the starting index of that segment
 *
 * Second row is the length of that segment
 *
 * If cumulative, the ac is calculated in the following format:
 *
 * |-------|
 * 
 * |--------------|
 * 
 * |-------------------------|
 *
 * ...
 *
 * Else, the ac is calculated as :
 *
 * |-------|
 * 
 *         |------|
 * 
 *                 |--------|
 *
 * ...
 */
void write_auto_corr_file_from_data(std::string autocorr_filename,/**<Name of the file to write the autocorrelation to */
	double **data,/**< Input chains*/
	int length, /**< length of input data*/
	int dimension, /**< dimension of data*/
	int num_segments, /**< number of segements to compute the auto-corr length*/
	double target_corr, /**< Autocorrelation for which the autocorrelation length is defined (lag of autocorrelation for which it equals the target_corr)*/
	int num_threads, /**< Total number of threads to use*/
	bool cumulative /**< Boolean to calculate the autocorrelation cumulatively*/
	)
{
	double **autocorrout= allocate_2D_array(dimension+2, num_segments);//double just because thats what my internal functions are written for (write_file)
	int **autocorr= (int **)malloc(sizeof(int *)*dimension);
	for(int i =0 ; i<dimension; i++)
		autocorr[i] = (int *)malloc(sizeof(int)*num_segments);

	auto_corr_from_data(data, length, dimension, autocorr, 	
			num_segments, target_corr, num_threads,cumulative);

	int step = length/(num_segments);
	for(int i; i<num_segments; i++){
		if(cumulative){
			autocorrout[0][i] = 0;	
			autocorrout[1][i] = step*(i+1);	
		}
		else{
			autocorrout[0][i] = (i)*step;	
			autocorrout[1][i] = step;	
		}
	}

	for(int i =0 ; i<dimension; i++){
		for(int j =0 ; j<num_segments; j++){
			autocorrout[i+2][j] = autocorr[i][j];
		}
	}
	
	write_file(autocorr_filename, autocorrout, dimension+2, num_segments);

	for(int i =0 ; i<dimension; i++){
		free(autocorr[i]);
	}
	free(autocorr);
	deallocate_2D_array(autocorrout, dimension+2, num_segments);
}
				

/*! \brief Calculates the autocorrelation length for a set of data for a number of segments for each dimension -- completely host code, utilitizes FFTW3 for longer chuncks of the chains
 *
 * Takes in the data from a sampler, shape data[N_steps][dimension]
 *
 * Outputs lags that correspond to the target_corr -- shape output[dimension][num_segments]
 *
 * If cumulative, the ac is calculated in the following format:
 *
 * |-------|
 * 
 * |--------------|
 * 
 * |-------------------------|
 *
 * ...
 *
 * Else, the ac is calculated as :
 *
 * |-------|
 * 
 *         |------|
 * 
 *                 |--------|
 *
 * ...
 */
void auto_corr_from_data(double **data, /**<Input data */
			int length, /**< length of input data*/
			int dimension, /**< dimension of data*/
			int **output, /**<[out] array that stores the auto-corr lengths -- array[num_segments]*/
			int num_segments, /**< number of segements to compute the auto-corr length*/
			double target_corr, /**< Autocorrelation for which the autocorrelation length is defined (lag of autocorrelation for which it equals the target_corr)*/
			int num_threads, /**< Total number of threads to use*/
			bool cumulative /**< Boolean to calculate the autocorrelation cumulatively*/
			)
{
	//transpose data
	double **data_transpose = allocate_2D_array(dimension, length);
	for(int i =0; i<dimension; i++){
		for(int j = 0; j <length; j++){
			data_transpose[i][j] = data[j][i];
		}
	}
	int serial_threads;
	int fft_threads;
	if(num_threads > 2){
		serial_threads = num_threads/2;
		fft_threads = num_threads/2 + num_threads%2;
	}
	else{
		if(length<MAX_SERIAL){
			serial_threads = 0;
			fft_threads = 1;
		}
		else{
			serial_threads = 1;
			fft_threads = 0;
		}
	}
	
	int step = length/(num_segments);
	int lengths[num_segments];
	int startids[num_segments];
	int endids[num_segments];
	int fftw_lengths[num_segments];
	fftw_outline *plans_forward= (fftw_outline *)malloc(sizeof(fftw_outline)*num_segments);
	fftw_outline *plans_reverse= (fftw_outline *)malloc(sizeof(fftw_outline)*num_segments);
	threaded_ac_jobs_serial jobs_s[num_segments*dimension];
	threaded_ac_jobs_fft jobs_f[num_segments*dimension];
	for(int i =0 ; i<num_segments; i++){
		if(cumulative){
			lengths[i] = (i+1)*step;
			startids[i] = 0;
			endids[i] = lengths[i];
		}
		else{
			lengths[i] = step;
			startids[i] = (i)*step;
			endids[i] = (i+1)*step;
		}
	}
	{	
		{
			threadPool<threaded_ac_jobs_serial,comparator_ac_serial> serial_jobs(num_threads,threaded_ac_serial);
			for(int j =0 ; j<dimension; j++){
				for(int i =0 ; i< num_segments; i++)
				{
					if(lengths[i]<=MAX_SERIAL){
						jobs_s[j*num_segments + i].data = data_transpose;
						jobs_s[j*num_segments + i].length = &lengths[i];
						jobs_s[j*num_segments + i].target = &target_corr;		
						jobs_s[j*num_segments + i].dimension = j;		
						jobs_s[j*num_segments + i].start = &startids[i];		
						jobs_s[j*num_segments + i].end = &endids[i];
						jobs_s[j*num_segments + i].lag = &output[j][i];
						serial_jobs.enqueue(jobs_s[j*num_segments + i]);
					}
				}
			}
			for(int i =0 ; i<num_segments; i++){
				if(lengths[i]>MAX_SERIAL){
					fftw_lengths[i] = pow(2, std::ceil(std::log2(lengths[i])));	
					allocate_FFTW_mem_forward(&plans_forward[i],fftw_lengths[i]);
					allocate_FFTW_mem_reverse(&plans_reverse[i],fftw_lengths[i]);
				}
			}	
		}

		threadPool<threaded_ac_jobs_fft,comparator_ac_fft> fftw_jobs(num_threads,threaded_ac_spectral);
		for(int j =0 ; j<dimension; j++){
			for(int i =0 ; i< num_segments; i++)
			{
				if(lengths[i]>MAX_SERIAL){
					jobs_f[j*num_segments + i].data = data_transpose;
					jobs_f[j*num_segments + i].length = &lengths[i];
					jobs_f[j*num_segments + i].target = &target_corr;		
					jobs_f[j*num_segments + i].dimension = j;		
					jobs_f[j*num_segments + i].start = &startids[i];		
					jobs_f[j*num_segments + i].end =&endids[i];
					jobs_f[j*num_segments + i].lag = &output[j][i];
					jobs_f[j*num_segments + i].planforward = &plans_forward[i];
					jobs_f[j*num_segments + i].planreverse = &plans_reverse[i];
					fftw_jobs.enqueue(jobs_f[j*num_segments + i]);
				}
			}
		}
	}
	for(int i =0 ; i<num_segments; i++){
		if(lengths[i]>MAX_SERIAL){
			deallocate_FFTW_mem(&plans_forward[i]);
			deallocate_FFTW_mem(&plans_reverse[i]);
		}
	}	
	deallocate_2D_array(data_transpose,dimension, length);
}

/*! \brief Internal routine to calculate an spectral autocorrelation job
 *
 * Allows for a more efficient use of the threadPool class
 */
void threaded_ac_spectral(int thread, threaded_ac_jobs_fft job)
{
	double *ac = (double *)malloc(sizeof(double)*(*job.length));	
	auto_correlation_spectral(job.data[job.dimension], *job.length, *job.start,ac, job.planforward, job.planreverse);
	for(int i =0; i<*job.length; i++){
		if(ac[i]<*job.target){
			*job.lag = i;
			break;
		}	
	}
	free(ac);
}
/*! \brief Internal routine to calculate a serial autocorrelation job
 *
 * Allows for a more efficient use of the threadPool class
 */
void threaded_ac_serial(int thread, threaded_ac_jobs_serial job)
{
	*job.lag = auto_correlation_serial(job.data[job.dimension], *job.length, *job.start, *job.target);
}

/*! \brief Calculates the autocorrelation of a chain with the brute force method
 */
double auto_correlation_serial(double *arr, /**< input array*/
		int length, /**< Length of input array*/
		int start, /**< starting index (probably 0)*/
		double target/**< Target autocorrelation for which ``length'' is defined*/
		){
	double sum =0;
	int k;	
	double end = start + length;
	for(k =start ; k< end; k++){
		sum+= arr[k];
	}
	double ave = sum/length;
	double gamma_0_sum = 0;
	for (k=start;k<end;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	double rho = 1;
	double gamma_sum, gamma;
	int counter = 0;
	int h = 1;
	//int h = start;
	while(rho>target){	
		h++;
		gamma_sum=0;
		for(k=start;k<(end-h);k++){
			gamma_sum += (arr[k+h] - ave)*(arr[k]-ave);
		}

		gamma = gamma_sum/(length-h);
		rho = gamma/gamma_0;
	}	
	return h;

}

/*! \brief Wrapper function for convience -- assumes the data array starts at 0
 */
void auto_correlation_spectral(double *chain, int length, double *autocorr, fftw_outline *plan_forw, fftw_outline *plan_rev)
{
	auto_correlation_spectral(chain, length, 0, autocorr, plan_forw, plan_rev);
}
/*! \brief Faster approximation of the autocorrelation of a chain. Implements FFT/IFFT -- accepts FFTW plan as argument for plan reuse and multi-threaded applications
 *
 * Based on the Wiener-Khinchin Theorem.
 *
 * Algorithm used from https://lingpipe-blog.com/2012/06/08/autocorrelation-fft-kiss-eigen/
 *
 * *NOTE* the length used in initializing the fftw plans should be L = pow(2, std::ceil( std::log2(length) ) ) -- the plans are padded so the total length is a power of two
 *
 * Option to provide starting index for multi-dimension arrays in collapsed to one dimension
 *
 * length is the length of the segment to be analyzed, not necessarily the dimension of the chain
 */
void auto_correlation_spectral(double *chain, int length, int start, double *autocorr, fftw_outline *plan_forw, fftw_outline *plan_rev)
{
	int stop = start+length;
	//Normalize
	double *x_cent = (double *)malloc(sizeof(double)*length);
	//Calculate Average
	double ave = 0;
	for(int i =start; i<stop; i++)
		ave+= chain[i];
	ave /= length;
	
	//Create normalized vector
	for(int i = start ; i<stop; i++){
		x_cent[i-start] = chain[i]-ave;
	}

	//Padded length
	int L = pow(2, std::ceil( std::log2(length) ) );	

	//Padded Vector
	double *x_pad = (double *)malloc(sizeof(double)*L);

	//Copy centered vector
	for(int i = 0 ; i < length; i++){
		x_pad[i] = x_cent[i];
	}

	//Add padding
	for(int i = length ; i < L; i++){
		x_pad[i] = 0;
	}

	//Allocate FFTW3 memory
	double *norm = (double *)malloc(sizeof(double)*L);
	
	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
	for(int i =0 ; i<L; i++){
		in[i][0] = x_pad[i];
		in[i][1] = 0;
	}

	//Execute Forward Transform
	fftw_execute_dft(plan_forw->p, in, out);

	//Take norm^2 of the output
	for(int i =0 ; i<L; i++){
		norm[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
	}
	//Execute Reverse Transform
	for(int i =0 ; i<L; i++){
		in[i][0] = norm[i];
		in[i][1] = 0;
	}
	fftw_execute_dft(plan_rev->p, in, out);
	
	//acov is the result
	double *acov = (double *)malloc(sizeof(double)*length);
	for(int i =0 ; i< length; i++){
		acov[i] = out[i][0];	
	}

	//adjust the cov
	double *mask = (double *)malloc(sizeof(double)*L);
	//first length elements are 1
	for(int i = 0 ; i < length; i++){
		mask[i]=1;
	}
	// last L-length elements are 0
	for(int i = length ; i < L; i++){
		mask[i]=0;
	}
	for(int i =0 ; i<L; i++){
		in[i][0] = mask[i];
		in[i][1] = 0;
	}

	//execute fft
	fftw_execute_dft(plan_forw->p, in ,out);
	
	//output vector -- will be trimmed to length
	double *normadj = (double *)malloc(sizeof(double)*length);
	//trimmed output
	for(int i =0 ; i< length; i++){
		normadj[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
	}
	
	
	double var = acov[0];

	for(int i = 0 ; i< length; i++)
		autocorr[i] = acov[i]/var;

	//Free memory
	free(norm);
	free(mask);
	free(normadj);
	fftw_free(in);
	fftw_free(out);
	free(x_cent);
	free(x_pad);
}
/*! \brief Faster approximation of the autocorrelation of a chain. Implements FFT/IFFT 
 *
 * Based on the Wiener-Khinchin Theorem.
 *
 * Algorithm used from https://lingpipe-blog.com/2012/06/08/autocorrelation-fft-kiss-eigen/
 *
 */
void auto_correlation_spectral(double *chain, int length, double *autocorr)
{
	fftw_outline forward;
	fftw_outline reverse;
	int fftw_length = pow(2, std::ceil(std::log2(length)));	
	allocate_FFTW_mem_forward(&forward,fftw_length);
	allocate_FFTW_mem_reverse(&reverse,fftw_length);
	auto_correlation_spectral(chain, length,autocorr, &forward, &reverse);
	
	deallocate_FFTW_mem(&forward);
	deallocate_FFTW_mem(&reverse);
	////Normalize
	//double *x_cent = (double *)malloc(sizeof(double)*length);
	////Calculate Average
	//double ave = 0;
	//for(int i =0; i<length; i++)
	//	ave+= chain[i];
	//ave /= length;
	//
	////Create normalized vector
	//for(int i = 0 ; i<length; i++){
	//	x_cent[i] = chain[i]-ave;
	//}

	////Padded length
	//int L = pow(2, std::ceil( std::log2(length) ) );	

	////Padded Vector
	//double *x_pad = (double *)malloc(sizeof(double)*L);

	////Copy centered vector
	//for(int i = 0 ; i < length; i++){
	//	x_pad[i] = x_cent[i];
	//}

	////Add padding
	//for(int i = length ; i < L; i++){
	//	x_pad[i] = 0;
	//}

	////Allocate FFTW3 memory
	//double *norm = (double *)malloc(sizeof(double)*L);
	//fftw_outline plan;
	//initiate_likelihood_function(&plan, L);
	//
	//fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
	//fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
	//for(int i =0 ; i<L; i++){
	//	in[i][0] = x_pad[i];
	//	in[i][1] = 0;
	//}

	////Execute Forward Transform
	//fftw_execute_dft(plan.p, in, out);

	////Take norm^2 of the output
	//for(int i =0 ; i<L; i++){
	//	norm[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
	//}
	////Execute Reverse Transform
	//fftw_outline plan_inv;
	//allocate_FFTW3_mem_inverse(&plan_inv, L);
	//for(int i =0 ; i<L; i++){
	//	in[i][0] = norm[i];
	//	in[i][1] = 0;
	//}
	//fftw_execute_dft(plan_inv.p, in, out);
	//
	////acov is the result
	//double *acov = (double *)malloc(sizeof(double)*length);
	//for(int i =0 ; i< length; i++){
	//	acov[i] = out[i][0];	
	//}

	////adjust the cov
	//double *mask = (double *)malloc(sizeof(double)*L);
	////first length elements are 1
	//for(int i = 0 ; i < length; i++){
	//	mask[i]=1;
	//}
	//// last L-length elements are 0
	//for(int i = length ; i < L; i++){
	//	mask[i]=0;
	//}
	//for(int i =0 ; i<L; i++){
	//	in[i][0] = mask[i];
	//	in[i][1] = 0;
	//}

	////execute fft
	//fftw_execute_dft(plan.p, in ,out);
	//
	////output vector -- will be trimmed to length
	//double *normadj = (double *)malloc(sizeof(double)*length);
	////trimmed output
	//for(int i =0 ; i< length; i++){
	//	normadj[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
	//}
	//
	////adjust the cov vector
	////for(int i =0 ; i<length ; i++){
	////	acov[i]/=normadj[i];
	////}
	//
	//double var = acov[0];

	//for(int i = 0 ; i< length; i++)
	//	autocorr[i] = acov[i]/var;

	////Free memory
	//deactivate_likelihood_function(&plan);
	//deactivate_likelihood_function(&plan_inv);
	//free(norm);
	//free(mask);
	//free(normadj);
	//fftw_free(in);
	//fftw_free(out);
	//free(x_cent);
	//free(x_pad);
}


//Calculate the autocorrelation of a chain - if the chain is >100,000,
//the program will use a box-search method to help with computation time
/*! \brief OUTDATED -- numerically finds autocorrelation length -- not reliable 
 */
double auto_correlation(double *arr, int length , double tolerance){
	
	//if the chain is short enough, its easier to just calculate it in serial
	if(length<=100000){return auto_correlation_serial_old(arr,length);}
	
	double sum =0;
	int k;	
	for(k =0 ; k< length; k++){
		sum+= arr[k];
	}
	double ave = sum/length;
	double gamma_0_sum = 0;

	for (k=0;k<length;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	double step_multiplier = 1. + .01*tolerance;
	double error_tol = .01*tolerance;	/*error tolerance to find stopping point*/
	
	double rho = 1;
	int h= (int)(0.1*length);
	int direction = 1;	/*Variable to track which direction h should change each iteration*/
	double gamma_sum, gamma;
	int counter = 0;

	while(rho>.01+error_tol || rho < .01-error_tol || rho <0.){
		if(counter%1000==0)std::cout<<"Rho: "<<rho<<std::endl;
		
		if(counter%10000 == 0)h = 1.12*h;
		
		/* Pick new h based on direction */
		if (direction > 0&& h< (int)(length/step_multiplier)){
			h = (int)(step_multiplier*h);
		}
		else{
			h = (int)(h/step_multiplier);
		}
	
		/*calculate new Gamma and rho*/
		gamma_sum=0;
	
		for(k=0;k<(length-h);k++){
			gamma_sum += (arr[k+h] - ave)*(arr[k]-ave);
		}

		gamma = gamma_sum/(length-h);
		rho = gamma/gamma_0;
		/*Update direction for next iteration*/
		if(rho - (.01+error_tol)>0){
			direction = 1;
		}
		else if ((0.01 - error_tol) - rho>0){
			direction = -1;
		} 
		counter++;
	}	
	//printf("Loops Required %i, rho: %f \n",counter,rho);
	return h;
}	

/*! \brief OUTDATED Calculates the autocorrelation -- less general version
 */
double auto_correlation_serial_old(double *arr, int length  ){

	double sum =0;
	int k;	
	for(k =0 ; k< length; k++){
		sum+= arr[k];
	}
	double ave = sum/length;
	double gamma_0_sum = 0;
	for (k=0;k<length;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	double rho = 1;
	double gamma_sum, gamma;
	int counter = 0;
	int h = 1;
	while(rho>.01){	
		h++;
		rho = auto_correlation_internal(arr, length, h, ave)/gamma_0;
	}	
	return h;
}

/*! \brief OUTDATED -- Grid search method of computing the autocorrelation -- unreliable
 *
 * Hopefully more reliable than the box-search method, which can sometimes get caught in a recursive loop when the stepsize isn't tuned, but also faster than the basic linear, serial search
 */
double auto_correlation_grid_search(double *arr, /**< Input array to use for autocorrelation*/
			int length  , /**< Length of input array*/
			int box_num, /**< number of boxes to use for each iteration, default is 10*/
			int final_length, /**< number of elements per box at which the grid search ends and the serial calculation begins*/
			double target_length /**< target correlation that corresponds to the returned lag*/
			)
{
	//if array isn't long enough, just calculate serial
	if(length < final_length*2) 
		return auto_correlation_serial_old(arr, length);
	//#######################################################	
	//Zero lag variance
	double sum =0;
	int k;	
	for(k =0 ; k< length; k++){
		sum+= arr[k];
	}

	double ave = sum/length;
	double gamma_0_sum = 0;
	for (k=0;k<length;k++){
		 gamma_0_sum += (arr[k] - ave) * (arr[k] - ave);
	}
	double gamma_0 = gamma_0_sum/length;
	

	//#######################################################	
	int lag_previous = 0;
	int lag = length-1;	
	int count = 0;
	double rho;
	bool more_bins = true;
	double rho_final = 1;
	int start_final = 0;
	int stop_final = 0;
	int success_iteration = 0;
	while(lag != lag_previous && success_iteration < 5){
		count ++;
		lag_previous = lag;
		int start=1, stop=lag, loop_length, index;
		while(stop-start > final_length){
			int boundary_num = box_num +1;
			double auto_corrs[boundary_num];
			int lags[boundary_num];
			loop_length = stop-start;
			double lag_step = ( (double)(loop_length) ) / (boundary_num-1);
			for (int i =0 ; i<boundary_num; i++){
				lags[i] = start + (int)(lag_step * i );
				
			}
			for(int j =0 ; j<boundary_num ; j++){
				auto_corrs[j] = auto_correlation_internal(
							arr, length,lags[j], ave)/gamma_0;	
			}
			for (int j =0 ; j<box_num; j++){
				if(auto_corrs[j+1]<target_length){
					start = lags[j];
					//index = j;
					rho_final = auto_corrs[j];
					start_final = start;
					stop = lags[j+1];
					stop_final = stop;
					more_bins=false;
					break;
				}
			}
			if (more_bins) {
				std::cout<<"SAFETY "<<count<<" "<<box_num<<" "<<lag<< " "<<auto_corrs[boundary_num-1]<<std::endl;
				for(int j =0 ; j<boundary_num ; j++){
					std::cout<<auto_corrs[j]<<std::endl;	
				}
				//break;
				box_num +=5;
				if(final_length < 3*box_num){final_length*=2;}
					
			}
			more_bins = true;
			
		}
		//loop_length = stop-start;
		//rho = auto_corrs[index];
		rho = rho_final;
		lag = start_final;
		while (rho>target_length && lag<stop_final){
			lag ++;
			rho = auto_correlation_internal(arr, length, lag,ave)/gamma_0;
		}
		if(lag == lag_previous){
			success_iteration ++;
			box_num +=5;
			if(final_length < 3*box_num){final_length*=2;}
		}
		else{
			success_iteration = 0;
		}
	}	
	//std::cout<<rho<<std::endl;
	return lag;
}

/*! \brief Internal function to compute the auto correlation for a given lag
 *
 */
double auto_correlation_internal(double *arr, int length, int lag, double ave)
{
		double gamma_sum=0;
		for(int k=0;k<(length-lag);k++){
			gamma_sum += (arr[k+lag] - ave)*(arr[k]-ave);
		}
		return gamma_sum / (length-lag);

}

/*! \brief Function that computes the autocorrelation length on an array of data at set intervals to help determine convergence
 * 
 * outdated version -- new version uses FFTs
 */
void auto_corr_intervals_outdated(double *data, /**<Input data */
			int length, /**< length of input data*/
			double *output, /**<[out] array that stores the auto-corr lengths -- array[num_segments]*/
			int num_segments, /**< number of segements to compute the auto-corr length*/
			double accuracy /**< longer chains are computed numerically, this specifies the tolerance*/
			)
{
	double stepsize = (double)length/num_segments;
	int lengths[num_segments];
	for (int i =0; i<num_segments;i++)
		lengths[i]=(int)(stepsize*(1. + i));
	double *temp = (double *)malloc(sizeof(double)*length);		
	for(int l =0; l<num_segments; l++){
		for(int j =0; j< lengths[l]; j++){
			temp[j] = data[j];
		}
		//output[l]=auto_correlation(data,lengths[l], accuracy);
		output[l]=auto_correlation_serial_old(data,lengths[l]);
		//output[l]=auto_correlation_grid_search(data,lengths[l], 10, 100, .01);
	}
	free(temp);
	
}
/*! \brief OUTDATED -- writes autocorrelation lengths for a data array, but only with the serial method and only for a target correlation of .01
 */
void write_auto_corr_file_from_data(std::string auto_corr_filename, 
				double **output,
				int intervals, 
				int dimension, 
				int N_steps)
{
	double **ac = allocate_2D_array(dimension+1, intervals);
	//First row is the step size for the given auto-corr length
	double stepsize = (double)N_steps/intervals;
	for (int i =0 ; i<intervals; i++)
		ac[0][i] = (int)(stepsize*(1.+i));
	
	#pragma omp parallel for 
	for (int i = 0 ; i< dimension; i++)
	{
		auto_corr_intervals_outdated(output[i],N_steps, ac[i+1], intervals, 0.01);
	}	

	write_file(auto_corr_filename, ac, dimension+1, intervals);

	deallocate_2D_array(ac,dimension+1, intervals);
}

/*! \brief OUTDATED -- writes autocorrelation lengths for a data file, but only with the serial method and only for a target correlation of .01
 */
void write_auto_corr_file_from_data_file(std::string auto_corr_filename, 
				std::string output_file,
				int intervals, 
				int dimension, 
				int N_steps)
{
	double **output = allocate_2D_array(N_steps,dimension);
	read_file(output_file,output, N_steps, dimension);
	double **temp = (double **) malloc(sizeof(double*)*N_steps);
	for (int i = 0 ; i< dimension; i++){
		temp[i] = (double *)malloc(sizeof(double)*N_steps);
		for(int j =0; j< N_steps; j++){
			temp[i][j] = output[j][i];
		}
	}
	write_auto_corr_file_from_data(auto_corr_filename, temp, intervals, dimension, N_steps);
	

	deallocate_2D_array(output,N_steps,dimension);
	for (int i = 0 ; i< dimension; i++){
		free(temp[i]);
	}
	free(temp);
}
