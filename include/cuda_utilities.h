#ifndef CUDA_UTILITIES_H
#define CUDA_UTILITIES_H

/*!\file
 */

#define THREADS_PER_BLOCK 512

static double **chains_internal;
static double target_corr_internal;
static int dimension_internal;
static int chain_length_internal;
static int num_segments_internal;
static double **autocorr_internal;

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, int num_segments, double target_corr, double **autocorr);

int launch_ac_gpu(int device, int element, double **data, int length, int dimension, double target_corr, int num_segments);

void ac_gpu_wrapper(int thread, int job_id);
#endif

