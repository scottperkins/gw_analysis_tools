#ifndef CUDA_UTILITIES_H
#define CUDA_UTILITIES_H
#include <string>

/*!\file
 *
 * Header file for CUDA accelerated algorithms
 *
 * Currently, no algorithms are used in any other parts of the project, so if CUDA or CUDA-enabled devices are not available, this file can be skipped in compilation by commenting out the OBJECTSCUDA line in the makefile
 */

#define THREADS_PER_BLOCK 512

static double **chains_internal;
static double target_corr_internal;
static int dimension_internal;
static int chain_length_internal;
static int num_segments_internal;
static double **autocorr_internal;

void write_file_auto_corr_from_data_file_accel(std::string acfile, std::string chains_file, int dimension, int N_steps, int num_segments, double target_corr);

void write_file_auto_corr_from_data_accel(std::string acfile, double **output, int dimension, int N_steps, int num_segments, double target_corr);

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, int num_segments, double target_corr, double **autocorr);

void launch_ac_gpu(int device, int element, double **data, int length, int dimension, double target_corr, int num_segments);

void ac_gpu_wrapper(int thread, int job_id);
#endif

