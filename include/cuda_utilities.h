#ifndef CUDA_UTILITIES
#define CUDA_UTILITIES

void auto_corr_from_data_accel(double **output, int dimension, int N_steps, int num_segments, double target_corr, double **autocorr);


#endif

