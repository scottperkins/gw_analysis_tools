#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H
#include <string>
#include "util.h"

/*! \brief Class to contain spectral method jobs
 */
class threaded_ac_jobs_fft
{
public:
	double **data;/*! Read only -- Data to use -- full chain*/
	int *length;/*! Read only -- length of total data*/
	int *start; /*! Read only -- start index*/
	int *end;/*! Read only -- end index*/
	int dimension;/*! Read only -- dimension being analyzed*/
	fftw_outline *planforward;/*! fftw plan to use for spectral method*/
	fftw_outline *planreverse;/*! fftw plan to use for spectral method*/
	int *lag;/*! READ AND WRITE -- final lag */
	double *target;/*! Read only -- target correlation*/
};
/*! \brief Class to contain serial method jobs
 */
class threaded_ac_jobs_serial
{
public:
	double **data;/*! Read only -- Data to use -- full chain*/
	int *length;/*! Read only -- length of total data*/
	int *start; /*! Read only -- start index*/
	int *end;/*! Read only -- end index*/
	int dimension;/*! Read only -- dimension being analyzed*/
	int *lag;/*! READ AND WRITE -- final lag */
	double *target;/*! Read only -- target correlation*/
};
/*! \brief comparator to sort ac-jobs
 *
 * Starts with the longest jobs, then works down the list
 */
class comparator_ac_fft
{
public:
	bool operator()(threaded_ac_jobs_fft t, threaded_ac_jobs_fft k)
	{
		int t_length = t.end - t.start;
		int k_length = k.end - k.start;
		//return false;
		if(t_length<k_length){ return true;}
		else{ return true;}
	}	
};
/*! \brief comparator to sort ac-jobs
 *
 * Starts with the longest jobs, then works down the list
 */
class comparator_ac_serial
{
public:
	bool operator()(threaded_ac_jobs_serial t, threaded_ac_jobs_serial k)
	{
		int t_length = t.end - t.start;
		int k_length = k.end - k.start;
		//return false;
		if(t_length<k_length){ return true;}
		else{ return true;}
	}	
};
void auto_corr_from_data(double **data, int length, int dimension, int **output, int num_segments,  double accuracy, int num_threads);
void threaded_ac_spectral(int thread, threaded_ac_jobs_fft job);
void threaded_ac_serial(int thread, threaded_ac_jobs_serial job);
double auto_correlation_serial(double *arr, int length, int start, double target);

void auto_correlation_spectral(double *chain, int length, double *autocorr, fftw_outline *plan_forw, fftw_outline *plan_rev);
void auto_correlation_spectral(double *chain, int length, int start, double *autocorr, fftw_outline *plan_forw, fftw_outline *plan_rev);
void auto_correlation_spectral(double *chain, int length, double *autocorr);

double auto_correlation(double *arr, int length, double tolerance);

double auto_correlation_serial_old(double *arr, int length);

double auto_correlation_grid_search(double *arr, int length, int box_num=10, int final_length= 50, double target_length=.01);

double auto_correlation_internal(double *arr, int length, int lag, double ave);

void auto_corr_intervals_outdated(double *data, int length,double *output, int num_segments,  double accuracy);

void write_auto_corr_file_from_data(std::string autocorr_filename, double **output,int intervals, int dimension, int N_steps);

void write_auto_corr_file_from_data_file(std::string autocorr_filename, std::string output_file,int intervals, int dimension, int N_steps);
#endif
