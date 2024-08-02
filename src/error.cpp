#include <error.h>
#include <math.h>
#include <iostream>
#include "fisher.h"
#include "util.h"
#include "detector_util.h"
#include "waveform_generator.h"
#include "waveform_util.h"

using namespace std; 

/* This is a file to calculate statistical and systematic error. See fisher.cpp (or fisher.h) for a description of the function fisher_numerical (which is what we should use in these calculations). See util.h for a description of the gen_params_base and all of the different parameters it contains (this is what we will use to pass parameters throughout the code). The file tests/src/test_fishers.cpp gives a good example of how to call the fisher function and how to specify the waveform template used with the string "generation_method".  */
void calculate_systematic_error(
    double *frequency,
	std::complex<double>* h_model,
	std::complex<double>* h_true, 
	int length,/**< if 0, standard frequency range for the detector is used*/ 
	string generation_method, 
	string detector, 
	string reference_detector, 
	double *output,/**< double [dimension][dimension]*/
	int dimension, 
	gen_params_base<double> *parameters,
	int order,/**< Order of the numerical derivative (2 or 4)**/
	//double *parameters,
	double *noise
)
{
	//populate noise and frequency
	double *internal_noise = new double[length];
	if (noise)
	{
		for(int i = 0 ; i < length;i++)
		{
			internal_noise[i] = noise[i];
		}
	}
	else{
		populate_noise(frequency,detector, internal_noise,length);
		for (int i =0; i<length;i++)
		        internal_noise[i] = internal_noise[i]*internal_noise[i];	
	}

	std::cout<<__LINE__<<"psd input:"<<noise[0]<<std::endl;
	//Get Fisher (basic test if stuff is running or not)
	double **fisher = allocate_2D_array(dimension,dimension);
	double **fisher_inverse = allocate_2D_array(dimension,dimension);
	for(int i = 0 ; i<dimension; i++){
		for(int j = 0 ; j<dimension; j++){
			fisher[i][j]= 0;
			fisher_inverse[i][j]= 0;
		}
	}

	//fisher_numerical(frequency, length, generation_method, detector, reference_detector, fisher, dimension, parameters, order);
	fisher_numerical(frequency, length, generation_method, detector, reference_detector, fisher, dimension, parameters, order);
	
	std::cout<<__LINE__<<"Fisher input:"<<fisher[0][0]<<std::endl;

	for(int i = 0 ; i<dimension; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dimension; j++){
		  std::cout<<fisher[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	//Get derivative of waveform
	std::complex<double> **response_deriv = new std::complex<double>*[dimension];
	for (int i = 0 ; i<dimension; i++){
		response_deriv[i] = new std::complex<double>[length];
	}
	//std::complex<double> **dhConj = new std::complex<double>*[dimension];

	calculate_derivatives(response_deriv, 
			frequency,
			length, 
			dimension, 
			detector, 
			reference_detector, 
			generation_method,
			parameters,
			order);


	gsl_LU_matrix_invert(fisher, fisher_inverse, dimension);

	std::cout<<__LINE__<<"Fisher INVERSE input:"<<fisher_inverse[0][0]<<std::endl;
	

	//Get elements of systematic error
	output[0] = calculate_sys_err_elements(frequency, h_model, h_true, response_deriv, fisher_inverse, internal_noise, length, dimension, 0, 0);
	//double uwu = calculate_sys_err_elements(frequency, h_model, h_true, response_deriv, fisher_inverse, internal_noise, length, dimension, 0, 0);


	//Deallocation
	for (int i = 0 ; i<dimension; i++){
		//delete [] response_deriv[i];
		//delete [] dhConj[i];
	}
	//delete [] response_deriv;
	//delete [] dhConj;
	//delete[] internal_noise;
	deallocate_2D_array(fisher,dimension,dimension);
	deallocate_2D_array(fisher_inverse,dimension,dimension);
	
	
}

/**
 * TODO: Add a function that works like SysErrFunction[a_, b_] in the Mathematica notebook,
 *  calculating (and possibly integrating) systematic error for an element a and b
 */
double calculate_sys_err_elements(
	double *frequency,
	std::complex<double>* h_model,
	std::complex<double>* h_true,
	std::complex<double> **dh,
	double** fisher_inverse,
	double *psd, 
	int length,
	int dimension,
	int a,
	int b
){
	
	//Safety check for a and b
	if(a<0 || b<0 || a>=dimension || b>=dimension) return -100.0;

	std::cout<<"Fisher:"<<fisher_inverse[a][b]<<std::endl;
	std::cout<<"Noise:"<<psd[0]<<std::endl;

	double *integrand=new double[length];

	for(int i = 0; i < length; i++){
		//Calculate the systematic error integrand
		//integrand[i] = 10.0;
		integrand[i] = 
		real(
		fisher_inverse[a][b] 
	//	* (h_true[i] - h_model[i]) * std::conj(dh[b][i])
		)
		/
		psd[i]
		;
	}
	//Integrate and result results
	//return integrand[0];
	return 4*simpsons_sum(frequency[1]-frequency[0], length, integrand);
}