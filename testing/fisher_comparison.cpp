#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include <string>
#include "waveform_generator.h"
#include "IMRPhenomD.h"
#include "mcmc_routines.h"
#include "noise_util.h"
#include "util.h"
#include "waveform_util.h"
#include <adolc/adouble.h>
#include "fisher.h"
#include "ppE_IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "waveform_generator_C.h"
#include "mcmc_sampler.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_rng.h"
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/taping.h"
#include "limits"

void write_data();

int main()
{
	std::cout<<"TEST"<<std::endl;
	write_data();
	return 1;
}


void write_data(){

	//output files:
	std::string numerical_file = "fisher_data/numerical.csv";
	std::string autodiff_file = "fisher_data/autodiff.csv";
	std::string param_file = "fisher_data/parameters.csv";

	int length = 5000;
	double fhigh =300;
	double flow =15;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	cout<<"Freq spacing "<<df<<endl;
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;
	int dimension = 7;

	double **output = (double **)malloc(dimension * sizeof(**output));	
	double **output2 = (double **)malloc(dimension * sizeof(**output2));	

	for (int i = 0;i<dimension;i++){
		output[i] = (double *)malloc(dimension*sizeof(double));
		output2[i] = (double *)malloc(dimension*sizeof(double));
	}

	int num_samples = 1000;	
	double mass1s[num_samples];
	double mass2s[num_samples];
	double spin1s[num_samples];
	double spin2s[num_samples];
	double DLs[num_samples];


	gen_params params;
	double chirpm = 49.78;
	double eta =.21;
	params.mass1 = calculate_mass1(chirpm,eta);
	params.mass2 = calculate_mass2(chirpm,eta);
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	//double amp[length];
	//double phaseout[length];
	//complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .4;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = .0;
	params.tc = -.0;
	params.Luminosity_Distance = 410.;
	//params.betappe = new double[1] ;
	//params.betappe[0]=-100.;
	//params.bppe  =new int[1];
	//params.bppe[0] = -3;
	//params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;
	

	clock_t start7,end7;
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	cout.precision(5);
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}
	start7 = clock();
	fisher_autodiff(freq, length, "IMRPhenomD","Hanford_O1_fitted", output2, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER autodiff: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output2[i][j]<<"   ";
		cout<<endl;
	}
	std::cout<<"fractional DIFF: "<<std::endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<(output2[i][j]-output[i][j])/output2[i][j]<<"   ";
		cout<<endl;
	}
	free(output);
	free(output2);
}

