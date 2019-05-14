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
	write_data();
	return 1;
}


void write_data(){
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T= gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_uniform(r);




	//output files:
	std::string numerical_file = "testing/fisher_data/numerical.csv";
	std::string autodiff_file = "testing/fisher_data/autodiff.csv";
	std::string param_file = "testing/fisher_data/parameters.csv";
	std::string time_file = "testing/fisher_data/timing.csv";

	int num_samples = 1000;	
	int dimension = 7;
	std::string method= "IMRPhenomD";

	int length = 3000;
	double fhigh =1000;
	double flow =15;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	double ***outputauto = allocate_3D_array(num_samples,dimension,dimension);	
	double ***outputnum = allocate_3D_array(num_samples,dimension,dimension);	
	double *times_auto = (double *)malloc(sizeof(double)*num_samples);
	double *times_num = (double *)malloc(sizeof(double)*num_samples);

	double mass1s[num_samples];
	double mass2s[num_samples];
	double spin1s[num_samples];
	double spin2s[num_samples];
	double DLs[num_samples];
	double tcs[num_samples];
	double phics[num_samples];
	
	double maxm = 100., minm=3, maxspin=.9, minspin=-.9, 
		maxDL=1000, minDL=10,maxtc = 10, mintc=0, maxphic=2*M_PI, minphic=0;

	double m1temp, m2temp;
	double alpha;

	for(int i =0; i<num_samples;i++)
	{
		alpha = gsl_rng_uniform(r);
		m1temp = alpha*(maxm-minm)+minm;	

		alpha = gsl_rng_uniform(r);
		m2temp = alpha*(maxm-minm)+minm;	

		if(m1temp>m2temp){ mass1s[i] = m1temp; mass2s[i]=m2temp;}
		else{ mass2s[i] = m1temp; mass1s[i]=m2temp;}

		alpha = gsl_rng_uniform(r);
		DLs[i]= alpha*(maxDL-minDL)+minDL;

		alpha = gsl_rng_uniform(r);
		spin1s[i]= alpha*(maxspin-minspin)+minspin;

		alpha = gsl_rng_uniform(r);
		spin2s[i]= alpha*(maxspin-minspin)+minspin;

		alpha = gsl_rng_uniform(r);
		tcs[i]= alpha*(maxtc-mintc)+mintc;

		alpha = gsl_rng_uniform(r);
		phics[i]= alpha*(maxphic-minphic)+minphic;
	}

	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted",noise,length);
	for(int i =0; i<length;i++)
		noise[i] = noise[i]* noise[i];

	clock_t start7,end7;
	gen_params params;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.sky_average=true;

	//omp_set_num_threads(10);
	//#pragma omp parallel for
	for (int i =0; i< num_samples; i++)
	{
		params.mass1 = mass1s[i];//Sol Mass
		params.mass2 = mass2s[i];//Sol Mass
		params.spin1[0] = 0;
		params.spin1[1] = 0;
		params.spin1[2] = spin1s[i];
		params.spin2[0] = 0;
		params.spin2[1] = 0;
		params.spin2[2] = spin2s[i];
		params.phic = phics[i];
		params.tc = tcs[i];
		params.Luminosity_Distance = DLs[i];//MPC

		//params.betappe = new double[1] ;
		//params.betappe[0]=-100.;
		//params.bppe  =new int[1];
		//params.bppe[0] = -3;
		//params.Nmod = 1;
		

		
		start7 = clock();
		fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", 
				outputnum[i], dimension, &params,NULL,NULL, noise);

		end7 = clock();
		times_num[i]=(double)(end7-start7)/CLOCKS_PER_SEC;
		//cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
		start7 = clock();
		fisher_autodiff(freq, length, "IMRPhenomD","Hanford_O1_fitted", 
				outputauto[i], dimension, &params,NULL,NULL, noise);

		end7 = clock();
		times_auto[i]=(double)(end7-start7)/CLOCKS_PER_SEC;

		printProgress((double)i/num_samples);

		//cout<<"TIMING: FISHER autodiff: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;

		//cout.precision(5);
		//cout<<"PARAMETERS: m1, m2, spin1, spin2, DL, phic, tc"<<std::endl;	
		//cout<<params.mass1<<" , "<<params.mass2<<" , "<<params.spin1[2]<<" , "<<params.spin2[2]<<" , "<<params.Luminosity_Distance<<" , "<<params.phic<<" , "<<params.tc<<std::endl;	
		//cout<<"FISHER NUM: "<<std::endl;
		//for (int k = 0;k <dimension;k++)
		//{
		//	for (int j=0;j <dimension; j++)
		//		cout<<outputnum[i][k][j]<<"   ";
		//	cout<<endl;
		//}
		//cout<<endl;
		//cout<<"FISHER AUTO: "<<std::endl;
		//for (int k = 0;k <dimension;k++)
		//{
		//	for (int j=0;j <dimension; j++)
		//		cout<<outputauto[i][k][j]<<"   ";
		//	cout<<endl;
		//}
		//cout<<endl;
		//std::cout<<"fractional DIFF: "<<std::endl;
		//for (int k = 0;k <dimension;k++)
		//{
		//	for (int j=0;j <dimension; j++)
		//		cout<<(outputauto[i][k][j]-outputnum[i][k][j])/outputauto[i][k][j]<<"   ";
		//	cout<<endl;
		//}
	}

	double **outputnum_trans = allocate_2D_array(num_samples, dimension*dimension);
	double **outputauto_trans = allocate_2D_array(num_samples, dimension*dimension);
	double **params_out = allocate_2D_array(num_samples, dimension);
	double **time_out = allocate_2D_array(num_samples, 2);
	for (int i =0;i<num_samples; i++){
		for(int j =0; j<dimension; j++){
			for(int k =0; k<dimension; k++){
				outputnum_trans[i][j*dimension+k] = outputnum[i][j][k];
				outputauto_trans[i][j*dimension+k] = outputauto[i][j][k];
			}
		}
		params_out[i][0] =mass1s[i];
		params_out[i][1] =mass2s[i];
		params_out[i][2] =spin1s[i];
		params_out[i][3] =spin2s[i];
		params_out[i][4] =DLs[i];
		params_out[i][5] =tcs[i];
		params_out[i][6] =phics[i];

		time_out[i][0] = times_num[i];
		time_out[i][1] = times_auto[i];
	}

	write_file(numerical_file,outputnum_trans, num_samples, dimension*dimension);
	write_file(autodiff_file,outputauto_trans, num_samples, dimension*dimension);
	write_file(param_file,params_out, num_samples, dimension);
	write_file(time_file,time_out, num_samples, 2);



	deallocate_3D_array(outputauto,num_samples, dimension,dimension);
	deallocate_3D_array(outputnum,num_samples, dimension,dimension);

	deallocate_2D_array(outputauto_trans,num_samples, dimension*dimension);
	deallocate_2D_array(outputnum_trans,num_samples, dimension*dimension);
	deallocate_2D_array(params_out,num_samples, dimension);
	deallocate_2D_array(time_out,num_samples, dimension);

	gsl_rng_free(r);
	free(times_auto);
	free(times_num);
}

