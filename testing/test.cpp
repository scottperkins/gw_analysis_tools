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

using namespace std;

void test_1();
void test_2();
void test_3();
void test_4();
void test_5();



int main(){

	test_1();	
	return 0;
}
void test_5()
{

}
void test_4()
{

	cout.precision(15);

	gen_params params;

	int length = 16000;
	double fhigh =20;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);
	double *freqnew = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++){
		freq[i]=flow+i*df;
		//freqnew[i] = freq[i]-.1;
		freqnew[i] = freq[i];
	}
	//for(int i=0;i<length;i++){
	//	cout<<freqnew[i]<<" "<<freq[i] <<endl;
	//}
	
	double chirpmass = 20;
	double eta = .2;
	params.mass1 = calculate_mass1(chirpmass,eta);
	params.mass2 = calculate_mass2(chirpmass,eta);
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	string method= "IMRPhenomPv2";
	complex<double> *waveformout = (complex<double> *)malloc(sizeof(complex<double>) * length);
	complex<double> *waveformoutcross = (complex<double> *)malloc(sizeof(complex<double>) * length);
	params.spin1[0] = 0.01;
	params.spin1[1] = 0;
	params.spin1[2] = .1;
	params.spin2[0] = 0.01;
	params.spin2[1] = 0;
	params.spin2[2] = .3;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 1.0;
	params.tc = .0;
	params.Luminosity_Distance = 100.2;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 1.2;
	params.theta =3.4;
	params.incl_angle = 1.3;
	//params.f_ref=10;
	//params.phiRef=0.;
	
	
	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,waveformoutcross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;

	double *noise = (double *)malloc(sizeof(double)*length);
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];
	
	//params.mass1 = 100;
	//params.mass2 = 5;
	fftw_outline plan;
	initiate_likelihood_function(&plan,length);
	
	int masslen = 10;
	//double chirp = calculate_chirpmass(params.mass1,params.mass2);
	//double eta = calculate_eta(params.mass1,params.mass2);
	double masses[masslen];
	for (int i =0; i <masslen;i++)
		masses[i] = (i+1.)*.1*chirpmass;
	double mass2 = params.mass2;
	cout<<"Mass2: "<<mass2<<std::endl;
	


	double *real = (double *)malloc(sizeof(double)*length);
	double *imag = (double *)malloc(sizeof(double)*length);
	for ( int i =0; i<length;i++)
	{
		real[i]=(waveformout[i]).real();
		imag[i]=(waveformout[i]).imag();
	}
	
	complex<double> *hplus_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *detector_response = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hcross_new = (complex<double> *)malloc(sizeof(complex<double>)*length);
	complex<double> *hplus_old = (complex<double> *)malloc(sizeof(complex<double>)*length);
	double llnew, llold,sum;
	complex<double> q;
	for (int i =0;i<masslen;i++)
	{
		params.mass1 = calculate_mass1(masses[i],eta);
		params.mass2 = calculate_mass2(masses[i],eta);

		
		start = clock();
		llnew = maximized_coal_Log_Likelihood(waveformout, noise,freqnew,length, 
					&params,"Hanford","IMRPhenomPv2",&plan);
		//llnew = maximized_coal_Log_Likelihood(real,imag, noise,freqnew,length, 
					//&params,"Hanford","IMRPhenomD",&plan);
		end = clock();
		cout<<"logl new  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		start = clock();
		fourier_detector_response(freqnew,length, detector_response,"Hanford",
					"IMRPhenomPv2",&params);
		end = clock();
		cout<<"waveform Pv2  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		start = clock();
		llold = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, real,imag, 
					noise,  calculate_chirpmass(params.mass1,params.mass2), 
					calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
		end = clock();
		
		start = clock();
		q = Q(params.theta,params.phi,params.incl_angle);
		fourier_waveform(freq,length, hplus_old,
					"IMRPhenomD",&params);
		for (int i = 0; i<length;i++)
			hplus_old[i] = q * hplus_old[i];
		end = clock();
		cout<<"waveform old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
		
		//cout<<"LOGLnew: "<<llnew<<endl;
		//cout<<"LOGLold: "<<llold<<endl;
		cout<<"chirpmass: "<<masses[i]<<" new: "<<llnew<<" old: "<<llold<<" diff ll: "<<(llold-llnew)/llold<<endl;
		for (int i =0 ; i< length;i ++)
			sum += abs((detector_response[i]-hplus_old[i])/hplus_old[i]);
		cout<<"Average Diff waveform plus: "<<sum/length<<endl;
	}



	free(hplus_new);
	free(hcross_new);
	free(hplus_old);
	free(real);
	free(imag);
	free(waveformout);
	free(waveformoutcross);
	free(freq);
	free(freqnew);
	free(detector_response);
	free(noise);
	delete [] params.betappe;
	delete [] params.bppe;

	deactivate_likelihood_function(&plan);	
}
void test_3()
{
	gen_params params;
	IMRPhenomPv2<double> modeld;
	IMRPhenomPv2<adouble> modela;
	int length = 900;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomPv2";
	//string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout_plus[length];
	complex<double> waveformout_cross[length];
	params.spin1[0] = .0;
	params.spin1[1] = .0;
	params.spin1[2] = -.2;
	params.spin2[0] = .0;
	params.spin2[1] = 0.0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=1.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = M_PI/3.;
	params.theta = M_PI/3;
	params.incl_angle = 0;
	
	double freq[length];
	for(int i=0;i<length;i++)
		freq[i]=10.+i*1e-1;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout_plus,waveformout_cross,method,&params);
	end=clock();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	delete [] params.betappe;
	delete [] params.bppe;
	
}
void test_2()
{
	double alpha, epsilon;
	IMRPhenomPv2<double> modelP;
	alpha = modelP.alpha(2,3,1./2,.75);
	epsilon = modelP.epsilon(2,3,1./2,.75);
	cout<<alpha<<endl; 
	cout<<epsilon<<endl; 
	long factorial_num = factorial(15);
	cout<<factorial_num<<endl;

	double d = modelP.d(2,1,0,.4);
	cout<<d<<endl;
}
void test_1()
{

	gen_params params;
	IMRPhenomD<double> modeld;
	IMRPhenomD<adouble> modela;
	int length = 1000;
	params.mass1 = 200;
	params.mass2 = 50;
	string method= "IMRPhenomD";
	//string method= "ppE_IMRPhenomD_Inspiral";
	double amp[length];
	double phaseout[length];
	complex<double> waveformout[length];
	params.spin1[0] = 0;
	params.spin1[1] = 0;
	params.spin1[2] = -.2;
	params.spin2[0] = 0;
	params.spin2[1] = 0;
	params.spin2[2] = .9;
	double *spin1  = params.spin1;
	double *spin2= params.spin2;
	params.phic = 2.0;
	params.tc = 8.0;
	params.Luminosity_Distance = 800.;
	params.betappe = new double[1] ;
	params.betappe[0]=10.;
	params.bppe  =new int[1];
	params.bppe[0] = -1;
	params.Nmod = 1;
	params.Nmod = 1;
	params.NSflag = false;
	params.phi = 0;
	params.theta = 0;
	params.incl_angle = 0;
	params.f_ref = 100;
	params.phiRef = 1.0;
	
	double fhigh =100;
	double flow =17;
	double df = (fhigh-flow)/(length-1);
	double *freq = (double *)malloc(sizeof(double) * length);

	cout<<"Freq spacing "<<df<<endl;

	for(int i=0;i<length;i++)
		freq[i]=flow+i*df;

	clock_t  start, end;
	start = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	end=clock();
	
	clock_t  start2, end2;
	start2 = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	end2=clock();

	clock_t  start3, end3;
	start3 = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	end3=clock();

	ofstream ampfile;
	ampfile.open("testing/data/amplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	ofstream phasefile;
	phasefile.open("testing/data/phase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++)
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
	phasefile.close();

	ofstream wavefilereal;
	wavefilereal.open("testing/data/real_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	ofstream wavefileimag;
	wavefileimag.open("testing/data/imag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	cout<<"TIMING waveform: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING amp: "<<(double)(end2-start2)/CLOCKS_PER_SEC<<endl;
	cout<<"TIMING phase: "<<(double)(end3-start3)/CLOCKS_PER_SEC<<endl;
	
	double noise[length];
	populate_noise(freq,"Hanford_O1_fitted", noise,length);
	for (int i =0; i<length;i++)
		noise[i] = noise[i]*noise[i];

	double snr = calculate_snr("Hanford_O1_fitted", waveformout, freq, length);
	cout<<"SNR: "<<snr<<endl;
	snr = data_snr_maximized_extrinsic(freq,length,waveformout,"Hanford_O1_fitted","IMRPhenomD",
			params);
	cout<<" FULL SNR: "<<snr<<endl;

	double logl;
	fftw_outline plan;
	clock_t  start4, end4;
	initiate_likelihood_function(&plan,length);
	start4 = clock();
	logl = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	end4 = clock();
	double logl2 = maximized_coal_log_likelihood_IMRPhenomD(freq, length, waveformout, 
				noise, snr, 1.2*calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false,&plan);
	//params.mass1=300;
	//params.mass2=100;
	//params.spin1[2] = .9;
	//params.spin2[2] = -.2;
	start = clock();
	double logl3 = maximized_coal_Log_Likelihood(waveformout, noise,freq,length, 
				&params,"Hanford","IMRPhenomD",&plan);
	end = clock();
	cout<<"logl new TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	start = clock();
	double logl4 = maximized_coal_log_likelihood_IMRPhenomD_Full_Param(freq, length, waveformout, 
				noise,  calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), params.spin1[2], params.spin2[2],params.Luminosity_Distance,params.theta,params.phi,params.incl_angle, false,&plan);
	end = clock();
	cout.precision(15);
	cout<<"logl old old TIMING: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<"logl old  TIMING: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	cout<<logl<<endl;
	cout<<logl2<<endl;
	cout<<"LOGLnew: "<<logl3<<endl;
	cout<<"LOGLold: "<<logl4<<endl;
	deactivate_likelihood_function(&plan);	



	double real_data[length];
	double imag_data[length];
	double loglpy;
	for (int i = 0; i<length; i++)
	{
		real_data[i] = real(waveformout[i]);
		imag_data[i] = imag(waveformout[i]);
	}
	start4 = clock();
	loglpy = maximized_coal_log_likelihood_IMRPhenomD(freq, length, real_data, imag_data, 
				noise, snr, calculate_chirpmass(params.mass1,params.mass2), 
				calculate_eta(params.mass1,params.mass2), spin1[2], spin2[2], false);
	end4 = clock();
	cout<<"logl TIMING with setup: "<<(double)(end4-start4)/CLOCKS_PER_SEC<<endl;
	cout<<loglpy<<endl;

//###################################################################################################
	
	method = "ppE_IMRPhenomD_Inspiral";
	clock_t  startppe, endppe;
	startppe = clock(); 
	fourier_waveform(freq, length, waveformout,method,&params);
	endppe=clock();
	cout<<"TIMING waveform ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;
	
	startppe = clock(); 
	fourier_amplitude(freq, length, amp,method,&params);
	endppe=clock();
	cout<<"TIMING amplitude ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;

	startppe = clock(); 
	fourier_phase(freq, length, phaseout,method,&params);
	endppe=clock();
	cout<<"TIMING phase ppE: "<<(double)(endppe-startppe)/CLOCKS_PER_SEC<<endl;

	//ofstream ampfile;
	ampfile.open("testing/data/ppeamplitude_output.csv");
	ampfile.precision(15);
	for(int i = 0;i<length;i++)
		ampfile<<freq[i]<<','<<amp[i]<<endl;
	ampfile.close();
	
	//ofstream phasefile;
	phasefile.open("testing/data/ppephase_output.csv");
	phasefile.precision(15);
	for(int i = 0;i<length;i++)
		phasefile<<freq[i]<<','<<phaseout[i]<<endl;
	phasefile.close();

	//ofstream wavefilereal;
	wavefilereal.open("testing/data/ppereal_waveform_output.csv");
	wavefilereal.precision(15);
	for(int i = 0;i<length;i++)
		wavefilereal<<freq[i]<<','<<real(waveformout[i])<<endl;
	wavefilereal.close();

	//ofstream wavefileimag;
	wavefileimag.open("testing/data/ppeimag_waveform_output.csv");
	wavefileimag.precision(15);
	for(int i = 0;i<length;i++)
		wavefileimag<<freq[i]<<','<<imag(waveformout[i])<<endl;
	wavefileimag.close();
	
	
	
	
	
	int dimension = 7;
	
	double parameters[dimension] = {params.mass1,params.mass2,params.Luminosity_Distance,spin1[2],spin2[2],params.phic,params.tc};
	double **amp_derivative = (double**) malloc(dimension * sizeof(**amp_derivative));
	for (int i = 0; i<dimension;i++)
		amp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **phase_derivative = (double**) malloc(dimension * sizeof(**phase_derivative));
	for (int i = 0; i<dimension;i++)
		phase_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double spin1vec[3] = {0,0,parameters[3]};
	double spin2vec[3] = {0,0,parameters[4]};
	source_parameters<double> source_params;
	source_params = source_params.populate_source_parameters(parameters[0],
			parameters[1],parameters[2],spin1vec,spin2vec,parameters[5],
			parameters[6]);

	lambda_parameters<double> lambda;
	modeld.assign_lambda_param(&source_params, &lambda);
	modeld.post_merger_variables(&source_params);
	source_params.f1 = 0.014/(source_params.M);
	source_params.f3 = modeld.fpeak(&source_params, &lambda);
	source_params.f1_phase = 0.018/(source_params.M);
	source_params.f2_phase = source_params.fRD/2.;
	source_params.bppe = new int[1];
	source_params.bppe[0] =-1;
	source_params.betappe = new double[1];
	source_params.betappe[0] = 10;
	source_params.Nmod=1;

	double A0 = source_params.A0; 
	double tc = source_params.tc;
	double phic = source_params.phic;
	double chirpmass = source_params.chirpmass;
	double symm = source_params.eta;
	double chi_s = source_params.chi_s;
	double chi_a = source_params.chi_a;


	IMRPhenomD<adouble> model;
	int amptapes[3] = {10,11,12};
	int phasetapes[3] = {13,14,15};
	model.amplitude_tape(&source_params, amptapes);
	model.phase_tape(&source_params, phasetapes);



	clock_t start5,end5;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modela.construct_amplitude_derivative(freq,length,dimension,amp_derivative, &source_params); 
	
	modela.construct_phase_derivative(freq,length,dimension,phase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream derivA;
	derivA.open("testing/data/deriv_amp.csv");
	derivA.precision(15);
	for(int i = 0;i<length;i++)
		derivA<<freq[i]<<','<<A0*amp_derivative[0][i]<<','<<amp_derivative[1][i]<<','<<amp_derivative[2][i]<<','<<chirpmass*amp_derivative[3][i]<<','<<symm*amp_derivative[4][i]<<','<<amp_derivative[5][i]<<','<<amp_derivative[6][i]<<endl;
	derivA.close();
	ofstream derivp;
	derivp.open("testing/data/deriv_phase.csv");
	derivp.precision(15);
	for(int i = 0;i<length;i++)
		derivp<<freq[i]<<','<<A0*phase_derivative[0][i]<<','<<phase_derivative[1][i]<<','<<phase_derivative[2][i]<<','<<chirpmass*phase_derivative[3][i]<<','<<symm*phase_derivative[4][i]<<','<<phase_derivative[5][i]<<','<<phase_derivative[6][i]<<endl;
	derivp.close();

	
	clock_t start7,end7;
	double **output = (double **)malloc(dimension * sizeof(**output));	
	for (int i = 0;i<dimension;i++)
		output[i] = (double *)malloc(dimension*sizeof(double));
	
	start7 = clock();
	fisher(freq, length, "IMRPhenomD","Hanford_O1_fitted", output, dimension, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimension;i++)
	{
		for (int j=0;j <dimension; j++)
			cout<<output[i][j]<<"   ";
		cout<<endl;
	}

	int dimensionppe = dimension +1;
	double **outputppe = (double **)malloc(dimensionppe * sizeof(**outputppe));	
	for (int i = 0;i<dimensionppe;i++)
		outputppe[i] = (double *)malloc(dimensionppe*sizeof(double));
	start7 = clock();
	fisher(freq, length, "ppE_IMRPhenomD_IMR","Hanford_O1_fitted", outputppe, dimensionppe, 
				&params );

	end7 = clock();
	cout<<"TIMING: FISHER ppE: "<<(double)(end7-start7)/CLOCKS_PER_SEC<<endl;
	for (int i = 0;i <dimensionppe;i++)
	{
		for (int j=0;j <dimensionppe; j++)
			cout<<outputppe[i][j]<<"   ";
		cout<<endl;
	}


	double **ppeamp_derivative = (double**) malloc(dimensionppe * sizeof(**ppeamp_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppeamp_derivative[i] = (double *)malloc(length * sizeof(double)); 
	double **ppephase_derivative = (double**) malloc(dimensionppe * sizeof(**ppephase_derivative));
	for (int i = 0; i<dimensionppe;i++)
		ppephase_derivative[i] = (double *)malloc(length * sizeof(double)); 


	ppE_IMRPhenomD_Inspiral<double> modelppe;
	start5=clock();
	//for (int i = 0; i<100;i++)
	modelppe.construct_amplitude_derivative(freq,length,dimensionppe,ppeamp_derivative, &source_params); 
	
	modelppe.construct_phase_derivative(freq,length,dimensionppe,ppephase_derivative, &source_params); 
	end5=clock();
	
	cout<<"TIMING: 2 grad ppE: "<<(double)(end5-start5)/CLOCKS_PER_SEC<<endl;
	ofstream ppederivA;
	ppederivA.open("testing/data/ppederiv_amp.csv");
	ppederivA.precision(15);
	for(int i = 0;i<length;i++)
		ppederivA<<freq[i]<<','<<A0*ppeamp_derivative[0][i]<<','<<ppeamp_derivative[1][i]<<','<<ppeamp_derivative[2][i]<<','<<chirpmass*ppeamp_derivative[3][i]<<','<<symm*ppeamp_derivative[4][i]<<','<<ppeamp_derivative[5][i]<<','<<ppeamp_derivative[6][i]<<','<<ppeamp_derivative[7][i]<<endl;
	ppederivA.close();
	ofstream ppederivp;
	ppederivp.open("testing/data/ppederiv_phase.csv");
	ppederivp.precision(15);
	for(int i = 0;i<length;i++)
		ppederivp<<freq[i]<<','<<A0*ppephase_derivative[0][i]<<','<<ppephase_derivative[1][i]<<','<<ppephase_derivative[2][i]<<','<<chirpmass*ppephase_derivative[3][i]<<','<<symm*ppephase_derivative[4][i]<<','<<ppephase_derivative[5][i]<<','<<ppephase_derivative[6][i]<<','<<ppephase_derivative[7][i]<<endl;
	ppederivp.close();
	
	
	
	
	for (int i =0;i<dimension;i++)
	{
		free( amp_derivative[i]);
		free( phase_derivative[i]);
		free(output[i]);
	}
	for (int i =0;i<dimensionppe;i++)
	{
		free(outputppe[i]);
		free( ppeamp_derivative[i]);
		free( ppephase_derivative[i]);
	}
	free(amp_derivative);
	free(freq);
	free(phase_derivative);
	free(ppeamp_derivative);
	free(ppephase_derivative);
	free(output);
	free(outputppe);
	delete [] params.betappe;
	delete [] params.bppe;
	delete [] source_params.betappe;
	delete [] source_params.bppe;
}

