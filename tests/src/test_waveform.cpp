#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/io_util.h"
#include <iostream>


#include <lal/LALSimulation.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/Sequence.h>


int time(int argc, char *argv[]);
int LALSuite_vs_GWAT_WF(int argc, char *argv[]);
void RT_ERROR_MSG();
const double MPC_M=3.08567758128e22;

int main(int argc, char *argv[])
{
	std::cout<<"TESTING WAVEFORM CALCULATIONS"<<std::endl;
		
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = stoi(argv[1]);	
	if(runtime_opt == 0){
		return time(argc,argv);
	}
	if(runtime_opt == 1){
		return LALSuite_vs_GWAT_WF(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int time(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .3;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phic = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	
	//params.mass1 = 5e6;
	//params.mass2 = 5e6;
	//params.f_ref = 20;
	//params.psi = .1;
	//double fmin = 20;
	//double fmax = 2048;
	//params.incl_angle = .1;
	//params.tc = 10;
	//params.equatorial_orientation = false;
	//double T = 32;
	params.mass1 = 9e1;
	params.mass2 = 4e1;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 0;
	params.equatorial_orientation = true;
	double fmin = 3e-2;
	double fmax = 1e-1;
	double T = T_year/12;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}

	std::string method = "IMRPhenomPv2";

	double *times_N = new double[length];
	double *times_AD = new double[length];
	double *times_0PN = new double[length];
	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	time_phase_corrected(times_N, length, frequency, &params, method, false);
	double chirpm = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	for(int i = 0 ; i<length; i++){
		times_0PN[i]=t_0PN(frequency[i],chirpm);
	}
		
	double **output = allocate_2D_array(length,2);
	for(int i = 0 ; i<length; i++){
		output[i][0] = frequency[i];
		output[i][1] = times_N[i];
	}
	write_file("data/times_N.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_AD[i];
	}
	write_file("data/times_AD.csv",output,length,2);
	for(int i = 0 ; i<length; i++){
		output[i][1] = times_0PN[i];
	}
	write_file("data/times_0PN.csv",output,length,2);
	
	deallocate_2D_array(output,length,2);
	delete [] times_N;
	delete [] times_AD;
	delete [] times_0PN;
	
	delete [] frequency;
	return 0;
}
int LALSuite_vs_GWAT_WF(int argc, char *argv[])
{
	std::cout.precision(15);
	bool P = true;
	//###############################################################################
	COMPLEX16FrequencySeries *hptilde=NULL;
	COMPLEX16FrequencySeries *hctilde=NULL;
	const REAL8 s1x = .7, s1y=.1,s1z=.2;
	const REAL8 s2x = -.2, s2y=.0,s2z=.8;
	//const REAL8 incl = M_PI/5.;
	const REAL8 incl = .8;
	REAL8 chi1_l  ;
	REAL8 chi2_l  ;
	REAL8 chip ;
	REAL8 thetaJ ;
	REAL8 zeta_polariz ;
	const REAL8 m1_SI = 7.1*LAL_MSUN_SI;
	const REAL8 m2_SI = 7*LAL_MSUN_SI;
	const REAL8 distance = 100e23;
	REAL8 alpha0 ;
	const REAL8 phiRef = 40;
	REAL8 phi_aligned;
	//double deltaf = .1;
	//const REAL8 f_min = 1e1;
	const REAL8 f_min = .002*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
	//const REAL8 f_max = 7e2;
	const REAL8 f_max = .2*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
	//int length = 4016;
	int length = 14016;
	//int length = 10;
	double deltaf = (f_max-f_min)/(length-1);
	const REAL8 f_ref = (f_max-f_min)/2.;
	IMRPhenomP_version_type  version = IMRPhenomPv2_V;
	LALDict *extraParams = NULL;

	NRTidal_version_type tidalType= NoNRT_V;

	REAL8Sequence *freqs = XLALCreateREAL8Sequence(length);
	double *frequencies = new double[length];
	for(int i = 0 ; i<length; i++){
		freqs->data[i] = f_min + i * deltaf;
		frequencies[i] = f_min+i*deltaf;
	}
	clock_t start,end;
	start = clock();
	XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJ, &alpha0, &phi_aligned, &zeta_polariz, m1_SI, m2_SI, f_ref, phiRef, incl, s1x,s1y,s1z,s2x,s2y,s2z,version);
	//XLALSimIMRPhenomP(&hptilde,&hctilde,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,deltaf,f_min,f_max,f_ref, version, extraParams);
	if(P)
		XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, version, tidalType,extraParams);
	else
		XLALSimIMRPhenomDFrequencySequence(&hptilde,freqs,phiRef,f_ref,m1_SI,m2_SI,s1z,s2z,distance, extraParams,tidalType);
	end = clock();
	std::cout<<"LAL timing: "<<(double)(end-start)/(CLOCKS_PER_SEC)<<std::endl;
	//###############################################################################
	gen_params param;
	param.mass1 = m1_SI/LAL_MSUN_SI;	
	param.mass2 = m2_SI/LAL_MSUN_SI;	
	param.Luminosity_Distance = distance/MPC_M;
	param.incl_angle = incl;
	param.phiRef = phiRef;
	param.f_ref = f_ref;
	param.spin1[0] = s1x;
	param.spin1[1] = s1y;
	param.spin1[2] = s1z;
	param.spin2[0] = s2x;
	param.spin2[1] = s2y;
	param.spin2[2] = s2z;
	param.NSflag1=false;
	param.NSflag2=false;
	param.sky_average=false;
	param.shift_time = true;
	param.shift_phase = true;
	//param.tc = -13 ;
	param.tc = .0 ;
	std::complex<double> *hplus = new std::complex<double>[length];
	std::complex<double> *hcross = new std::complex<double>[length];
	double *pplus = new double[length];
	double *pcross = new double[length];
	std::string method ;
	if(P)
		method = "IMRPhenomPv2";
	else 
		method = "IMRPhenomD";
	//param.Nmod_alpha = 1;
	//param.delta_alpha=new double[1];
	//param.delta_alpha[0]=1;
	//param.alphai=new int[1];
	//param.alphai[0]=2;
	start =clock();
	fourier_waveform(frequencies, length, hplus,hcross, method,&param);
	end =clock();
	fourier_phase(frequencies, length, pplus,pcross, method,&param);
	std::cout<<"GWAT timing: "<<(double)(end-start)/(CLOCKS_PER_SEC)<<std::endl;
	//########################################################################
	double **shplus = allocate_2D_array(length, 2);
	double **shcross = allocate_2D_array(length, 2);
	double **LALhplus = allocate_2D_array(length, 2);
	double **LALhcross = allocate_2D_array(length, 2);
	for(int i = 0 ; i<length; i++){
		double LALreal = GSL_REAL((hptilde->data->data)[i]);
		double LALimag = GSL_IMAG((hptilde->data->data)[i]);
		double sreal = std::real(hplus[i]);
		double simag = std::imag(hplus[i]);
		shplus[i][0] = sreal;
		shplus[i][1] = simag;
		LALhplus[i][0] = LALreal;
		LALhplus[i][1] = LALimag;
		if(P){
			LALreal = GSL_REAL((hctilde->data->data)[i]);
			LALimag = GSL_IMAG((hctilde->data->data)[i]);
			sreal = std::real(hcross[i]);
			simag = std::imag(hcross[i]);
			shcross[i][0] = sreal;
			shcross[i][1] = simag;
			LALhcross[i][0] = LALreal;
			LALhcross[i][1] = LALimag;
		}
	}
	write_file("data/ppgwat.csv",pplus, length);
	write_file("data/pcgwat.csv",pcross, length);

	std::string hpgwat = "data/hpgwat.csv";
	std::string hcgwat = "data/hcgwat.csv";
	std::string hpLAL = "data/hpLAL.csv";
	std::string hcLAL = "data/hcLAL.csv";
	std::string freq = "data/freqs.csv";
	write_file(hpgwat,shplus, length,2);
	if(P)
		write_file(hcgwat,shcross, length,2);
	write_file(hpLAL,LALhplus, length,2);
	if(P)
		write_file(hcLAL,LALhcross, length,2);
	write_file(freq,frequencies,length);
	deallocate_2D_array(shplus,length,2);
	deallocate_2D_array(shcross,length,2);
	deallocate_2D_array(LALhplus,length,2);
	deallocate_2D_array(LALhcross,length,2);
	delete [] hplus;
	delete [] hcross;
	delete [] pplus;
	delete [] pcross;
	XLALDestroyREAL8Sequence(freqs);

	return 1; 
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Test phase-time relationship"<<std::endl;
	std::cout<<"1 --- Compare LALSuite vs GWAT"<<std::endl;
}
