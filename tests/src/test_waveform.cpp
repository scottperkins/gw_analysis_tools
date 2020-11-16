#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
//#include "gwat/IMRPhenomD.h"
#include "gwat/io_util.h"
#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


#define _LAL
#ifdef _LAL
	#include <lal/LALSimulation.h>
	#include <lal/LALDatatypes.h>
	#include <lal/LALSimIMR.h>
	#include <lal/LALConstants.h>
	#include <lal/FrequencySeries.h>
	#include <lal/LALAtomicDatatypes.h>
	#include <lal/Sequence.h>
	#include <lal/LALDetectors.h>
	#include <lal/LALDatatypes.h>
	#include <lal/DetResponse.h>
	#include <lal/Units.h>
#endif
#ifndef _LAL
	#include "gwat/gIMRPhenomD.h"
	#include "gwat/gIMRPhenomP.h"
#endif

int LALSuite_vs_GWAT_WF(int argc, char *argv[]);
int PNSeries_testing(int argc, char *argv[]);
int tc_comparison(int argc, char *argv[]);
int gIMR_testing(int argc, char *argv[]);
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
		#ifdef _LAL
			return LALSuite_vs_GWAT_WF(argc,argv);
		#endif
	}
	if(runtime_opt == 1){
		#ifdef _LAL
			return tc_comparison(argc,argv);
		#endif
	}
	if(runtime_opt == 2){
		#ifndef _LAL
			return gIMR_testing(argc,argv);
		#endif
	}
	if(runtime_opt == 3){
			return PNSeries_testing(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int PNSeries_testing(int argc, char *argv[])
{
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .1;
	params.spin1[1] = .3;
	params.spin2[1] = .1;
	params.spin1[0] = .3;
	params.spin2[0] = .1;
	//params.chip = .07;
	//params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 20;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 3.6;
	params.mass2 = 2.8;
	params.tc = 0;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;
	params.gmst=3;

	params.Nmod = 8;
	params.bppe = new double[8];
	params.bppe[0] = -7;
	params.bppe[1] = -5;
	params.bppe[2] = -4;
	params.bppe[3] = -3;
	params.bppe[4] = -2;
	params.bppe[5] = -1;
	params.bppe[6] = 1;
	params.bppe[7] = 2;
	params.betappe = new double[8];
	params.betappe[0] = .001;
	params.betappe[1] = 5;
	params.betappe[2] = -7;
	params.betappe[3] = 8;
	params.betappe[4] = -2;
	params.betappe[5] = 8;
	params.betappe[6] = 3;
	params.betappe[7] = -15;
	
	int length = 10000;
	double freqs[length];
	double fhigh =10* pow(6,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	double flow = pow(100,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	std::cout<<flow<< " " <<fhigh<<std::endl;
	double delta_f = (fhigh - flow)/length;
	for(int i = 0 ; i<length; i++){
		freqs[i] = flow + delta_f*i;
	}
	std::complex<double> response[length];
	double **output = new double*[2];
	output[0] = new double[length];
	output[1] = new double[length];


	fourier_detector_response(freqs, length, response, "Hanford", "ppE_IMRPhenomPv2_Inspiral",&params, (double*)NULL);
	for(int i =0 ; i<length ; i++){
		output[0][i] = std::real(response[i]);
		output[1][i] = std::imag(response[i]);
	}
	write_file("data/PNSeries_Raw_ppE.csv",output, 2,length);
	//################################
	for(int i = 1 ; i<8; i++){
		params.betappe[i] = params.betappe[i]/params.betappe[0];
		std::cout<<params.betappe[i]<<std::endl;
	}
	fourier_detector_response(freqs, length, response, "Hanford", "PNSeries_ppE_IMRPhenomPv2_Inspiral",&params, (double*)NULL);
	for(int i =0 ; i<length ; i++){
		output[0][i] = std::real(response[i]);
		output[1][i] = std::imag(response[i]);
	}
	write_file("data/PNSeries_PNS_ppE.csv",output, 2,length);
	
	
	

	delete [] params.bppe; delete[] params.betappe;delete [] output[0];delete [] output[1];delete [] output;
	return 0;
}
#ifndef _LAL
int gIMR_testing(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .1;
	params.spin1[1] = .3;
	params.spin2[1] = .1;
	params.spin1[0] = .3;
	params.spin2[0] = .1;
	//params.chip = .07;
	//params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 20;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 36;
	params.mass2 = 28;
	params.tc = 0;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;

	params.Nmod_phi = 1;
	params.phii = new int[1];
	params.phii[0]=4;
	params.delta_phi = new double[1];
	params.delta_phi[0]=-2;

	double beta;
	int b;
	gIMRPhenomD<double> model;
	model.ppE_gIMR_mapping(&params, 4, &beta, &b);
	std::cout<<"Beta: "<<beta<<std::endl;
	std::cout<<"b: "<<b<<std::endl;
	params.Nmod = 1;
	params.bppe = new int[1];
	params.bppe[0] = b;
	params.betappe = new double[1];
	params.betappe[0] = beta;


	double fmin = 10;
	double fmax = 80;
	double T = 4;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}
	std::complex<double> hpg[length];
	std::complex<double> hcg[length];
	std::complex<double> hpppE[length];
	std::complex<double> hcppE[length];
	fourier_waveform(frequency, length, hpg, hcg, "gIMRPhenomPv2",&params);
	fourier_waveform(frequency, length, hpppE, hcppE, "ppE_IMRPhenomPv2_Inspiral",&params);
	//fourier_waveform(frequency, length, hpg, hcg, "gIMRPhenomD",&params);
	//fourier_waveform(frequency, length, hpppE, hcppE, "ppE_IMRPhenomD_Inspiral",&params);
	for(int i = 0 ; i<length; i++){
		//std::cout<<hpg[i]<<std::endl;;
		std::cout<<( hpg[i] - hpppE[i])*2./( abs(hpg[i]) + abs(hpppE[i]))<<std::endl;;
	}
	delete [] params.delta_phi;
	delete [] params.phii;
	delete [] params.bppe;
	delete [] params.betappe;
	delete [] frequency;
	return 0;
}
#endif
#ifdef _LAL
int tc_comparison(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .3;
	params.spin2[2] = .1;
	params.chip = .7;
	params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 1e-5;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 9e6;
	params.mass2 = 4e5;
	params.theta_l = 1;
	params.phi_l = 2;
	params.tc = 0;
	params.equatorial_orientation = true;
	double fmin = 3e-5;
	double fmax = 2e-3;
	double T = T_year;

	int length = T*(fmax-fmin);
	double *frequency = new double[length];
	
	for(int i = 0 ; i<length; i++){
		frequency[i]=fmin + (double)i /T;
	}

	double fpeak,frd,fdamp;
	std::string method = "IMRPhenomD";

	postmerger_params(&params,method, &fpeak,&fdamp,&frd);
	std::cout<<"Peak frequency: "<<fpeak<<std::endl;

	double *times_AD = new double[length];
	std::complex<double> *response = new std::complex<double>[length];
	double **output = allocate_2D_array(length,4);

	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	fourier_detector_response(frequency, length,response, "LISA",method, &params, times_AD);
	for(int i = 0 ; i<length; i++){
		output[i][0] = frequency[i];
		output[i][1] = real(response[i]);
		output[i][2] = imag(response[i]);
		output[i][3] = times_AD[i];
	}
	write_file("data/response_0tc.csv",output,length,4);

	params.tc = T*3./4.;
	//params.tc = T;
	//params.tc = T_year;
	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	fourier_detector_response(frequency, length,response, "LISA",method, &params, times_AD);
	for(int i = 0 ; i<length; i++){
		output[i][0] = frequency[i];
		output[i][1] = real(response[i]);
		output[i][2] = imag(response[i]);
		output[i][3] = times_AD[i];
	}
	write_file("data/response_nonzero_tc.csv",output,length,4);

	
	params.tc = T;
	time_phase_corrected_autodiff(times_AD, length, frequency, &params, method, false,NULL);
	fourier_detector_response(frequency, length,response, "LISA",method, &params, times_AD);
	for(int i = 0 ; i<length; i++){
		output[i][0] = frequency[i];
		output[i][1] = real(response[i]);
		output[i][2] = imag(response[i]);
		output[i][3] = times_AD[i];
	}
	write_file("data/response_cyclic_tc.csv",output,length,4);
	deallocate_2D_array(output,length,4);

	delete [] times_AD;
	delete [] response;
	
	delete [] frequency;
	return 0;
}
int LALSuite_vs_GWAT_WF(int argc, char *argv[])
{
	LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO;	
	std::cout.precision(15);
	bool P = true;
	gsl_rng_env_setup();	
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	gsl_rng_set(r,10);
	int iterations = 50;
	double times[iterations][2];
	//###############################################################################
	for(int k = 0 ; k<iterations ; k++){
		int d  ;
		if(k%3 == 0 ) {d = 2;}
		else if(k%3 == 1 ) {d = 5;}
		else if(k%3 == 2 ) {d = 6;}
		LALDetector LALD = lalCachedDetectors[d];
		std::string DETECTOR = "";
		if(k%3 == 0){ DETECTOR = "Virgo";}
		else if(k%3 == 1){ DETECTOR = "Hanford";}
		else if(k%3 == 2){ DETECTOR = "Livingston";}
		COMPLEX16FrequencySeries *hptilde=NULL;
		COMPLEX16FrequencySeries *hctilde=NULL;
		double alpha[15];
		for (int j = 0 ; j<15; j++){
			alpha[j] = gsl_rng_uniform(r);
		}
		const REAL8 s1x = -.1+alpha[0]*.2, s1y=-.2+alpha[1]*.3,s1z=-.4+alpha[2]*.6;
		const REAL8 s2x = -.1+alpha[3]*.2, s2y=-.2+alpha[4]*.3,s2z=-.4+alpha[5]*.6;
		//const REAL8 incl = M_PI/5.;
		const REAL8 incl = M_PI * alpha[6];
		double RA = 2*M_PI * alpha[7];
		double DEC =-M_PI/2.+ M_PI * alpha[8];
		double psi =2*M_PI * alpha[9];
		REAL8 chi1_l  ;
		REAL8 chi2_l  ;
		REAL8 chip ;
		REAL8 thetaJ ;
		REAL8 zeta_polariz ;
		double tempm1=1+5*alpha[10],tempm2 = 1+5*alpha[11];
			
		REAL8 m1_SI;
		REAL8 m2_SI;
		if(tempm1>tempm2){
			 m1_SI= tempm1*LAL_MSUN_SI;
			 m2_SI= tempm2*LAL_MSUN_SI;
		}
		else{
			 m1_SI= tempm2*LAL_MSUN_SI;
			 m2_SI= tempm1*LAL_MSUN_SI;
		}
		const REAL8 distance = 100e23 + 100e24*alpha[12] ;
		REAL8 alpha0 ;
		const REAL8 phiRef = 2*M_PI * alpha[13];
		double gmst = 2*M_PI*alpha[14];
		REAL8 phi_aligned;
		const REAL8 f_min = .002*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		const REAL8 f_max = .19*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		int length = 4016;
		//int length = 131072;
		double deltaf = (f_max-f_min)/(length-1);
		const REAL8 f_ref = (f_max-f_min)/2.;
		IMRPhenomP_version_type  version = IMRPhenomPv2_V;
		LALDict *extraParams = NULL;

		//NRTidal_version_type tidalType= NoNRT_V;

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
		if(P){
			//XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, version, tidalType,extraParams);
			XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, version,extraParams);
			for(int i = 0 ; i<length; i++){
				gsl_complex tempPlus = (hptilde->data->data)[i];	
				gsl_complex tempCross = (hctilde->data->data)[i];	
				gsl_complex f1 = gsl_complex_rect(cos(2.*zeta_polariz),0.);
				gsl_complex f2 = gsl_complex_rect(sin(2.*zeta_polariz),0.);
				gsl_complex f3 = gsl_complex_rect(-sin(2.*zeta_polariz),0.);
				(hptilde->data->data)[i] = gsl_complex_add(gsl_complex_mul(f1,tempPlus)
						,gsl_complex_mul(f2,tempCross));
				(hctilde->data->data)[i] = gsl_complex_add(gsl_complex_mul(gsl_complex_rect(cos(2.*zeta_polariz),0.),tempCross)
						,gsl_complex_mul(f3,tempPlus));

			}
		}
		else{
			//XLALSimIMRPhenomDFrequencySequence(&hptilde,freqs,phiRef,f_ref,m1_SI,m2_SI,s1z,s2z,distance, extraParams,tidalType);
			XLALSimIMRPhenomDFrequencySequence(&hptilde,freqs,phiRef,f_ref,m1_SI,m2_SI,s1z,s2z,distance, extraParams);
			hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, freqs->data[1]-freqs->data[0], &lalStrainUnit, length);
			for(int i = 0 ; i<length; i++){
				gsl_complex f2 = gsl_complex_rect(0.,-cos(incl));
				(hctilde->data->data)[i] = gsl_complex_mul(f2,(hptilde->data->data)[i]);
				gsl_complex f1 = gsl_complex_rect(0.5*(1. + cos(incl)*cos(incl)),0);
				(hptilde->data->data)[i] = gsl_complex_mul(f1,(hptilde->data->data)[i]);

			}
		}
		COMPLEX16FrequencySeries *det = XLALCreateCOMPLEX16FrequencySeries("det: FD waveform", &ligotimegps_zero, 0.0, freqs->data[1]-freqs->data[0], &lalStrainUnit, length);
		double fplus,fcross;
		double fplusG,fcrossG;
		XLALComputeDetAMResponse(&fplus,&fcross, LALD.response, RA,DEC,psi,gmst);
		detector_response_functions_equatorial(DETECTOR,RA,DEC,psi,gmst, &fplusG,&fcrossG);
		std::cout<<"Fractional error on F+/Fx: "<<(fplus-fplusG)/fplus<<" "<<(fcross-fcrossG)/fcross<<std::endl;
		for(int i = 0 ; i<length ; i++){
			(det->data->data)[i]=gsl_complex_add(
			gsl_complex_mul(gsl_complex_rect(fplus,0.),(hptilde->data->data)[i]),
			gsl_complex_mul(gsl_complex_rect(fcross,0.0) , (hctilde->data->data)[i]));
		}
		end = clock();
		times[k][0] = (double)(end-start)/(CLOCKS_PER_SEC);
		//std::cout<<"LAL timing: "<<(double)(end-start)/(CLOCKS_PER_SEC)<<std::endl;
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
		param.equatorial_orientation=false;
		param.horizon_coord=false;
		param.RA = RA;
		param.DEC = DEC;
		param.psi = psi;
		param.gmst = gmst;
		param.tc = .0 ;
		std::complex<double> *response = new std::complex<double>[length];
		std::string method ;
		if(P)
			method = "IMRPhenomPv2";
		else 
			method = "IMRPhenomD";
		start =clock();

		fourier_detector_response(frequencies, length, response,DETECTOR, method,&param,(double*)NULL);
		end =clock();
		times[k][1] = (double)(end-start)/(CLOCKS_PER_SEC);
		//std::cout<<"GWAT timing: "<<(double)(end-start)/(CLOCKS_PER_SEC)<<std::endl;
		//########################################################################
		double *phaseLALprep = new double[length];
		double *phaseGWATprep = new double[length];
		double *phaseLAL = new double[length];
		double *phaseGWAT = new double[length];
		for(int i = 0 ; i<length ; i++){
			phaseLALprep[i]= std::atan2(GSL_IMAG((det->data->data)[i]),GSL_REAL((det->data->data)[i]));
			phaseGWATprep[i]= std::atan2(std::imag(response[i]),std::real(response[i]));
		}
		unwrap_array(phaseLALprep, phaseLAL,length);
		unwrap_array(phaseGWATprep, phaseGWAT,length);
		double **output = allocate_2D_array(length,7);
		for(int i = 0 ; i<length; i++){
			output[i][0] = frequencies[i];
			output[i][1] = GSL_REAL((det->data->data)[i]);
			output[i][2] = GSL_IMAG((det->data->data)[i]);
			output[i][3] = std::real(response[i]);
			output[i][4] = std::imag(response[i]);
			output[i][5] = phaseLAL[i];
			output[i][6] = phaseGWAT[i];
		}
		write_file("data/response_"+std::to_string(k)+".csv",output,length,7);
		deallocate_2D_array(output,length,7);
		delete [] response;
		delete [] phaseLALprep;
		delete [] phaseGWATprep;
		delete [] phaseLAL;
		delete [] phaseGWAT;
		delete [] frequencies;
		XLALDestroyREAL8Sequence(freqs);
		XLALDestroyCOMPLEX16FrequencySeries(hptilde);
		XLALDestroyCOMPLEX16FrequencySeries(hctilde);
		XLALDestroyCOMPLEX16FrequencySeries(det);
	}
	double lal_sum = 0 ;
	double gwat_sum = 0 ;
	for(int i = 0 ; i<iterations; i++){
		lal_sum +=times[i][0];
		gwat_sum +=times[i][1];

	}
	std::cout<<"Average LAL time: "<<lal_sum / iterations<<std::endl;
	std::cout<<"Average GWAT time: "<<gwat_sum / iterations<<std::endl;
	gsl_rng_free(r);
	return 1; 
}
#endif
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Compare LALSuite vs GWAT"<<std::endl;
	std::cout<<"1 --- Compare the effect of tc on LISA"<<std::endl;
	std::cout<<"2 --- test gIMR waveforms"<<std::endl;
	std::cout<<"2 --- test PNSeries waveforms"<<std::endl;
}
