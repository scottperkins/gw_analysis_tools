#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/pn_waveform_util.h"
#include "gwat/ppE_utilities.h"
//#include "gwat/waveform_generator.h"
//#include "gwat/IMRPhenomD.h"
#include "gwat/io_util.h"
#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>


//#define _LAL
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
	#include "gwat/TaylorT2.h"
#endif

int time_domain_testing(int argc, char *argv[]);
int LALSuite_vs_GWAT_WF(int argc, char *argv[]);
int PNSeries_testing(int argc, char *argv[]);
int tc_comparison(int argc, char *argv[]);
int gIMR_testing(int argc, char *argv[]);
int polarization_testing(int argc, char *argv[]);
int BHEvaporation_test(int argc, char *argv[]);
int EA_fully_restricted_test(int argc, char *argv[]);
int EA_fully_restricted_parameterization_test(int argc, char *argv[]);
int EA_fully_restricted_consistency_test(int argc, char *argv[]);
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
	if(runtime_opt == 4){
		return polarization_testing(argc,argv);
	}
	if(runtime_opt == 5){
		return BHEvaporation_test(argc,argv);
	}
	if(runtime_opt == 6){
		return time_domain_testing(argc,argv);
	}
	if(runtime_opt == 7){
		return EA_fully_restricted_test(argc,argv);
	}
	if(runtime_opt == 8){
		return EA_fully_restricted_parameterization_test(argc,argv);
	}
	if(runtime_opt == 9){
		return EA_fully_restricted_consistency_test(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}

int EA_fully_restricted_consistency_test(int argc, char *argv[])
{
	std::cout<<"EA CONSISTENCY TEST"<<std::endl;
	gen_params params;	
	params.spin1[1] = .0;
	params.spin2[1] = .0;
	params.spin1[0] = .0;
	params.spin2[0] = .0;
	//params.chip = .07;
	//params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 20;
	params.NSflag1 = true;
	params.NSflag2 = true;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.tc = 6;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;
	params.gmst=3;

	params.Nmod = 4;
	params.bppe = new double[4];
	//These don't matter, don't worry about them -- overwritten by prep_source_parameters
	params.bppe[0] = -13;
	params.bppe[1] = -13;
	params.bppe[2] = -13;
	params.bppe[3] = -13;
	params.betappe = new double[4];
	//params.betappe[0] = .001;

	source_parameters<double> sp ;


	int iterations = 1;
	int samples = 8032;
	double **output = allocate_2D_array(samples, 6);
		

	double FMIN = 5;
	//double FMAX = 2048;
	double FMAX = 100;
	double deltaf = (FMAX-FMIN)/samples;

	double *freqs= new double[samples];
	for (int i = 0 ; i<samples; i++){
		freqs[i] = FMIN + deltaf*i;
	}

	const gsl_rng_type *T;
	gsl_rng *r ;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (int i = 0 ; i<iterations; i++){
		//Make these random numbers
		params.betappe[0] = 1e-1; // c1
		params.betappe[1] = 2e-1; // c2
		params.betappe[2] = 3e-1; // c3
		params.betappe[3] = 4e-1; // c4
		params.mass1 = gsl_rng_uniform(r) +1;
		params.mass2 = gsl_rng_uniform(r) +1;
		if(params.mass2>params.mass1){
			double temp = params.mass2;
			params.mass2 = params.mass1;
			params.mass1 = temp;
		}

		params.spin1[2] = gsl_rng_uniform(r)*.05 -.025;
		params.spin2[2] = gsl_rng_uniform(r)*.05 -.025;

		params.tidal1 = gsl_rng_uniform(r)*100+5;
		params.tidal2 = gsl_rng_uniform(r)*100+5;
		
		std::complex<double> *responseEA =  new std::complex<double>[samples];
		std::complex<double> *responseGR =  new std::complex<double>[samples];

		fourier_detector_response(freqs, samples, responseEA, "Hanford", "EA_fully_restricted_v1_IMRPhenomD_NRT", &params, (double *) NULL);
		fourier_detector_response(freqs, samples, responseGR, "Hanford", "IMRPhenomD_NRT", &params, (double *) NULL);

		double *phase_EA = new double[samples];
		double *phase_GR = new double[samples];
		double *phase_EA_unwrap = new double[samples];
		double *phase_GR_unwrap = new double[samples];
		for(int i = 0 ; i<samples ; i++){
			phase_EA[i]= std::atan2(std::imag(responseEA[i]),std::real(responseEA[i]));
			phase_GR[i]= std::atan2(std::imag(responseGR[i]),std::real(responseGR[i]));
		}
		unwrap_array(phase_EA, phase_EA_unwrap,samples);
		unwrap_array(phase_GR, phase_GR_unwrap,samples);
	
		for(int j = 0 ; j < samples ; j++){
			output[j][0] = std::real(responseEA[j]);
			output[j][1] = std::imag(responseEA[j]);
			output[j][2] = std::real(responseGR[j]);
			output[j][3] = std::imag(responseGR[j]);
			output[j][4] = phase_EA_unwrap[j];
			output[j][5] = phase_GR_unwrap[j];
		}
		//write_file("data/EA_GR_COMP_"+std::to_string(i)+".csv", output, samples, 6);

		delete [] responseEA;
		delete [] responseGR;
		delete [] phase_EA;
		delete [] phase_GR;
		delete [] phase_EA_unwrap;
		delete [] phase_GR_unwrap;
	}
	gsl_rng_free(r);
	delete [] freqs;
	deallocate_2D_array(output, samples, 4);
	delete [] params.betappe;
	delete [] params.bppe;
	return 0;
}
int EA_fully_restricted_parameterization_test(int argc, char *argv[])
{
	gen_params params;	
	params.spin1[2] = .03;
	params.spin2[2] = .01;
	params.spin1[1] = .0;
	params.spin2[1] = .0;
	params.spin1[0] = .0;
	params.spin2[0] = .0;
	//params.chip = .07;
	//params.phip = 0.1;
	params.Luminosity_Distance = 100;
	params.phiRef = 1;
	params.RA = 2.;
	params.DEC = -1.1;
	params.f_ref = 20;
	params.NSflag1 = true;
	params.NSflag2 = true;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 1.4;
	params.mass2 = 1.3;
	params.tc = 6;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;
	params.gmst=3;

	params.Nmod = 4;
	params.bppe = new double[4];
	//These don't matter, don't worry about them -- overwritten by prep_source_parameters
	params.bppe[0] = -13;
	params.bppe[1] = -13;
	params.bppe[2] = -13;
	params.bppe[3] = -13;
	params.betappe = new double[4];
	//params.betappe[0] = .001;
	params.tidal1 = 10;
	params.tidal2 = 10;

	source_parameters<double> sp ;

	//Random number declaration and seeding
	const gsl_rng_type *T;
	gsl_rng *r; 
	
	gsl_rng_env_setup();
	
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, time(NULL)); //seeding the random number generator with time
	
	int iterations = 1000;
	int num_bad_points = 0;
	int num_infspeeds = 0;
	int dim = 13;// Increase for every factor you want to output
	double **output = allocate_2D_array(iterations, dim);

	double a_bound = 4.; 
	double theta_bound = 3.; 
	double w_bound = 2.;
	double sigma_bound = 15.;

	double superluminal = 1 - 9*pow(10, -15.); //This ensures that c_I^2 > 1 - order(10^(-15)).

	for (int i = 0 ; i<iterations; i++){
	  //Uses gsl to get a random number (uniformly distributed between -i_bound and i_bound)
	  //params.betappe[0] = gsl_rng_uniform(r)*2*pow(10, -a_bound) - pow(10, - a_bound); // ca (allowing it to be negative - violates positive energy constraint)
	  params.betappe[0] = gsl_rng_uniform(r)*pow(10, -a_bound); // ca (this parameter has to be positive because of energy constraints)
	  params.betappe[1] = 3.*params.betappe[0]*(1 + (2*gsl_rng_uniform(r) - 1.)*pow(10, -theta_bound)); // ctheta (Note special form here)
	  params.betappe[2] = (2*gsl_rng_uniform(r) - 1.)*pow(10, w_bound -1); // cw
	  params.betappe[3] = (2*gsl_rng_uniform(r) - 1.)*pow(10, -sigma_bound); // csigma
	  

	  //Random numbers with more variation of order
	  /*
	  params.betappe[0] = gsl_rng_uniform(r)*pow(10, - (gsl_rng_uniform(r)*16. + a_bound)); // ca (this parameter has to be positive because of energy constraints)
	  //params.betappe[1] = 3.*params.betappe[0]*(1 + gsl_rng_uniform(r)*2*pow(10, -(gsl_rng_uniform(r)*17 + theta_bound)) - pow(10, -(gsl_rng_uniform(r)*17 + theta_bound))); // ctheta
	  params.betappe[1] = (2*gsl_rng_uniform(r) - 1.)*pow(10, -(gsl_rng_uniform(r)*17 + theta_bound)); // ctheta
	  params.betappe[2] = (2*gsl_rng_uniform(r) - 1.)*pow(10, -(gsl_rng_uniform(r)*(20 + w_bound) - w_bound)); // cw
	  params.betappe[3] = (2*gsl_rng_uniform(r) - 1.)*pow(10, -(gsl_rng_uniform(r)*5 + sigma_bound)); // csigma
	  */

	  prep_source_parameters(&sp, &params,"EA_fully_restricted_v1_IMRPhenomD_NRT");//This should also run pre_calculate_EA_factors
	  //std::cout<<params.betappe[0]<<"\t"<<params.betappe[1]<<"\t"<<params.betappe[2]<<"\t"<<params.betappe[3]<<std::endl;

	  if(params.betappe[2] < (-params.betappe[3]/(1. - params.betappe[3])))
	    {
	      //Throws out points with cw < -csigma/(1 - csigma) because these violate the positive energy condition for the spin-1 mode (vector mode?)
	      //Note that the positive energy condition for the spin-2 modes is always satisfied and for the spin-0 mode is satisfied by requiring
	      //ca > 0 which is done when drawing random values of ca.
	      //Equation 40 of arXiv:gr-qc/0507059v3
	      num_bad_points++;
	      i--;
	      continue;
	    }
	  else if (isinf(sp.cV_EA) || isinf (sp.cS_EA))
	    {
	      //Pull GR limit points out of the data set so that they don't mess with my plotting
	      //Keep track of how many we pull
	      num_infspeeds++;
	      i--;
	      continue; 
	    }
	  else if(sp.cTsq_EA < superluminal || sp.cVsq_EA < superluminal || sp.cSsq_EA < superluminal)
	    {
	      //Throws out points with speeds not greater than or equal to the speed of light (these would produce Cherenkov-type radiation that is not observed)
	      //arXiv:hep-ph/0505211
	      num_bad_points++;
	      i--;
	      continue;
	    }
	  else if(sp.cT_EA - 1. < -3*pow(10, -15.) || sp.cT_EA -1. > 7*pow(10, -16.))
	    {
	      //Throws out points that don't obey the cT constraint from GW170817 and GRB170817A
	      //arXiv:1710.05834
	      num_bad_points++;
	      i--;
	      //std::cout<<"cT incorrect: "<<sp.cT_EA - 1<<std::endl;
	      //I wonder why all the incorrect points it prints have the same value? 
	      continue;
	    }
	  else if(abs(sp.alpha1_EA) > pow(10, -4.) || abs(sp.alpha2_EA) > pow(10, -7.))
	    {
	      //Throws out points that do not obey observational solar system constraints on alpha1 and alpha2
	      //arXiv:1403.7377 and arXiv:gr-qc/0509114
	      num_bad_points++;
	      i--;
	      continue; 
	    }

	  
	  output[i][0] = params.betappe[0];
	  output[i][1] = params.betappe[1];
	  output[i][2] = params.betappe[2];
	  output[i][3] = params.betappe[3];
	  output[i][4] = sp.cT_EA;
	  output[i][5] = sp.cS_EA;
	  output[i][6] = sp.cV_EA;
	  output[i][7] = sp.alpha_ppE_2T_0_EA;
	  output[i][8] = sp.gb1_EA;
	  output[i][9] = sp.abL_EA;
	  output[i][10] = sp.gX1_EA;
	  output[i][11] = sp.alpha1_EA;
	  output[i][12] = sp.alpha2_EA;
		
	  //std::cout<<sp.abL_EA<<std::endl; 
	  
	  cleanup_source_parameters(&sp,"EA_fully_restricted_v1_IMRPhenomD_NRT");
	}
	write_file("data/EA_parameter_MC.csv",output, iterations, dim);
	std::cout<<"Wrote data to 'data/EA_parameter_MC.csv'"<<std::endl;
	std::cout<<"Threw out "<<num_bad_points<<" unphysical data points."<<std::endl;
	std::cout<<"Found cS or cV to be infinite "<<num_infspeeds<<" times. Removed these cases from data."<<std::endl;
	gsl_rng_free(r); 
	deallocate_2D_array(output, iterations, dim);
	delete [] params.betappe;
	delete [] params.bppe;
	return 0;
}
int EA_fully_restricted_test(int argc, char *argv[])
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
	params.NSflag1 = true;
	params.NSflag2 = true;
	params.horizon_coord = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 1.4;
	params.mass2 = 1.3;
	params.tc = 6;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;
	params.gmst=3;

	params.Nmod = 4;
	params.bppe = new double[4];
	params.bppe[0] = -13;
	params.bppe[1] = -13;
	params.bppe[2] = -13;
	params.bppe[3] = -13;
	params.betappe = new double[4];
	//params.betappe[0] = .001;
	params.betappe[0] = 1;
	params.betappe[1] = 1;
	params.betappe[2] = 1;
	params.betappe[3] = 1;
	params.tidal1 = 10;
	params.tidal2 = 10;
	
	int length = 10000;
	double freqs[length];
	//double fhigh =pow(6,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	double fhigh =4000;//pow(6,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	double flow = pow(100,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	std::cout<<flow<< " " <<fhigh<<std::endl;
	double delta_f = (fhigh - flow)/length;
	for(int i = 0 ; i<length; i++){
		freqs[i] = flow + delta_f*i;
	}
	std::complex<double> responseED[length];
	double **output = new double*[2];
	output[0] = new double[length];
	output[1] = new double[length];

	//std::complex<double> *hplus = new std::complex<double>[length];
	//std::complex<double> *hcross = new std::complex<double>[length];
	
	waveform_polarizations<double> wp;
	assign_polarizations("EA_fully_restricted_v1_IMRPhenomD_NRT",&wp);
	wp.allocate_memory(length);
	//wp.hplus = hplus;	
	//wp.hcross = hcross;	


	//fourier_detector_response(freqs, length, responseED, "Hanford", "EA_fully_restricted_v1_IMRPhenomD_NRT",&params, (double*)NULL);
	fourier_waveform(freqs, length, &wp, "EA_fully_restricted_v1_IMRPhenomD_NRT",&params);
	//fourier_detector_response(freqs, length, responseED, "Hanford", "EA_fully_restricted_v1_IMRPhenomD",&params, (double*)NULL);
	//std::string dets[2] = {"Hanford","Livingston"};
	//std::complex<double> **responses= new std::complex<double>*[2];
	//responses[0] = new std::complex<double>[length];
	//responses[1] = new std::complex<double>[length];
	//create_coherent_GW_detection_reuse_WF(dets, 2, freqs, length, &params, "EA_fully_restricted_v1_IMRPhenomD",responses);
	//delete[]responses[0];
	//delete[]responses[1];
	//delete [] responses;

	//for(int i =0 ; i<length ; i++){
	//	output[0][i] = std::real(responseED[i]);
	//	output[1][i] = std::imag(responseED[i]);
	//}
	//write_file("data/EA_fully_restricted_v1_waveform.csv",output, 2,length);
	//################################
	//delete [] params.bppe; delete[] params.betappe;delete [] output[0];delete [] output[1];delete [] output;
	delete[] params.betappe;delete [] params.bppe;delete [] output[0];delete [] output[1];delete [] output;
	wp.deallocate_memory();
	//delete [] hplus; delete [] hcross;
	return 0;
}
int time_domain_testing(int argc, char *argv[])
{

	//#############################################################################
	gen_params p ;
	p.mass1 = 20;
	p.mass2 = 10;
	p.Luminosity_Distance = 100;
	p.x0 = pow((p.mass1+p.mass2)*MSOL_SEC*10,2./3.);
	p.incl_angle = 0;
	p.RA = .1;
	p.DEC = .1;
	p.psi = .1;
	p.phi = .1;
	p.theta = .1;
	p.horizon_coord=false;
	p.gmst= .21;
	
	int length = 1000;
	std::complex<double> *hplus = new std::complex<double>[length];
	std::complex<double> *hcross = new std::complex<double>[length];
	std::complex<double> *response = new std::complex<double>[length];
	
	waveform_polarizations<double> wp;
	wp.hplus = hplus;	
	wp.hcross = hcross;	
	double times[length]; 
	for(int i=0 ; i<length; i++){
		times[i] = -100 + 90*i/length;
	}
	time_waveform(times, length, &wp, "TaylorT2", &p);	
	time_detector_response(times, length, response, "Hanford", "TaylorT2", &p);	
	
	delete[] hplus;
	delete[] hcross;
	delete[] response;

	return 0;
}
int BHEvaporation_test(int argc, char *argv[])
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

	params.Nmod = 1;
	params.bppe = new double[1];
	params.bppe[0] = -13;
	params.betappe = new double[1];
	//params.betappe[0] = .001;
	params.betappe[0] = 100;
	
	int length = 10000;
	double freqs[length];
	double fhigh =10* pow(6,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	double flow = pow(100,-3./2)/(M_PI * (params.mass1+params.mass2)*MSOL_SEC);
	std::cout<<flow<< " " <<fhigh<<std::endl;
	double delta_f = (fhigh - flow)/length;
	for(int i = 0 ; i<length; i++){
		freqs[i] = flow + delta_f*i;
	}
	std::complex<double> responseED[length];
	std::complex<double> responseBHE[length];
	double **output = new double*[2];
	output[0] = new double[length];
	output[1] = new double[length];


	fourier_detector_response(freqs, length, responseED, "Hanford", "ExtraDimension_IMRPhenomPv2",&params, (double*)NULL);
	for(int i =0 ; i<length ; i++){
		output[0][i] = std::real(responseED[i]);
		output[1][i] = std::imag(responseED[i]);
	}
	write_file("data/ED_waveform.csv",output, 2,length);
	//################################
	double Z = Z_from_DL(params.Luminosity_Distance,"PLANCK15");
	double ten_micrometer = 10.e-6/c;
	params.betappe[0] = -2.8e-7 * ( pow_int((1+Z) / params.mass1,2)+pow_int((1+Z) / params.mass2,2)) * params.betappe[0]/pow_int(ten_micrometer,2) * MSOL_SEC/T_year;

	fourier_detector_response(freqs, length, responseBHE, "Hanford", "BHEvaporation_IMRPhenomPv2",&params, (double*)NULL);
	for(int i =0 ; i<length ; i++){
		output[0][i] = std::real(responseBHE[i]);
		output[1][i] = std::imag(responseBHE[i]);
	}
	std::complex<double> ave=0;
	for (int i = 0 ; i<length ; i++){
		if (std::abs(responseED[i] ) > 0 || std::abs(responseBHE[i] ) > 0){
			ave+= std::abs(responseED[i] -  responseBHE[i])*2./(std::abs(responseED[i]) +  std::abs(responseBHE[i]));
		}
	}
	ave/=length;
	std::cout<<"Average fractional difference: "<<ave<<std::endl;
	write_file("data/BHE_waveform.csv",output, 2,length);
	
	
	

	delete [] params.bppe; delete[] params.betappe;delete [] output[0];delete [] output[1];delete [] output;
	return 0;
}
int polarization_testing(int argc, char *argv[])
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
	params.theta = 2.;
	params.phi = -1.1;
	params.f_ref = 20;
	params.horizon_coord =false;
	params.NSflag1 = false;
	params.NSflag2 = false;
	params.shift_time=true;
	params.shift_phase=true;
	
	params.mass1 = 3.6;
	params.mass2 = 2.8;
	params.tc = 0;
	params.equatorial_orientation = false;
	params.psi = 1.;
	params.incl_angle = M_PI/3.;
	params.gmst=3;

	
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
	
	fourier_detector_response(freqs, length, response, "Hanford", "polarization_test_IMRPhenomD",&params, (double*)NULL);
	for(int i =0 ; i<length ; i++){
		output[0][i] = std::real(response[i]);
		output[1][i] = std::imag(response[i]);
	}
	write_file("data/polarization_test.csv",output, 2,length);
	

	delete [] output[0];delete [] output[1];delete [] output;
	return 0;
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
	params.bppe = new double[1];
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
	waveform_polarizations<double> wp;
	wp.hplus = hpg;	
	wp.hcross = hcg;	
	fourier_waveform(frequency, length, &wp, "gIMRPhenomPv2",&params);
	wp.hplus = hpppE;	
	wp.hcross = hcppE;	
	fourier_waveform(frequency, length, &wp, "ppE_IMRPhenomPv2_Inspiral",&params);
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
	bool P = false;
	bool NRT = true;
	gsl_rng_env_setup();	
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	gsl_rng_set(r,10);
	int iterations = 100;
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
		double alpha[17];
		for (int j = 0 ; j<17; j++){
		  //alpha[j] = 0.1; 
		  alpha[j] = gsl_rng_uniform(r);
		}
		const REAL8 s1x = -.1+alpha[0]*.2, s1y=-.2+alpha[1]*.3,s1z=-.4+alpha[2]*.6;
		const REAL8 s2x = -.1+alpha[3]*.2, s2y=-.2+alpha[4]*.3,s2z=-.4+alpha[5]*.6; 
		//const REAL8 s1x = 0, s1y=0,s1z=-.4+alpha[2]*.6;
		//const REAL8 s2x = 0, s2y=0,s2z=-.4+alpha[5]*.6;
		//const REAL8 s1x = 0.0, s1y=0.0,s1z=0.0;
		//const REAL8 s2x =0.0, s2y=0.0,s2z=0.0;
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
		//double tempm1=1+5*alpha[10],tempm2 = 1+5*alpha[11];
		double tempm1,tempm2 ;
		if(NRT){
			tempm1 = 1+1*alpha[10];
			tempm2 = 1+1*alpha[10];
		}
		else{
			tempm1 = 1+15*alpha[10];
			tempm2 = 1+15*alpha[11];

		}
			
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
		//const REAL8 f_min = .0017*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		const REAL8 f_min = .002*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		const REAL8 f_max = .15*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		//const REAL8 f_max = .15*LAL_MSUN_SI/MSOL_SEC/(m1_SI+m2_SI);
		//const REAL8 f_max = 1600;
		int length = 4016;
		//int length = 131072;
		double deltaf = (f_max-f_min)/(length-1);
		IMRPhenomP_version_type  version = IMRPhenomPv2_V;
		LALDict *extraParams = XLALCreateDict();
		//alpha[15] = 0;
		//alpha[16] = 0; 
		//alpha[15] = 2;
		//alpha[16] = .8;
		REAL8 lambda1 = 100*fabs(alpha[15]) + 1.; //this prevents tidal deformability from being less than 1
		REAL8 lambda2 = 100*fabs(alpha[16]) + 1.;
		/*if (alpha[15]<alpha[16]){
		  lambda1= 100*fabs(alpha[15]) ;
		      lambda2	= 100*fabs(alpha[16]) ;

		}
		else{
		  lambda1= 100*fabs(alpha[16]) ;
		      lambda2	= 100*fabs(alpha[15]) ;
		      }*/
		//std::cout<<"TIDAL LOVE NUMBERS: "<<lambda1<<" "<<lambda2<<std::endl;
		NRTidal_version_type NRT_v=NRTidalv2_V;

		double q0 = 0.1940;
		double q1 = 0.09163;
		double q2 = 0.04812;
		double q3 = -0.004283;
		double q4 = 0.00012450;

		double quad1 = exp(q0 + q1*log(lambda1) + q2*pow(log(lambda1), 2.) + q3*pow(log(lambda1), 3.) + q4*pow(log(lambda1), 4.));
		double quad2 = exp(q0 + q1*log(lambda2) + q2*pow(log(lambda2), 2.) + q3*pow(log(lambda2), 3.) + q4*pow(log(lambda2), 4.));
		//double quad1 = 1;
		//double quad2 = 1;
	        //std::cout<<"quad1: "<<quad1<<"\t quad2: "<<quad2<<std::endl; 

		

		//NRTidal_version_type tidalType= NoNRT_V;

		REAL8Sequence *freqs = XLALCreateREAL8Sequence(length);
		double *frequencies = new double[length];
		for(int i = 0 ; i<length; i++){
			freqs->data[i] = f_min + i * deltaf;
			frequencies[i] = f_min+i*deltaf;
		}
		//const REAL8 f_ref = frequencies[0];
		const REAL8 f_ref = 20.;
		clock_t start,end;
		start = clock();
		XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJ, &alpha0, &phi_aligned, &zeta_polariz, m1_SI, m2_SI, f_ref, phiRef, incl, s1x,s1y,s1z,s2x,s2y,s2z,version);
		//XLALSimIMRPhenomP(&hptilde,&hctilde,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,deltaf,f_min,f_max,f_ref, version, extraParams);
		if(P){
			//XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, version, tidalType,extraParams);
			if(!NRT){
				NRT_v = NoNRT_V;
				XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, version,NRT_v,extraParams);
			}
			else{
				//XLALSimIMRPhenomPFrequencySequence(&hptilde,&hctilde,freqs,chi1_l,chi2_l,chip,thetaJ,m1_SI,m2_SI,distance,alpha0,phi_aligned,f_ref, NRT_v,extraParams);
			}
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
			if(!NRT){
				NRT_v = NoNRT_V;
				XLALSimIMRPhenomDFrequencySequence(&hptilde,freqs,phiRef,f_ref,m1_SI,m2_SI,s1z,s2z,distance, extraParams,NRT_v);
			}
			else{
				XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams, quad1-1);
				XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams, quad2-1);
				XLALSimIMRPhenomDNRTidalFrequencySequence(&hptilde,freqs,phiRef,f_ref,distance,m1_SI,m2_SI,s1z,s2z, lambda1,lambda2,extraParams,NRT_v);
			}
			hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, freqs->data[1]-freqs->data[0], &lalStrainUnit, length);
			for(int i = 0 ; i<length; i++){
				gsl_complex f2 = gsl_complex_rect(0.,-cos(incl));
				(hctilde->data->data)[i] = gsl_complex_mul(f2,(hptilde->data->data)[i]);
				gsl_complex f1 = gsl_complex_rect(0.5*(1. + cos(incl)*cos(incl)),0);
				(hptilde->data->data)[i] = gsl_complex_mul(f1,(hptilde->data->data)[i]);

			}
		}
		COMPLEX16FrequencySeries *det = XLALCreateCOMPLEX16FrequencySeries("det: FD waveform", &ligotimegps_zero, 0.0, freqs->data[1]-freqs->data[0], &lalStrainUnit, length);
		double fplus,fcross,fb,fl,fx,fy;
		double fplusG,fcrossG,fbG,flG,fxG,fyG;
		XLALComputeDetAMResponse(&fplus,&fcross, LALD.response, RA,DEC,psi,gmst);
		XLALComputeDetAMResponseExtraModes(&fplus,&fcross,&fb,&fl,&fx,&fy, LALD.response, RA,DEC,psi,gmst);
		end = clock();
		
		det_res_pat<double> r_pat;
		r_pat.Fplus = &fplusG;
		r_pat.Fcross = &fcrossG;
		r_pat.Fx = &fxG;
		r_pat.Fy = &fyG;
		r_pat.Fb = &fbG;
		r_pat.Fl = &flG;
		bool pol_arr[6] = {true, true, true, true, true, true};
		r_pat.active_polarizations = &pol_arr[0];
		
		detector_response_functions_equatorial(DETECTOR,RA,DEC,psi,gmst, &r_pat);
		std::cout<<"Fractional error on F+/Fx/Fx/Fy/Fb/Fl: "<<(fplus-fplusG)/fplus<<" "<<(fcross-fcrossG)/fcross<<(fx-fxG)/fx<<" "<<(fy-fyG)/fy<<" "<<(fb-fbG)/fb<<" "<<(fl-flG)/fl<<" "<<std::endl;
		for(int i = 0 ; i<length ; i++){
			(det->data->data)[i]=gsl_complex_add(
			gsl_complex_mul(gsl_complex_rect(fplus,0.),(hptilde->data->data)[i]),
			gsl_complex_mul(gsl_complex_rect(fcross,0.0) , (hctilde->data->data)[i]));
		}
		times[k][0] = (double)(end-start)/(CLOCKS_PER_SEC);
		//std::cout<<"LAL timing: "<<(double)(end-start)/(CLOCKS_PER_SEC)<<std::endl;
		//###############################################################################
		gen_params param;
		param.mass1 = m1_SI/LAL_MSUN_SI;	
		param.mass2 = m2_SI/LAL_MSUN_SI;
		//std::cout<<"mass1: "<<param.mass1<<std::endl;
		//std::cout<<"mass2: "<<param.mass2<<std::endl; 
		param.Luminosity_Distance = distance/MPC_M;
		//std::cout<<"Luminosity distance: "<<param.Luminosity_Distance<<std::endl; 
		param.incl_angle = incl;
		param.phiRef = phiRef;
		param.f_ref = f_ref;
		param.spin1[0] = s1x;
		param.spin1[1] = s1y;
		param.spin1[2] = s1z;
		param.spin2[0] = s2x;
		param.spin2[1] = s2y;
		param.spin2[2] = s2z;
		//std::cout<<"spin1[0]: "<<param.spin1[0]<<"\t spin1[1]:"<<param.spin1[1]<<"\t spin2[2]:"<<param.spin1[2]<<std::endl; 
		//std::cout<<"spin2[0]: "<<param.spin2[0]<<"\t spin2[1]:"<<param.spin2[1]<<"\t spin2[2]:"<<param.spin2[2]<<std::endl; 

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
		param.tidal1 =lambda1 ;
		param.tidal2 =lambda2 ;
		//std::cout<<"tidal1: "<<param.tidal1<<"\t tidal2: "<<param.tidal2<<std::endl;
		
		std::complex<double> *response = new std::complex<double>[length];
		std::string method ;
		if(P){
			if(NRT){
				method = "IMRPhenomPv2_NRT";
			}
			else{
				method = "IMRPhenomPv2";
			}
		}
		else {
			if(NRT){
				method = "IMRPhenomD_NRT";
			}
			else{
				method = "IMRPhenomD";
			}
		}
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
		//std::cout<<output[0][1]<<"\t"<< output[0][2]<<"\t"<<output[0][3]<<"\t"<<output[0][4]<<std::endl; 
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
	std::cout<<"3 --- test PNSeries waveforms"<<std::endl;
	std::cout<<"4 --- test polarizations waveforms"<<std::endl;
	std::cout<<"5 --- test BHEvaporation waveforms"<<std::endl;
	std::cout<<"6 --- test time domain waveforms"<<std::endl;
	std::cout<<"7 --- EA fully restricted waveform"<<std::endl;
	std::cout<<"8 --- EA fully restricted waveform parameterization"<<std::endl;
	std::cout<<"9 --- EA fully restricted consistency test"<<std::endl;
}
