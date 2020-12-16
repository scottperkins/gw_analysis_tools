#include <iostream>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>
#include <gwat/util.h>
#include <gwat/waveform_util.h>
#include <gsl/gsl_rng.h>


int test_response_functions(int argc, char *argv[]);
void RT_ERROR_MSG();

int main(int argc , char* argv[]){
	std::cout<<"TESTING RESPONSE CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"Response function testing"<<std::endl;
		return test_response_functions(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}

int test_response_functions(int argc, char *argv[]){
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	std::string detectors[3] = {"Hanford","Livingston", "Virgo"};
	int iterations = 50;
	double **output= new double*[iterations];
	bool active_polarizations[6] = {true, true, false,false,false,false};
	int dim = 19;//RA,DEC,PSI,IOTA,gmst,THETA_S,PHI_S,THETA_L,PHI_L,THETA_S_ecl, PHI_S_ecl, THETA_L_ecl, PHI_L_ECL,(F+,Fx)x3
	for (int i = 0 ; i< iterations; i++){
		output[i] = new double[dim];
		double ra = 0.01 + gsl_rng_uniform(r)*1.9*M_PI;
		double dec = M_PI/2 - .1 + gsl_rng_uniform(r)*.9*M_PI;
		double phi_l =  .1 + gsl_rng_uniform(r)*1.9*M_PI;
		double theta_l =  .1 + gsl_rng_uniform(r)*.9*M_PI;
		//double gmst = i*1.9*M_PI/iterations+.01;	
		double gmst = 0;	
		gen_params p;
		p.RA = ra;
		p.DEC = dec;
		p.theta_l = theta_l;
		p.phi_l	 = phi_l;
		p.gmst = gmst;
		p.phiRef = 0;
		transform_orientation_coords(&p,"IMRPhenomD","");
		double incl = p.incl_angle;
		double psi = p.psi;
		
		double theta_s = M_PI/2. - dec;
		double phi_s = ra;
		double theta_s_ecl, phi_s_ecl, theta_l_ecl, phi_l_ecl;
		ecl_from_eq(theta_s, phi_s, &theta_s_ecl, &phi_s_ecl);
		ecl_from_eq(theta_l, phi_l, &theta_l_ecl, &phi_l_ecl);
		
		output[i][0] = ra;
		output[i][1] = dec;
		output[i][2] = psi;
		output[i][3] = incl;
		output[i][4] = gmst;
		output[i][5] = theta_s;
		output[i][6] = phi_s;
		output[i][7] = theta_l;
		output[i][8] = phi_l;
		output[i][9] = theta_s_ecl;
		output[i][10] = phi_s_ecl;
		output[i][11] = theta_l_ecl;
		output[i][12] = phi_l_ecl;
		double Fplus;
		double Fcross;
		det_res_pat<double> r_pat ;
		r_pat.Fplus = &Fplus;
		r_pat.Fcross = &Fcross;
		r_pat.active_polarizations= &active_polarizations[0];
		for(int j = 0 ; j<3; j++){
			detector_response_functions_equatorial(detectors[j], ra, dec, psi, gmst, &r_pat);
			output[i][j*2+13] = Fplus;
			output[i][j*2+14] = Fcross;
		}

	}
	write_file("data/response_functions.csv",output,iterations,dim);	
	for(int i = 0 ; i<iterations; i++){
		delete [] output[i];	
	}	
	delete [] output;	
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Response function testing"<<std::endl;
}
