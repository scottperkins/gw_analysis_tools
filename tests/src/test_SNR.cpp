#include <gwat/waveform_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/io_util.h>
#include <iostream>


void RT_ERROR_MSG();
int compare_GSL_simps_gl_MBH(int argc, char *argv[]);
int compare_GSL_simps_gl_stellar(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	std::cout<<"TESTING SNR CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"COMPARISON MBH"<<std::endl;
		return compare_GSL_simps_gl_MBH(argc,argv);
	}
	else if(runtime_opt == 1){
		std::cout<<"COMPARISON Stellar"<<std::endl;
		return compare_GSL_simps_gl_stellar(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int compare_GSL_simps_gl(int argc, char *argv[])
{
		
	return 0;

}

int compare_GSL_simps_gl_MBH(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;
	params.mass1 = 10e5;
	params.mass2 = 5e5;
	params.spin1[2] = .2;
	params.spin2[2] = .1;
	params.spin1[1] = .2;
	params.spin2[1] = .1;
	params.spin1[0] = .2;
	params.spin2[0] = .1;
	params.phiRef = 2.;
	params.tc = T_year;
	params.f_ref = 1e-5;
	params.NSflag1=false;
	params.NSflag2=false;
	params.shift_phase=false;
	params.shift_time=false;
	params.equatorial_orientation=true;
	params.horizon_coord=false;
	params.Luminosity_Distance  = 400;
	params.RA =1.;
	params.DEC =-1.;
	params.phi_l =2.;
	params.theta_l =2.;
	params.gmst=1.;
	params.LISA_phi0=0;
	params.LISA_alpha0=0;
	
	
	int length = T_year;
	double deltaf = 1./(T_year);
	double *freqs = new double[length];
	
	int lengthgl1 = 500;
	int lengthgl2 = 500;
	double *freqsgl1 = new double[lengthgl1];
	double *freqsgl2 = new double[lengthgl2];
	double *w1 = new double[lengthgl1];
	double *w2 = new double[lengthgl2];

	gauleg(1e-5, 1, freqsgl1, w1, lengthgl1);
	gauleg(log10(1e-5), log10(1), freqsgl2, w2, lengthgl2);
	for(int i =0; i<lengthgl2; i++){
		freqsgl2[i] = pow(10.,freqsgl2[i]);
	}

	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1);


	double snr_gsl, snr_simps, snr_gl1,snr_gl2;
	std::string noise_curve="LISA_CONF"	;
	std::string method= "IMRPhenomPv2";
	std::string detector="LISA";
	int iterations = 50;
	double **output = new double*[iterations];
	std::cout<<"Starting loop -- finished prep"<<std::endl;
	for(int i = 0 ; i<iterations; i++){
		output[i]=new double[10];
		double m1 = 10e4+1e5*gsl_rng_uniform(r);
		double m2 = 10e4+1e4*gsl_rng_uniform(r);
		if(m1>m2){
			params.mass1 = m1;
			params.mass2 = m2;
		}
		else{
			params.mass1 = m1;
			params.mass2 = m2;
		}
		params.spin1[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin2[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin1[1] = -.4 + gsl_rng_uniform(r)*.8;
		params.spin2[1] = .1;
		params.spin1[0] = .2;
		params.spin2[0] = 0;//-.4 + gsl_rng_uniform(r)*1.;
		params.Luminosity_Distance  = 50 + gsl_rng_uniform(r)*1000;
		params.RA =gsl_rng_uniform(r)*2*M_PI;
		params.DEC =-1.5 + gsl_rng_uniform(r) * M_PI;
		params.phi_l =gsl_rng_uniform(r)*2*M_PI;
		params.theta_l =gsl_rng_uniform(r)*M_PI;
		params.gmst=gsl_rng_uniform(r)*2*M_PI;
		transform_orientation_coords(&params,method, detector);

		double bounds[2];
		bool autodiff = true;
		clock_t start = clock();
		observation_bounds(1, T_year, detector,noise_curve, method, &params, bounds,autodiff);
		std::cout<<"Observation bounds ad: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;	

		for(int i = 0 ; i<length; i++){
			freqs[i] = bounds[0] + i*deltaf;
		}


		start = clock();
		calculate_snr_gsl(&snr_gsl, noise_curve,detector,method, &params, bounds[0], bounds[1], 1.e-12);
		double gsl_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		snr_simps=0;
		snr_simps = calculate_snr(noise_curve, detector,method, &params,freqs,length,"SIMPSONS",(double *)NULL,false);
		double simps_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		gauleg(bounds[0], bounds[1], freqsgl1, w1, lengthgl1);
		snr_gl1 = calculate_snr(noise_curve, detector,method, &params,freqsgl1,lengthgl1,"GAUSSLEG",w1,false);
		double gl1_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		gauleg(log10(bounds[0]), log10(bounds[1]), freqsgl2, w2, lengthgl2);
		for(int i =0; i<lengthgl2; i++){
			freqsgl2[i] = pow(10.,freqsgl2[i]);
		}
		snr_gl2 = calculate_snr(noise_curve, detector,method, &params,freqsgl2,lengthgl2,"GAUSSLEG",w2,true);
		double gl2_t = (double)(clock()-start)/CLOCKS_PER_SEC;

		output[i][0]=snr_gsl;
		output[i][1]=snr_simps;
		output[i][2]=snr_gl1;
		output[i][3]=snr_gl2;
		output[i][4]=gsl_t;
		output[i][5]=simps_t;
		output[i][6]=gl1_t;
		output[i][7]=gl2_t;
		
		autodiff = false;
		start = clock();
		observation_bounds(1, T_year, detector,noise_curve, method, &params, bounds,autodiff);
		std::cout<<"Observation bounds N: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<std::endl;	
		gauleg(bounds[0], bounds[1], freqsgl1, w1, lengthgl1);
		snr_gl1 = calculate_snr(noise_curve, detector,method, &params,freqsgl1,lengthgl1,"GAUSSLEG",w1,false);
		gauleg(log10(bounds[0]), log10(bounds[1]), freqsgl2, w2, lengthgl2);
		for(int i =0; i<lengthgl2; i++){
			freqsgl2[i] = pow(10.,freqsgl2[i]);
		}
		snr_gl2 = calculate_snr(noise_curve, detector,method, &params,freqsgl2,lengthgl2,"GAUSSLEG",w2,true);

		output[i][8] = snr_gl1;
		output[i][9] = snr_gl2;
		printProgress((double)i/iterations);
	}
	std::cout<<std::endl;
	
	write_file("data/snr_comp_MBH.csv",output,iterations, 10);
	delete [] freqsgl1;
	delete [] freqsgl2;
	delete [] w1;
	delete [] w2;
	for(int i = 0 ; i<iterations; i++){
		delete [] output[i];	
	}
	delete [] output;
	return 0;
}
int compare_GSL_simps_gl_stellar(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;
	params.mass1 = 10;
	params.mass2 = 5;
	params.spin1[2] = .2;
	params.spin2[2] = .1;
	params.spin1[1] = .2;
	params.spin2[1] = .1;
	params.spin1[0] = .2;
	params.spin2[0] = .1;
	params.phiRef = 2.;
	params.tc = 16.;
	params.f_ref = 20;
	params.NSflag1=false;
	params.NSflag2=false;
	params.shift_phase=false;
	params.shift_time=false;
	params.equatorial_orientation=false;
	params.horizon_coord=false;
	params.Luminosity_Distance  = 400;
	params.RA =1.;
	params.DEC =-1.;
	//params.psi =2.;
	//params.incl_angle =2.;
	params.theta_l =2.;
	params.phi_l =2.;
	params.gmst=1.;
	
	double Tobs = 4;
	int length =(2048 - 10)*Tobs ;
	double deltaf = 1./(Tobs);
	double *freqs = new double[length];
	for(int i = 0 ; i<length; i++){
		freqs[i] = 10 + i*deltaf;
	}
	
	int lengthgl1 = 500;
	int lengthgl2 = 500;
	double *freqsgl1 = new double[lengthgl1];
	double *freqsgl2 = new double[lengthgl2];
	double *w1 = new double[lengthgl1];
	double *w2 = new double[lengthgl2];

	gauleg(10, 2048, freqsgl1, w1, lengthgl1);
	gauleg(log10(10), log10(2048), freqsgl2, w2, lengthgl2);
	for(int i =0; i<lengthgl2; i++){
		freqsgl2[i] = pow(10.,freqsgl2[i]);
	}
	

	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1);


	double snr_gsl, snr_simps, snr_gl1,snr_gl2;
	std::string noise_curve="aLIGO_analytic"	;
	std::string method= "IMRPhenomPv2";
	std::string detector="Hanford";
	int iterations = 100;
	double **output = new double*[iterations];
	std::cout<<"Starting loop -- finished prep"<<std::endl;
	for(int i = 0 ; i<iterations; i++){
		output[i]=new double[8];
		double m1 = 10+1e2*gsl_rng_uniform(r);
		double m2 = 10+1e2*gsl_rng_uniform(r);
		if(m1>m2){
			params.mass1 = m1;
			params.mass2 = m2;
		}
		else{
			params.mass1 = m1;
			params.mass2 = m2;
		}
		params.spin1[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin2[2] = -.3+gsl_rng_uniform(r)*.5;
		params.spin1[1] = -.4 + gsl_rng_uniform(r)*.8;
		params.spin2[1] = .1;
		params.spin1[0] = .2;
		params.spin2[0] = 0;//-.4 + gsl_rng_uniform(r)*1.;
		params.Luminosity_Distance  = 50 + gsl_rng_uniform(r)*1000;
		params.RA =gsl_rng_uniform(r)*2*M_PI;
		params.DEC =-1.5 + gsl_rng_uniform(r) * M_PI;
		params.phi_l =gsl_rng_uniform(r)*2*M_PI;
		params.theta_l =gsl_rng_uniform(r)*M_PI;
		params.gmst=gsl_rng_uniform(r)*2*M_PI;
		transform_orientation_coords(&params,method, "");

		clock_t start = clock();
		calculate_snr_gsl(&snr_gsl, noise_curve,detector,method, &params, 10, 2048, 1.e-12);
		double gsl_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		snr_simps = calculate_snr(noise_curve, detector,method, &params,freqs,length,"SIMPSONS",(double *)NULL,false);
		double simps_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		snr_gl1 = calculate_snr(noise_curve, detector,method, &params,freqsgl1,lengthgl1,"GAUSSLEG",w1,false);
		double gl1_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		start = clock();
		snr_gl2 = calculate_snr(noise_curve, detector,method, &params,freqsgl2,lengthgl2,"GAUSSLEG",w2,true);
		double gl2_t = (double)(clock()-start)/CLOCKS_PER_SEC;
		output[i][0]=snr_gsl;
		output[i][1]=snr_simps;
		output[i][2]=snr_gl1;
		output[i][3]=snr_gl2;
		output[i][4]=gsl_t;
		output[i][5]=simps_t;
		output[i][6]=gl1_t;
		output[i][7]=gl2_t;
		printProgress((double)i/iterations);
	}
	std::cout<<std::endl;
	
	write_file("data/snr_comp_stellar.csv",output,iterations, 8);
	delete [] freqsgl1;
	delete [] freqsgl2;
	delete [] w1;
	delete [] w2;
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
	std::cout<<"0 --- Compare SNR methods MBH"<<std::endl;
	std::cout<<"1 --- Compare SNR methods stellar mass"<<std::endl;
}
