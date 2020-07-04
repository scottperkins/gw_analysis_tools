#include <gwat/waveform_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_sampler.h>
#include <gwat/mcmc_gw.h>
#include <gwat/detector_util.h>
#include <gwat/IMRPhenomP.h>
#include <gwat/fisher.h>
#include <iostream>


void RT_ERROR_MSG();
int mcmc_standard_test(int argc, char *argv[]);
int mcmc_injection(int argc, char *argv[]);
int mcmc_real_data(int argc, char *argv[]);
int mcmc_rosenbock(int argc, char *argv[]);
int mcmc_output_class(int argc, char *argv[]);
double log_test (double *c,mcmc_data_interface *interface,void *parameters);
double log_test_prior (double *c,mcmc_data_interface *interface,void *parameters);
double log_rosenbock (double *c,mcmc_data_interface *interface,void *parameters);
double log_rosenbock_prior (double *c,mcmc_data_interface *interface,void *parameters);
void fisher_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters);
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters);
double standard_log_prior_P(double *pos, mcmc_data_interface *interface,void *parameters);
double chirpmass_eta_jac(double chirpmass, double eta);
double T_mcmc_gw_tool ;
double standard_log_prior_dCS(double *pos, mcmc_data_interface *interface,void *parameters);
//void fisher_rosenbock(double *c,double **fisher,  mcmc_data_interface *interface,void *parameters);

int main(int argc, char *argv[])
{
	std::cout<<"TESTING MCMC CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"MCMC student t"<<std::endl;
		return mcmc_standard_test(argc,argv);
	}
	else if(runtime_opt == 1){
		std::cout<<"MCMC Injected GW"<<std::endl;
		return mcmc_injection(argc,argv);
	}
	else if(runtime_opt == 2){
		std::cout<<"MCMC data GW"<<std::endl;
		return mcmc_real_data(argc,argv);
	}
	else if(runtime_opt == 3){
		std::cout<<"MCMC N-dimensional Rosenbock"<<std::endl;
		return mcmc_rosenbock(argc,argv);
	}
	else if(runtime_opt == 4){
		std::cout<<"MCMC output class testing"<<std::endl;
		return mcmc_output_class(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}

int mcmc_output_class(int argc, char *argv[])
{
	int chain_N = 3;
	int dim = 3;
	mcmc_sampler_output output(chain_N,dim);
	double chain_temps[chain_N] = {1.,2.,1.};
	output.populate_chain_temperatures(chain_temps);
	for(int i = 0 ; i<chain_N; i++){
		std::cout<<chain_temps[i]<<std::endl;
	}
	std::cout<<"Cold chains: "<<output.cold_chain_number<<std::endl;
	std::cout<<"Cold chains id 1: "<<output.cold_chain_ids[0]<<std::endl;
	std::cout<<"Cold chains id 2: "<<output.cold_chain_ids[1]<<std::endl;
	double ***init_output = new double**[chain_N];
	double ***init_LL_LP = new double**[chain_N];
	int steps = 5;
	int init_positions[chain_N] = {5,3,5};
	for(int i = 0 ; i<chain_N; i++){
		init_output[i] = new double*[steps];
		init_LL_LP[i] = new double*[steps];
		for(int j =0 ; j<steps; j++){
			init_output[i][j] = new double[dim];
			init_LL_LP[i][j] = new double[2];
			if(j < init_positions[i]){
				init_LL_LP[i][j][0] = 10;
				init_LL_LP[i][j][1] = 1;
				for(int k =0 ; k<dim; k++){
					init_output[i][j][k] = i+j+k;
				}
			}
			else{
				init_LL_LP[i][j][0] = 0;
				init_LL_LP[i][j][1] = 0;
				for(int k =0 ; k<dim; k++){
					init_output[i][j][k] = 0;
				}
			}
		}
	}
	output.populate_initial_output(init_output, init_LL_LP,init_positions);
	std::cout<<"Initial output"<<std::endl;
	for(int i = 0 ; i<chain_N; i++){
		for(int j =0 ; j<init_positions[i]; j++){
			for(int k =0 ; k<dim; k++){
				std::cout<<output.output[i][j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
	output.trim_lengths[0]=2;
	output.trim_lengths[1]=2;
	output.trim_lengths[2]=2;


	std::cout<<"Creating data dump"<<std::endl;
	output.create_data_dump(false,true, "./data/testT.hdf5");
	output.create_data_dump(false,false, "./data/test.hdf5");
	//exit(1);

	double ***app_output = new double**[chain_N];
	double ***app_LL_LP = new double**[chain_N];
	int app_steps = 8;
	int app_positions[chain_N] = {8,5,8};
	for(int i = 0 ; i<chain_N; i++){
		app_output[i] = new double*[app_steps];
		app_LL_LP[i] = new double*[app_steps];
		for(int j =0 ; j<app_steps; j++){
			app_output[i][j] = new double[dim];
			app_LL_LP[i][j] = new double[2];
			if(j < app_positions[i]){
				app_LL_LP[i][j][0] = 10*(i+j);
				app_LL_LP[i][j][1] = 11*(i+j);
				for(int k =0 ; k<dim; k++){
					app_output[i][j][k] = 10*(i+j+k);
				}
			}
			else{
				app_LL_LP[i][j][0] = 0;
				app_LL_LP[i][j][1] = 0;
				for(int k =0 ; k<dim; k++){
					app_output[i][j][k] = 0;
				}
			}
		}
	}
	chain_temps[1]=5;
	output.append_to_output(app_output, app_LL_LP,app_positions);
	output.populate_chain_temperatures(chain_temps);
	
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Appended output"<<std::endl;
	for(int i = 0 ; i<chain_N; i++){
		for(int j =0 ; j<output.chain_lengths[i]; j++){
			for(int k =0 ; k<dim; k++){
				std::cout<<output.output[i][j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}

	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Appending new data to data dump"<<std::endl;
	output.append_to_data_dump("./data/testT.hdf5");
	output.append_to_data_dump("./data/test.hdf5");


	output.calc_ac_vals(true);
	for(int i = 0 ; i<output.cold_chain_number; i++){
		for(int j = 0 ; j<output.dimension; j++){
			std::cout<<output.ac_vals[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	output.count_indep_samples(false);
	std::cout<<"Indep samples "<<output.indep_samples<<std::endl;


	output.write_flat_thin_output( "./data/test_flat",false,false);
	//############################################################
	//Cleanup
	//############################################################
	for(int i = 0 ; i<chain_N; i ++){
		for(int j = 0 ; j<steps; j++){
			delete [] init_output[i][j];
			delete [] init_LL_LP[i][j];
		}
		delete [] init_output[i];
		delete [] init_LL_LP[i];
	}
	delete [] init_output;
	delete [] init_LL_LP;
	init_output = NULL;
	init_LL_LP = NULL;
	for(int i = 0 ; i<chain_N; i ++){
		for(int j = 0 ; j<app_steps; j++){
			delete [] app_output[i][j];
			delete [] app_LL_LP[i][j];
		}
		delete [] app_output[i];
		delete [] app_LL_LP[i];
	}
	delete [] app_output;
	delete [] app_LL_LP;
	app_output = NULL;
	app_LL_LP = NULL;
	
	return 0;
}
struct rosenbock_param
{
	double a;
	double mu;
	int n1;
	int n2;
	int n;
	double *b;
	double *bounds_high;
	double *bounds_low;
};
//https://arxiv.org/pdf/1903.09556.pdf
int mcmc_rosenbock(int argc, char *argv[])
{
	double a = 1./20.;
	//int n1 = 3;
	//int n2 = 2;
	int n1 = 5;
	int n2 = 3;
	int n = (n1-1)*n2 + 1;
	double b[n];
	for(int i = 0 ; i<n ; i++){
		b[i]=100./20.;
	}
	double mu = 1.;
	double bounds_high[n];
	double bounds_low[n];
	for(int i=0 ; i<n; i++){
		bounds_high[i]=10000;
		bounds_low[i]=-10000;
	}
	double out_param[5+n];
	out_param[0]=a;
	out_param[1]=mu;
	out_param[2]=n;
	out_param[3]=n1;
	out_param[4]=n2;
	for(int i = 0 ; i<n; i++){
		out_param[5+i]=b[i];
	}
	write_file("data/rosenbock_parameters.csv",out_param,5+n);



	double initial_pos[n];
	double seeding_var[n] ;
	for(int i = 0 ; i<n ; i++){
		initial_pos[i]=pow(-1,i)*500.;
		seeding_var[i]=10.;
	}
	int N_steps = 20000;
	int chain_N= 100;
	int max_chain_N= 50;
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 4;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	//chain_temps[0] = 1.;
	//double c = 1.8;
	//for(int i =1; i < chain_N/2;  i ++)
	//	chain_temps[i] =  chain_temps[i-1] * c;
	//chain_temps[chain_N/2] = 1;
	//for(int i =chain_N/2+1; i < chain_N;  i ++)
	//	chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output_RB.hdf5";
	std::string statfilename = "data/mcmc_statistics_RB.txt";
	std::string checkpointfile = "data/mcmc_checkpoint_RB.csv";
	//std::string LLfile = "data/mcmc_LL_RB.csv";
	std::string LLfile = "";
	
	int numThreads = 10;
	bool pool = true;
	bool show_progress = true;
	
	double **output;
	output = allocate_2D_array(  N_steps, n );
	int t0 = 50000;
	int nu = 100;
	double corr_threshold = 0.01;
	int corr_segments = 5;
	double corr_convergence_thresh = 0.01;
	double corr_target_ac = 0.01;
	std::string chain_distribution_scheme="double";
	int max_chunk_size = 1000000;
	
	rosenbock_param **param = new rosenbock_param*[chain_N];
	for(int i = 0 ; i<chain_N; i++){
		param[i]=new rosenbock_param;
		param[i]->a = a;
		param[i]->b = b;
		param[i]->bounds_high = bounds_high;
		param[i]->bounds_low = bounds_low;
		param[i]->mu = mu;
		param[i]->n1 = n1;
		param[i]->n2 = n2;
		param[i]->n = n;
	}
	
	mcmc_sampler_output sampler_output(chain_N,n);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, n, N_steps, chain_N, max_chain_N,initial_pos,seeding_var,chain_temps, swp_freq, t0,nu,corr_threshold, corr_segments, corr_convergence_thresh,corr_target_ac, max_chunk_size,chain_distribution_scheme, log_rosenbock_prior, log_rosenbock,NULL,(void **)param,numThreads, pool,show_progress, statfilename,chainfile, LLfile,checkpointfile );	
	sampler_output.calc_ac_vals(true);
	for(int i = 0 ; i<sampler_output.cold_chain_number; i++){
		std::cout<<sampler_output.max_acs[i]<<std::endl;
	}
	sampler_output.count_indep_samples(true);
	std::cout<<"Indep samples "<<sampler_output.indep_samples<<std::endl;
	sampler_output.create_data_dump(true,false, chainfile);
	sampler_output.create_data_dump(true,true,"data/mcmc_output_RB_trimmed.hdf5" );
	//sampler_output.write_flat_thin_output( "./data/test.hdf5",true,true);
	deallocate_2D_array(output, N_steps, n);
	for(int i = 0 ; i<chain_N; i++){
		delete param[i];	
	}
	delete param;	
	std::cout<<"ENDED"<<std::endl;


	return 0;
}

//void fisher_rosenbock(double *c,int dim, double **fisher,  void *parameters)
//{
//	double epsilon = 1e-3;
//	double temp[dim];
//	for (int i = 0 ; i<dim ; i++){
//		temp[i]=c[i];	
//	}
//	for (int i = 0 ; i<dim ; i++){
//		for (int j = 0 ; j<i ; j++){
//			temp[i] = c[i]+epsilon;
//			temp[j] = c[j]+epsilon;
//			double plusplus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]-epsilon;
//			temp[j] = c[j]-epsilon;
//			double minusminus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]+epsilon;
//			temp[j] = c[j]-epsilon;
//			double plusminus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i]-epsilon;
//			temp[j] = c[j]+epsilon;
//			double minusplus = log_rosenbock(temp,dim, parameters);
//			temp[i] = c[i];
//			temp[j] = c[j];
//			fisher[i][j]= -(plusplus -plusminus - minusplus + minusminus)/(4.*epsilon*epsilon);
//			//fisher[i][j]= 0;
//		}
//		temp[i]= c[i]+epsilon;
//		double plus = log_rosenbock(temp,dim,parameters);	
//		temp[i]= c[i]-epsilon;
//		double minus = log_rosenbock(temp,dim,parameters);	
//		temp[i]=c[i];
//		double eval = log_rosenbock(temp,dim,parameters);	
//		fisher[i][i]= -(plus -2*eval +minus)/(epsilon*epsilon);
//		//fisher[i][i]= 1;
//	}
//	for (int i = 0 ; i<dim ; i++){
//		for (int j = i+1 ; j<dim ; j++){
//			fisher[i][j] = fisher[j][i];	
//		}
//	}
//	
//	return;
//} 
double log_rosenbock (double *c,mcmc_data_interface *interface,void *parameters)
{
	int dim = interface->max_dim;
	rosenbock_param *param = (rosenbock_param *)parameters;
	double LL = 0;
	LL-= param->a * pow_int(c[0] - param->mu,2)	;
	for(int j = 0 ; j<param->n2; j++){
		for(int i = 1 ; i<param->n1; i++){
			if(i == 1){
				LL -= param->b[(j)*(param->n1-1) + i]*pow_int(c[(j)*(param->n1-1) + i] - c[0]*c[0],2);
			}
			else{
				LL -= param->b[(j)*(param->n1-1) + i]*pow_int(c[(j)*(param->n1-1) + i] - c[(j)*(param->n1-1) + i-1]*c[(j)*(param->n1-1) + i-1],2);

			}
		}
	}
	return LL;
	//return 2;
}
double log_rosenbock_prior (double *c,mcmc_data_interface *interface,void *parameters)
{
	int dim =interface->max_dim;
	double a = -std::numeric_limits<double>::infinity();
	rosenbock_param *param = (rosenbock_param *)parameters;
	for(int i = 0 ; i<param->n; i++){
		if( c[i]<param->bounds_low[i] || c[i] > param->bounds_high[i]){return a;}
	}
	return 0;
}
int mcmc_real_data(int argc, char *argv[])
{
	double gps = 1126259462.4;
	
	int detect_number = 2;
	std::string detectors[2] = {"Hanford","Livingston"};
	T_mcmc_gw_tool= 32;
	int dat_length = 65537;
	double **unpack = allocate_2D_array(dat_length, 4);
	read_file("data/data_H1.csv", unpack,dat_length , 4);
	std::complex<double> **data = new std::complex<double>*[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[dat_length-1];
	}
	double **freq = allocate_2D_array( 2,dat_length- 1);
	double **psd = allocate_2D_array( 2,dat_length- 1);
	for(int i = 1 ; i < dat_length; i++){
		data[0][i-1] = std::complex<double>(unpack[i][2],unpack[i][3]);
		freq[0][i-1] = unpack[i][0];
		psd[0][i-1] = unpack[i][1];
	}
	read_file("data/data_L1.csv", unpack,dat_length , 4);
	for(int i = 1 ; i < dat_length; i++){
		data[1][i-1] = std::complex<double>(unpack[i][2],unpack[i][3]);
		freq[1][i-1] = unpack[i][0];
		psd[1][i-1] = unpack[i][1];
	}

	

	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=dat_length-1;	
	}
	//#############################################################
	//#############################################################
	
	double **whitened = allocate_2D_array(data_lengths[0],7);
	for(int i = 0 ; i<data_lengths[0]; i++){
		whitened[i][0]= freq[0][i];
		whitened[i][1]=psd[0][i];
		whitened[i][2]=psd[1][i];
		whitened[i][3]=real(data[0][i]);
		whitened[i][4]=imag(data[0][i]);
		whitened[i][5]=real(data[1][i]);
		whitened[i][6]=imag(data[1][i]);
	}	
	write_file("data/whitened_data.csv",whitened,data_lengths[0],7);
	deallocate_2D_array(whitened, data_lengths[0],7);
	
	//#############################################################
	//#############################################################
	

	int dim = 11;
	std::string recovery_method = "IMRPhenomD";
	int chains = 20;
	double temps[chains];
	double c = 1.2;
	temps[0]=1;
	for(int i = 1 ; i<chains; i++){
		temps[i]=temps[i-1]*c;
	}
	
	double initial_position[dim]= {1.8, sin(-1.2),1, -.9, 1, 32*3./4, log(410),log(32), .24, .0,0};
	//double initial_position[dim]= {injection.RA, sin(injection.DEC), cos(injection.incl_angle), injection.phiRef, tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	write_file("data/injections.csv",initial_position,dim);
	//initial_position[3] = -initial_position[3];
	//initial_position[1] = 0;
	//initial_position[0] = 2.;
	double *seeding = NULL;
	int swap_freq = 3;
	int threads = 10;
	bool pool = true;
	bool show_progress = true;
	std::string stat_file = "data/experiment_stat.txt";
	std::string output_file = "data/experiment_output.csv";
	std::string ll_file = "data/experiment_ll.csv";
	std::string checkpoint_file = "data/experiment_checkpoint.csv";
	//std::string ac_file = "data/injection_ac.csv";
	std::string ac_file = "";
	MCMC_modification_struct mod_struct;
	
	int samples = 50000;
	double ***output  = allocate_3D_array(chains, samples, dim);
	PTMCMC_MH_GW(output, dim , samples, chains, initial_position, seeding, temps, swap_freq, standard_log_prior_D, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ac_file,ll_file, checkpoint_file);
	deallocate_3D_array(output, chains, samples, dim);
	//int samples = 2000;
	//double **output  = allocate_2D_array( samples, dim);
	//double t0 = 1000;
	//double nu = 500;
	//double corr_threshold = .01;
	//double corr_segments = 10;
	//double corr_converge_thresh = .1;
	//double corr_target_ac = .01;
	//std::string chain_distribution="double";
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dim , samples, chains, chains,initial_position, seeding, temps, swap_freq, t0,nu,corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_distribution,standard_log_prior_D, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, Nmod, bppe, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//deallocate_2D_array(output,  samples, dim);

	

	return 0;
}
int mcmc_injection(int argc, char *argv[])
{
	gen_params injection;
	injection.mass1 = 20;
	injection.mass2 = 1.4;
	//injection.mass1 = 9;
	//injection.mass2 = 7;
	double chirpmass = calculate_chirpmass(injection.mass1,injection.mass2);
	double eta = calculate_eta(injection.mass1,injection.mass2);
	//double chirpmass = 8.49;
	//double eta = .22;
	//injection.mass1 = calculate_mass1(chirpmass,eta);
	//injection.mass2 = calculate_mass2(chirpmass,eta);
	injection.Luminosity_Distance =300;
	injection.psi = .2;
	injection.phiRef = 2.;
	injection.f_ref = 20.;
	injection.RA = 1.5;
	injection.DEC = -1.2;
	//injection.RA = 1.5;
	//injection.DEC = 1.;
	injection.spin1[2] = .41;
	injection.spin2[2] = -.01;
	injection.spin1[1] = .01;
	injection.spin2[1] = -.01;
	injection.spin1[0] = .01;
	injection.spin2[0] = -.01;
	//IMRPhenomPv2<double> model;
	//injection.chip = model.PhenomPv2_inplane_spin(&injection);
	//injection.phip = 1;
	//injection.incl_angle = -M_PI +.01;
	injection.incl_angle = -M_PI +.51;
	double gps = 1126259462.4;
	//double gps = 1180922494.5;
	injection.gmst = gps_to_GMST_radian(gps);
	injection.equatorial_orientation = false;
	injection.horizon_coord = false;
	injection.shift_time = true;
	injection.shift_phase = true;


	//injection.Nmod_phi = 1;
	//injection.delta_phi = new double[1];
	//injection.delta_phi[0] = 1;
	//injection.phii = new int[1];
	//injection.phii[0] = 5;
	
	int detect_number = 3;
	std::string detectors[4] = {"Hanford","Livingston","Virgo","Kagra"};
	//std::string detectors[3] = {"Livingston","Hanford","Virgo"};
	//std::string SN[4] = {"Hanford_O1_fitted","Hanford_O1_fitted","Hanford_O1_fitted","KAGRA_pess"};
	std::string SN[4] = {"AdLIGODesign","AdLIGODesign","AdLIGODesign","KAGRA_pess"};
	std::string injection_method = "IMRPhenomPv2";
	double fmin = 5;
	double fmax =1024;
	T_mcmc_gw_tool= 4;
	double tc_ref = T_mcmc_gw_tool*(1-3./4.);
	double deltaf = 1./T_mcmc_gw_tool;
	int length = (fmax-fmin)/deltaf;
	std::complex<double> **data = new std::complex<double>*[detect_number];
	double **psd = new double*[detect_number];
	double **freq = new double*[detect_number];
	double total_snr = 0;
	for(int i = 0 ; i<detect_number; i++){
		data[i]= new std::complex<double>[length];
		psd[i]= new double[length];
		freq[i]= new double[length];
		for(int j =0 ; j<length; j++){
			freq[i][j] = fmin + j*deltaf;
		}
		populate_noise(freq[i],SN[i],psd[i],length);
		for(int j =0 ; j<length; j++){
			psd[i][j] *= psd[i][j];
		}
		double deltat = DTOA_DETECTOR(injection.RA,injection.DEC,injection.gmst,detectors[0],detectors[i]);
		//injection.tc = tc_ref + deltat;
		injection.tc = tc_ref - deltat;
		fourier_detector_response(freq[i],length, data[i],detectors[i],injection_method, &injection, (double *)NULL);
		total_snr += pow_int( calculate_snr_internal(psd[i],data[i],freq[i],length, "SIMPSONS",(double*) NULL, false), 2);
	}
	std::cout<<"NETWORK SNR of injection: "<<sqrt(total_snr)<<std::endl;


	int data_lengths[detect_number];
	for(int i = 0 ; i<detect_number; i++){
		data_lengths[i]=length;	
	}
	//###################################################################
	//###################################################################
	double fisher_dim = 11;
	double **jac_spins = allocate_2D_array(fisher_dim,fisher_dim);
	for (int i = 0 ;i<fisher_dim; i++){
		for(int j =0 ;j<fisher_dim; j++){
			if(i == 9 and j ==10){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==9){
				jac_spins[i][j] = .5;
			}
			else if(i == 10 and j ==10){
				jac_spins[i][j] = -.5;
			}
			else if(i == 9 and j ==9){
				jac_spins[i][j] =.5;
			}
			else if(i != j ){
				jac_spins[i][j] =0;
			}
			else {
				jac_spins[i][j] =1;
			}
		}
	}
	double **temp_fisher = allocate_2D_array(fisher_dim,fisher_dim);
	double **fisher = allocate_2D_array(fisher_dim,fisher_dim);
	double **cov = allocate_2D_array(fisher_dim,fisher_dim);
	for(int i = 0 ; i<fisher_dim; i++){
		for(int j = 0 ; j<fisher_dim; j++){
			fisher[i][j]=0;	
		}
	}
	for(int i = 0 ; i<detect_number; i++){
		fisher_autodiff(freq[i],length, "IMRPhenomD", detectors[i],detectors[0],temp_fisher,fisher_dim,&injection, "SIMPSONS",(double *)NULL,false, psd[i],NULL,NULL);
		for(int k = 0 ; k<fisher_dim; k++){
			for(int j = 0 ; j<fisher_dim; j++){
				fisher[k][j]+=temp_fisher[k][j];	
			}
		}
	}
	matrix_multiply(fisher, jac_spins,temp_fisher,fisher_dim,fisher_dim,fisher_dim);
	matrix_multiply( jac_spins,temp_fisher,fisher,fisher_dim,fisher_dim,fisher_dim);
	//for(int k = 0 ; k<fisher_dim; k++){
	//	for(int j = 0 ; j<fisher_dim; j++){
	//		std::cout<<fisher[k][j]<<std::endl;	
	//	}
	//}
	gsl_LU_matrix_invert(fisher, cov, fisher_dim);
	std::cout<<"Fisher estimates of covariances: "<<std::endl;
	for(int i = 0  ;i<fisher_dim; i++){
		std::cout<<i<<" "<<sqrt(cov[i][i])<<" "<<std::endl;
	}
	std::cout<<std::endl;
	deallocate_2D_array(temp_fisher,fisher_dim,fisher_dim);
	deallocate_2D_array(fisher,fisher_dim,fisher_dim);
	deallocate_2D_array(cov,fisher_dim,fisher_dim);
	//###################################################################
	//###################################################################
	//#############################################################
	//#############################################################
	
	double **whitened = allocate_2D_array(data_lengths[0],7);
	for(int i = 0 ; i<data_lengths[0]; i++){
		whitened[i][0]= freq[0][i];
		whitened[i][1]=psd[0][i];
		whitened[i][2]=psd[1][i];
		whitened[i][3]=real(data[0][i]);
		whitened[i][4]=imag(data[0][i]);
		whitened[i][5]=real(data[1][i]);
		whitened[i][6]=imag(data[1][i]);
	}	
	write_file("data/whitened_data.csv",whitened,data_lengths[0],7);
	deallocate_2D_array(whitened, data_lengths[0],7);
	
	//#############################################################
	//#############################################################
	

	int dim = 15;
	std::string recovery_method = "dCS_IMRPhenomPv2";
	int chains = 20;
	double temps[chains];
	double c = 1.2;
	temps[0]=1;
	for(int i = 1 ; i<chains; i++){
		temps[i]=temps[i-1]*c;
	}
	
	double spin1sph[3];
	double spin2sph[3];
	transform_cart_sph(injection.spin1,spin1sph);
	transform_cart_sph(injection.spin2,spin2sph);
	double initial_position[dim]= {injection.RA, sin(injection.DEC),injection.psi, cos(injection.incl_angle), injection.phiRef, T_mcmc_gw_tool-tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, spin1sph[0],spin2sph[0],cos(spin1sph[1]),cos(spin2sph[1]),fabs(spin1sph[2]-spin2sph[2]),0.01};
	//double initial_position[dim]= {injection.RA, sin(injection.DEC),injection.psi, cos(injection.incl_angle), injection.phiRef, T_mcmc_gw_tool-tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2],.01};
	//double initial_position[dim]= {injection.RA, sin(injection.DEC), cos(injection.incl_angle), injection.phiRef, tc_ref, log(injection.Luminosity_Distance),log(chirpmass), eta, injection.spin1[2],injection.spin2[2]};
	write_file("data/injections.csv",initial_position,dim);
	//initial_position[3] = -initial_position[3];
	//initial_position[1] = 0;
	//initial_position[0] = 2.;
	double *seeding = NULL;
	int swap_freq = 3;
	int threads = 10;
	bool pool = true;
	bool show_progress = true;
	std::string stat_file = "data/injection_stat.txt";
	std::string output_file = "data/injection_output.csv";
	std::string ll_file = "data/injection_ll.csv";
	std::string checkpoint_file = "data/injection_checkpoint.csv";
	//std::string ac_file = "data/injection_ac.csv";
	std::string ac_file = "";
	MCMC_modification_struct mod_struct;
	mod_struct.ppE_Nmod = 1;
	
	int samples = 50000;
	double ***output  = allocate_3D_array(chains, samples, dim);
	PTMCMC_MH_GW(output, dim , samples, chains, initial_position, seeding, temps, swap_freq, standard_log_prior_dCS, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ac_file,ll_file, checkpoint_file);
	deallocate_3D_array(output, chains, samples, dim);
	//double **output  = allocate_2D_array( samples, dim);
	//int max_chains_ensemble = chains;
	//int t0 = 1000;
	//int nu = 100;
	//MCMC_modification_struct mod_struct;
	//mod_struct.ppE_Nmod = 1;
	//int b = -1;
	//mod_struct.bppe =&b ;
	//double ac_target = .01;
	//int correlation_thresh = 10;
	//int correlation_segs =10 ;
	//double correlation_convergence =.1 ;
	//std::string allocation_scheme = "double";
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dim , samples, chains, max_chains_ensemble,initial_position, seeding, temps, swap_freq,t0,nu,correlation_thresh, correlation_segs,correlation_convergence, ac_target,allocation_scheme, standard_log_prior_dCS, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, &mod_struct, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//deallocate_2D_array(output,  samples, dim);
	//
	//int samples = 2000;
	//double **output  = allocate_2D_array( samples, dim);
	//double t0 = 1000;
	//double nu = 500;
	//double corr_threshold = .01;
	//double corr_segments = 10;
	//double corr_converge_thresh = .1;
	//double corr_target_ac = .01;
	//std::string chain_distribution="double";
	//PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(output, dim , samples, chains, chains,initial_position, seeding, temps, swap_freq, t0,nu,corr_threshold, corr_segments, corr_converge_thresh, corr_target_ac,chain_distribution,standard_log_prior_D, threads, pool, show_progress, detect_number, data, psd, freq, data_lengths, gps, detectors, Nmod, bppe, recovery_method, stat_file, output_file, ll_file, checkpoint_file);
	//deallocate_2D_array(output,  samples, dim);

	
	for(int i = 0 ; i<detect_number; i++){
		delete [] freq[i];	
		delete [] psd[i];	
		delete [] data[i];	
	}
	delete [] freq;
	delete [] psd;
	delete [] data;
	deallocate_2D_array(jac_spins,  fisher_dim, fisher_dim);

	return 0;
}
int mcmc_standard_test(int argc, char *argv[])
{
	int dimension = 2;
	double initial_pos[2]={1,0.};
	double *seeding_var = NULL;
	int N_steps = 20000;
	int chain_N= 50;
	int max_chain_N= 10;
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 2;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.8;
	for(int i =1; i < chain_N/2;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	chain_temps[chain_N/2] = 1;
	for(int i =chain_N/2+1; i < chain_N;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output.csv";
	std::string statfilename = "data/mcmc_statistics.txt";
	std::string checkpointfile = "data/mcmc_checkpoint.csv";
	//std::string LLfile = "data/mcmc_LL.csv";
	std::string LLfile = "";
	
	int numThreads = 10;
	bool pool = true;
	bool show_progress = true;
	
	//double ***output;
	//output = allocate_3D_array( chain_N, N_steps, dimension );
	//PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	//deallocate_3D_array(output, chain_N, N_steps, dimension);
	double **output;
	output = allocate_2D_array(  N_steps, dimension );
	int t0 = 5000;
	int nu = 100;
	double corr_threshold = 0.01;
	int corr_segments = 5;
	double corr_convergence_thresh = 0.01;
	double corr_target_ac = 0.01;
	std::string chain_distribution_scheme="double";
	int max_chunk_size = 1000000;
	mcmc_sampler_output sampler_output(chain_N,dimension);
	PTMCMC_MH_dynamic_PT_alloc_uncorrelated(&sampler_output,output, dimension, N_steps, chain_N, max_chain_N,initial_pos,seeding_var,chain_temps, swp_freq, t0,nu,corr_threshold, corr_segments, corr_convergence_thresh,corr_target_ac, max_chunk_size,chain_distribution_scheme, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile, LLfile,checkpointfile );	
	deallocate_2D_array(output, N_steps, dimension);
	sampler_output.write_flat_thin_output( "./data/test_flat_standard",true,true);
	std::cout<<"ENDED"<<std::endl;


		
	return 0;

}
double log_test (double *c,mcmc_data_interface *interface,void *parameters)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_test_prior (double *c,mcmc_data_interface *interface,void *parameters)
{
	return 1.;
}

double fisher11(double *c)
{
	double x = c[0];
	double y = c[1];

	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(x,2) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(-2 - 128*pow(x,2) - 16*(9 + 4*pow(x,2) + 8*y)) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y),2));
}
double fisher12(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*x + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*x*(-2 + y) - 
     16*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(9 + 4*pow(x,2) + 8*y)*(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y)));

}
double fisher22(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) - 128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2)) + 
     128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(-2 + y,2) + 256*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(9 + 4*pow(x,2) + 8*y,2));

}
void fisher_test(double *c, double **fisher, mcmc_data_interface *interface, void *parameters)
{
	
	fisher[0][0] = fisher11(c);
	fisher[1][1] = fisher22(c);
	fisher[1][0] = fisher12(c);
	fisher[0][1] = fisher12(c);
	return;
} 
	
double standard_log_prior_D(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>60 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	double chi1 = pos[9]+pos[10];
	double chi2 = pos[9]-pos[10];
	if ((chi1)<-.95 || (chi1)>.95){return a;}//chi1 
	if ((chi2)<-.95 || (chi2)>.95){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}

}
double standard_log_prior_dCS(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	if (std::exp(pos[7])<.1 || std::exp(pos[7])>20 ){return a;}//chirpmass
	if ((pos[8])<.01 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chi2
	if ((pos[12])<-1 || (pos[12])>1){return a;}//chi2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//chi2
	if(pos[14]<0 || pos[14] > 100) { return a;}
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}
}
double standard_log_prior_P(double *pos, mcmc_data_interface *interface,void *parameters)
{
	double chirp = std::exp(pos[7]);
	double eta = pos[8];
	double a = -std::numeric_limits<double>::infinity();
	if ((pos[0])<0 || (pos[0])>2*M_PI){ return a;}//RA

	if ((pos[1])<-1 || (pos[1])>1){return a;}//sinDEC
	//if ((pos[1])<-M_PI/2 || (pos[1])>M_PI/2){return a;}//sinDEC

	if ((pos[2])<0 || (pos[2])>M_PI){return a;}//PSI
	if ((pos[3])<-1 || (pos[3])>1){return a;}//cos \iota
	if ((pos[4])<0 || (pos[4])>2*M_PI){return a;}//phiRef
	//if ((pos[5])<T_mcmc_gw_tool*3./4.-.1 || (pos[5])>T_mcmc_gw_tool*3./4.+.1){return a;}//tc
	if ((pos[5])<0 || (pos[5])>T_mcmc_gw_tool){return a;}//tc
	if (std::exp(pos[6])<10 || std::exp(pos[6])>1000){return a;}//DL
	if (std::exp(pos[7])<2 || std::exp(pos[7])>60 ){return a;}//chirpmass
	if ((pos[8])<.1 || (pos[8])>.249999){return a;}//eta
	if ((pos[9])<0 || (pos[9])>.95){return a;}//chi1 
	if ((pos[10])<0 || (pos[10])>.95){return a;}//chi2
	if ((pos[11])<-1 || (pos[11])>1){return a;}//chi2
	if ((pos[12])<-1 || (pos[12])>1){return a;}//chi2
	if ((pos[13])<0 || (pos[13])>2*M_PI){return a;}//chi2
	else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] ;}
	//else {return log(chirpmass_eta_jac(chirp,eta))+3*pos[6] -log(cos(asin(pos[1]))) ;}
}
//Uniform in m1 and m2, transformed to lnM and eta
double chirpmass_eta_jac(double chirpmass, double eta){
	return chirpmass*chirpmass/(sqrt(1. - 4.*eta)*pow(eta,1.2));
}

void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Standard"<<std::endl;
	std::cout<<"1 --- GW injection"<<std::endl;
	std::cout<<"2 --- GW experimentation"<<std::endl;
	std::cout<<"3 --- N-dimensional Rosenbock"<<std::endl;
	std::cout<<"4 --- MCMC output class test"<<std::endl;
}
