#include "psd_estimation.h"
#include "mcmc_sampler.h"
#include "util.h"
#include "io_util.h"
#include <iostream>
#include <complex>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <bits/stdc++.h>


/*Custom implementation of BayesLine (https://arxiv.org/pdf/1410.3852.pdf)
 */
void bayesline_psd_estimation(
	std::complex<double> *data_stream , 
	double *frequencies,
	int length,
	double T,
	int N_L_MIN, 
	int N_L_MAX, 
	int N_S_MIN, 
	int N_S_MAX,
	double *initial_pos,
	int *initial_status,
	double *seeding_var,
	std::string chain_allocation_scheme,
	int samples,	
	int chain_N,
	int max_chain_ensemble,
	double *chain_temps,
	int swap_freq,
	int t0,
	int nu,	
	int max_chunk_size,
	int threads,
	bool pool,
	bool show_prog,
	std::string chain_file,
	std::string stat_file,
	std::string checkpoint_file,
	PSD_output *output
	)
{
	bayesline_sampling_struct **parameters = new bayesline_sampling_struct*[chain_N];
	for(int i = 0 ; i<chain_N; i++){
		parameters[i] = new bayesline_sampling_struct;
		parameters[i]->N_L_MIN = N_L_MIN;
		parameters[i]->N_L_MAX = N_L_MAX;
		parameters[i]->N_S_MIN = N_S_MIN;
		parameters[i]->N_S_MAX = N_S_MAX;
		parameters[i]->data = data_stream;
		parameters[i]->data_mod_sq = new double[length];
		for(int j = 0 ; j<length; j++){
			parameters[i]->data_mod_sq[j] = pow_int(std::abs(data_stream[j]),2) ;
		}
		parameters[i]->signal_length = T;
		parameters[i]->data_length = length;
		parameters[i]->frequencies = frequencies;
		parameters[i]->deltaf_factor = 50;
		parameters[i]->MIN_FREQ = frequencies[0];
		parameters[i]->MAX_FREQ = frequencies[length-1];
		parameters[i]->MAX_SN_AMP = 1e-30;
		parameters[i]->MIN_SN_AMP = 1e-50;
		parameters[i]->MAX_LOGSN_AMP = log(parameters[i]->MAX_SN_AMP);
		parameters[i]->MIN_LOGSN_AMP = log(parameters[i]->MIN_SN_AMP);
		parameters[i]->MAX_L_AMP = 1e-20;
		parameters[i]->MIN_L_AMP = 1e-50;
		parameters[i]->MAX_LOGL_AMP = log(parameters[i]->MAX_L_AMP);
		parameters[i]->MIN_LOGL_AMP = log(parameters[i]->MIN_L_AMP);
		parameters[i]->MAX_Q = 1e15;
		parameters[i]->MIN_Q = 1e1;
	
	}
	

	int max_dim = N_L_MAX*3 + N_S_MAX*2+2;
	int min_dim = 2;
	double **temp_output = allocate_2D_array(samples, max_dim);
	int **parameter_status = allocate_2D_array_int(samples, max_dim);
	int *model_status = new int[samples];
	int nested_model_number = 0;
	int initial_model_status=0;
	if(!initial_status && !initial_pos){
		debugger_print(__FILE__,__LINE__,"Using generic initial position");
		initial_pos = new double[max_dim];
		initial_status= new int[max_dim];
		double step_size = 10;
		//int initial_spline_locs = int((parameters[0]->MAX_FREQ - parameters[0]->MIN_FREQ)/10);
		int initial_spline_locs = int((parameters[0]->MAX_FREQ - parameters[0]->MIN_FREQ)/step_size);
		//int initial_lorentzian_locs = 10;
		int initial_lorentzian_locs = 3;
		debugger_print(__FILE__,__LINE__,"INIT spline knots");	
		debugger_print(__FILE__,__LINE__,initial_spline_locs);	
		for(int i = 0 ; i<max_dim; i++){
			initial_status[i] = 0;	
			initial_pos[i] = 0;	
		}
		double mean_PSD;
		mean_list(parameters[0]->data_mod_sq, 100,&mean_PSD);	
		initial_pos[0] = log(mean_PSD);	
		mean_list(&(parameters[0]->data_mod_sq[parameters[0]->data_length-101]), 100,&mean_PSD);	
		initial_pos[1] = log(mean_PSD);	
		initial_status[0] = 1;
		initial_status[1] = 1;
		int ID_step = length/(initial_spline_locs-1);
		std::cout<<ID_step<<std::endl;
		//int window = 100;
		for(int i = 0 ; i<initial_spline_locs; i++){
			initial_pos[i*2+2] = parameters[0]->MIN_FREQ*(1.05)+step_size*i;
			if(i == 0){
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step]), 20,&mean_PSD);	
			}
			else if (i == initial_spline_locs-1){
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step-20]), 20,&mean_PSD);	

			}
			else{
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step-10]), 20,&mean_PSD);	
			}
			initial_pos[i*2+1+2] = log(mean_PSD);
			//debugger_print(__FILE__,__LINE__,i);
			//debugger_print(__FILE__,__LINE__,i*ID_step);
			//std::cout<<mean_PSD<<std::endl;
			//std::cout<<log(mean_PSD)<<std::endl;
			//std::cout<<initial_pos[i*2+2]<<std::endl;
			initial_status[i*2+2 ] = 1;
			initial_status[i*2+1+2 ] = 1;
		}
		write_file("data/initial_position.csv",initial_pos,max_dim);
		
		ID_step = length/(initial_lorentzian_locs-1);
		std::cout<<ID_step<<std::endl;
		for(int i = 0 ; i<initial_lorentzian_locs; i++){
			if(i == 0){
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step]), 20,&mean_PSD);	
			}
			else if (i == initial_lorentzian_locs-1){
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step-20]), 20,&mean_PSD);	

			}
			else{
				mean_list(&(parameters[0]->data_mod_sq[i*ID_step-10]), 20,&mean_PSD);	
			}
			initial_pos[2*N_S_MAX+2+i*3] = log(mean_PSD);
			initial_pos[2*N_S_MAX+2+i*3+1] = parameters[0]->MIN_FREQ*(1.05)+int((parameters[0]->MAX_FREQ - parameters[0]->MIN_FREQ)/initial_lorentzian_locs)*i;
			initial_pos[2*N_S_MAX+2+i*3+2] = 1e5;
			initial_status[2*N_S_MAX+2+i*3 ] = 1;
			initial_status[2*N_S_MAX+2+i*3+1 ] = 1;
			initial_status[2*N_S_MAX+2+i*3+2 ] = 1;
		}
	}
	

	mcmc_sampler_output sampler_output(chain_N, max_dim, nested_model_number);
	sampler_output.RJ= true;

	RJPTMCMC_MH_dynamic_PT_alloc_comprehensive(&sampler_output, temp_output, parameter_status, model_status, nested_model_number, max_dim, min_dim, samples, chain_N, max_chain_ensemble, initial_pos, initial_status, initial_model_status, seeding_var, (double **)NULL, (int **)NULL, (int*)NULL, chain_temps, swap_freq, t0, nu, max_chunk_size, chain_allocation_scheme, bayesline_prior, bayesline_likelihood, NULL, bayesline_RJ_proposal, (void**)parameters, threads, pool, show_prog, true, stat_file, chain_file, "", checkpoint_file);

	const gsl_rng_type *gslT = gsl_rng_default;	
	gsl_rng *r =gsl_rng_alloc(gslT);
	int average_num = 1000;
	for(int j = 0 ; j<length ; j++){
		output->SN[j] = 0;	
	}
	
	for (int i =0 ; i<average_num; i++){
		int alpha = int((gsl_rng_uniform(r)*chain_N)/max_chain_ensemble);//Pick a cold chain
		int chain = alpha * max_chain_ensemble;
		int step =int(gsl_rng_uniform(r)*(samples));// Pick a step
	
		double *pos= sampler_output.output[alpha][step];	
		int *status= sampler_output.status[alpha][step];	
		double SN[length];
		model_PSD(pos, status, 0, parameters[0], SN);
		for(int j = 0 ; j<length ; j++){
			output->SN[j] += SN[j];
		}
	
	}
	for(int j = 0 ; j<length ; j++){
		output->SN[j] /= average_num;
	}
	
	gsl_rng_free(r);

	deallocate_2D_array(temp_output, samples, max_dim);
	deallocate_2D_array(parameter_status, samples, max_dim);
	delete [] model_status;
	for(int i = 0 ; i<chain_N; i++){
		delete [] parameters[i]->data_mod_sq;	
		delete parameters[i];	
	}
	delete [] parameters;
}

double lorentzian(double pt, double amp, double q, double f,double deltaf){
	double z= 1;
	double diff = fabs(f-pt);
	if (diff > deltaf){
		z = std::exp(-(diff - deltaf) / (deltaf));
	}
	double lorentz = z* amp* pow_int(pt,4);
	lorentz/= ( pow_int( pt * f, 2) + q*q * pow_int(pt* pt - f * f, 2));
	return lorentz;
}
void lorentzian_component(double *pos, int NL, bayesline_sampling_struct *p, double *frequencies, int length, double *SN)
{
	double pts[NL];
	double q_facs[NL];
	double amps[NL];
	for(int i =0  ; i<NL*3; i+=3){
		amps[i/3] = exp(pos[i]);
		pts[i/3] = pos[i+1];
		q_facs[i/3] = pos[i+2];
	}
	for(int i = 0 ; i<length; i++){
		SN[i] = 0 ;
		for(int j=0; j<NL; j++){
			SN[i] += lorentzian(pts[j],amps[j],q_facs[j], frequencies[i],pts[j]/p->deltaf_factor);
			//SN[i] += 0;
		}
	}
}

void order_list(double *pos, int NS, double *pts, double *SNs)
{
	std::pair<double, double> pairs[NS];	
	for(int i= 0; i<NS*2; i+=2){
		pairs[i/2].first = pos[i+2];
		pairs[i/2].second = pos[i+1+2];
	}
	sort(pairs, pairs+NS);
	
	for(int i= 0; i<NS; i++){
		pts[i+1] = pairs[i].first;
		SNs[i+1] =pairs[i].second;
		//std::cout<<exp(pts[i+1])<<" "<<exp(SNs[i+1])<<std::endl;
	}

}
void smooth_component(double *pos,int NS, bayesline_sampling_struct *p,double *frequencies, int length, double *SN )
{
	double pts[NS+2];
	double SNs[NS+2];
	order_list(pos, NS, pts, SNs);
	pts[0] = frequencies[0];
	pts[NS+1] = frequencies[length-1];
	SNs[0] = pos[0];
	SNs[NS+1] = pos[1];
	

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	//gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,NS+2);
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,NS+2);
	gsl_spline_init(spline, pts, SNs, NS+2);
	
	for(int i = 0 ; i<length; i++){
		SN[i] = exp(gsl_spline_eval(spline, frequencies[i],acc));
		//std::cout<<frequencies[i]<<" "<<SN[i]<<std::endl;
	}
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}

void model_PSD(double *pos, int *status, int model, bayesline_sampling_struct *p, double *SN)
{

	int NS = 0 ;
	int NL = 0 ;
	for (int i = 2 ; i<2*p->N_S_MAX+2; i+=2){
		if(status[i] !=0){
			NS++;
		}
	}
	for (int i = 0 ; i<3*p->N_L_MAX; i+=3){
		if(status[p->N_S_MAX*2+2+i] !=0){
			NL++;
		}
	}
	double SN_S[p->data_length];
	double SN_L[p->data_length];
	for(int i = 0 ; i<p->data_length; i++){
		SN_S[i] = 0;
		SN_L[i] = 0;
	}
	smooth_component(pos, NS, p, p->frequencies, p->data_length, SN_S);
	lorentzian_component(&(pos[p->N_S_MAX*2+2]), NL, p, p->frequencies, p->data_length, SN_L);
	for(int i = 0; i<p->data_length; i++){
		SN[i] = (SN_S[i] + SN_L[i]);
		//SN[i] = (SN_S[i] );
	}

	return;
}

double bayesline_likelihood(double *pos, int *status, int model, mcmc_data_interface *interface, void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	bayesline_sampling_struct *p = (bayesline_sampling_struct *)parameters;
	double SN[p->data_length];
	model_PSD(pos, status, model, p, SN);
	double ll=0;
	for(int i = 0; i<p->data_length; i++){
		ll+= p->data_mod_sq[i] / SN[i];
		//std::cout<<p->frequencies[i]<<" "<<SN[i]<<" "<<p->data_mod_sq[i]<<std::endl;
	}
	ll*=(-2./p->signal_length );
	//debugger_print(__FILE__,__LINE__,ll);
	if(isnan(ll)){
		return a;
	}
	return ll;
}

double bayesline_prior(double *pos, int *status, int model, mcmc_data_interface *interface, void *parameters)
{
	double a = -std::numeric_limits<double>::infinity();
	bayesline_sampling_struct *p = (bayesline_sampling_struct *)parameters;
	int NL=0,NS=0;
	for (int i = 2 ; i<2*p->N_S_MAX; i+=2){
		if(status[i] !=0){
			NS++;
		}
	}
	for (int i = p->N_S_MAX*2+2 ; i<3*p->N_L_MAX; i+=3){
		if(status[i] !=0){
			NL++;
		}
	}
	double lp=1;
	for (int i = 0 ; i<2; i++){
		if(pos[i] < p->MIN_LOGSN_AMP|| pos[i] > p->MAX_LOGSN_AMP){return a;}
		lp-=log((p->MAX_SN_AMP - p->MIN_SN_AMP));
		lp-=pos[i];
	}
	for (int i = 2 ; i<2*NS+2; i+=2){
		if( pos[i] < p->MIN_FREQ|| pos[i] > p->MAX_FREQ){ return a;}
		else if(pos[i+1] < p->MIN_LOGSN_AMP|| pos[i+1] > p->MAX_LOGSN_AMP){return a;}
		
		lp-=log((p->MAX_FREQ - p->MIN_FREQ));
		lp-=log((p->MAX_SN_AMP - p->MIN_SN_AMP));
		lp-=pos[i+1];
	}
	for (int i = 0 ; i<3*NL; i+=3){
		if(pos[2*p->N_S_MAX+2+i] < p->MIN_LOGL_AMP|| pos[2*p->N_S_MAX+2+i] > p->MAX_LOGL_AMP){return a;}
		else if( pos[2*p->N_S_MAX+2+i+1] < p->MIN_FREQ|| pos[2*p->N_S_MAX+2+i+1] > p->MAX_FREQ){ return a;}
		else if( pos[2*p->N_S_MAX+2+i+2] < p->MIN_Q|| pos[2*p->N_S_MAX+2+i+2] > p->MAX_Q){ return a;}
		lp-=log((p->MAX_FREQ - p->MIN_FREQ));
		lp-=log((p->MAX_L_AMP - p->MIN_L_AMP));
		lp-=log((p->MAX_Q - p->MIN_Q));
		lp-=pos[2*p->N_S_MAX+2+i];
	}
	return lp;
}

void bayesline_RJ_proposal(double *current_pos, double *prop_pos,int *current_status, int *prop_status,int *current_model, int *prop_model, double *MH_correction,mcmc_data_interface *interface, void *parameters)
{
	for(int i = 0 ; i<interface->max_dim; i++){
		prop_pos[i] = current_pos[i];
		prop_status[i] = current_status[i];
	}
	*prop_model = *current_model;
	bayesline_sampling_struct *p = (bayesline_sampling_struct *)parameters;
	return;
}

