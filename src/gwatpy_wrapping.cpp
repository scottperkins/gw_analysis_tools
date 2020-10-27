#include "gwatpy_wrapping.h"
#include "util.h"
#include "fisher.h"
#include "waveform_generator.h"
#include "waveform_util.h"
#include "mcmc_gw.h"
#include "mcmc_sampler.h"
#include "mcmc_sampler_internals.h"
#include "fisher.h"
#include <bits/stdc++.h> 




mcmc_data_interface * mcmc_data_interface_py(
	int min_dim,
	int max_dim,
	int chain_id,
	int nested_model_number,
	int chain_number,
	double RJ_step_width,
	bool burn_phase)
{
	
	mcmc_data_interface *interface = new mcmc_data_interface;
	interface->min_dim = min_dim;
	interface->max_dim = max_dim;
	interface->chain_id = chain_id;
	interface->nested_model_number = nested_model_number;
	interface->chain_number = chain_number;
	interface->RJ_step_width = chain_number;
	interface->burn_phase = burn_phase;
	return interface;

}
void mcmc_data_interface_destructor_py(mcmc_data_interface *interface)
{
	delete interface;
}

void MCMC_modification_struct_py_destructor(MCMC_modification_struct *mod_struct)
{
	if(mod_struct->bppe){delete [] mod_struct->bppe;mod_struct->bppe = NULL;}
	if(mod_struct->gIMR_phii){delete [] mod_struct->gIMR_phii;mod_struct->gIMR_phii = NULL;}
	if(mod_struct->gIMR_sigmai){delete [] mod_struct->gIMR_sigmai;mod_struct->gIMR_sigmai = NULL;}
	if(mod_struct->gIMR_betai){delete [] mod_struct->gIMR_betai;mod_struct->gIMR_betai = NULL;}
	if(mod_struct->gIMR_alphai){delete [] mod_struct->gIMR_alphai;mod_struct->gIMR_alphai = NULL;}
}
MCMC_modification_struct * MCMC_modification_struct_py( 
	int ppE_Nmod, 
	double *bppe,
	int gIMR_Nmod_phi,
	int *gIMR_phii,
	int gIMR_Nmod_sigma,
	int *gIMR_sigmai,
	int gIMR_Nmod_beta,
	int *gIMR_betai,
	int gIMR_Nmod_alpha,
	int *gIMR_alphai,
	bool NSflag1,
	bool NSflag2
	)
{
	MCMC_modification_struct *mod_struct = new MCMC_modification_struct;
	mod_struct->ppE_Nmod = ppE_Nmod;
	mod_struct->bppe = NULL;
	if(mod_struct->ppE_Nmod !=0){
		mod_struct->bppe = new double[mod_struct->ppE_Nmod];
	}
	mod_struct->gIMR_Nmod_phi = gIMR_Nmod_phi;
	mod_struct->gIMR_phii = NULL;
	if(mod_struct->gIMR_Nmod_phi !=0){
		mod_struct->gIMR_phii = new int[mod_struct->gIMR_Nmod_phi];
		for(int i = 0 ;i<mod_struct->gIMR_Nmod_phi; i++){
			mod_struct->gIMR_phii[i] = gIMR_phii[i];
		}
	}
	mod_struct->gIMR_Nmod_sigma = gIMR_Nmod_sigma;
	mod_struct->gIMR_sigmai = NULL;
	if(mod_struct->gIMR_Nmod_sigma !=0){
		mod_struct->gIMR_sigmai = new int[mod_struct->gIMR_Nmod_sigma];
		for(int i = 0 ;i<mod_struct->gIMR_Nmod_sigma; i++){
			mod_struct->gIMR_sigmai[i] = gIMR_sigmai[i];
		}
	}
	mod_struct->gIMR_Nmod_beta = gIMR_Nmod_beta;
	mod_struct->gIMR_betai = NULL;
	if(mod_struct->gIMR_Nmod_beta !=0){
		mod_struct->gIMR_betai = new int[mod_struct->gIMR_Nmod_beta];
		for(int i = 0 ;i<mod_struct->gIMR_Nmod_beta; i++){
			mod_struct->gIMR_betai[i] = gIMR_betai[i];
		}
	}
	mod_struct->gIMR_Nmod_alpha = gIMR_Nmod_alpha;
	mod_struct->gIMR_alphai = NULL;
	if(mod_struct->gIMR_Nmod_alpha !=0){
		mod_struct->gIMR_alphai = new int[mod_struct->gIMR_Nmod_alpha];
		for(int i = 0 ;i<mod_struct->gIMR_Nmod_alpha; i++){
			mod_struct->gIMR_alphai[i] = gIMR_alphai[i];
		}
	}
	mod_struct->NSflag1 = NSflag1;
	mod_struct->NSflag2 = NSflag2;
	return mod_struct;
}


void pack_local_mod_structure_py(
	mcmc_data_interface *interface, 
	double *param, 
	int *status, 
	char *waveform_extended,
	//void * parameters, 
	MCMC_modification_struct *full_struct,
	MCMC_modification_struct *local_struct)
{

	//pack_local_mod_structure(interface, param, status, std::string(waveform_extended),parameters, full_struct, local_struct);
	pack_local_mod_structure(interface, param, status, std::string(waveform_extended),(void *)NULL, full_struct, local_struct);

	return;
}


char * MCMC_prep_params_py(
	double *param, 
	double *temp_params, 
	gen_params_base<double> *gen_params, 
	int dimension, 
	char * generation_method, 
	MCMC_modification_struct *mod_struct,
	bool save_gmst)
{
	double gmst;
	if (save_gmst){
		gmst = gen_params->gmst;
	}
	std::string str =  MCMC_prep_params(param, temp_params, gen_params, dimension, std::string(generation_method), mod_struct);
	char *char_arr = new char[  str.length()];
	strcpy(char_arr,str.c_str());

	if (save_gmst){
		gen_params->gmst=gmst;
	}

	return char_arr;
}

void repack_parameters_py(double *parameters, gen_params_base<double> *gen_param, char * generation_method, int dim )
{
	repack_parameters(parameters, gen_param, std::string(generation_method), dim , (gen_params_base<double> *)NULL);

}

int fourier_detector_response_py(double *frequencies,
	int length, 
	double  *response_real,
	double  *response_imaginary,
	char * detector,
	char *generation_method, 
	gen_params_base<double> *parameters)
{
	std::string gen_meth(generation_method);
	std::string det(detector);
	std::complex<double> *response = new std::complex<double>[length];
	int status =  fourier_detector_response(frequencies, length, response,det, gen_meth, parameters,(double*)NULL);
	for(int i = 0 ; i<length; i++){
		response_real[i] = std::real(response[i]);
		response_imaginary[i] = std::imag(response[i]);
	}
	delete [] response;	
	return status;
}

int fourier_waveform_py(double *frequencies,
	int length, 
	double  *wf_plus_real,
	double  *wf_plus_imaginary,
	double  *wf_cross_real,
	double  *wf_cross_imaginary,
	char *generation_method, 
	gen_params_base<double> *parameters)
{
	std::string gen_meth(generation_method);
	std::complex<double> *wf_p = new std::complex<double>[length];
	std::complex<double> *wf_c = new std::complex<double>[length];
	int status =  fourier_waveform(frequencies, length, wf_p, wf_c, gen_meth, parameters);
	for(int i = 0 ; i<length; i++){
		wf_plus_real[i] = std::real(wf_p[i]);
		wf_cross_real[i] = std::real(wf_c[i]);
		wf_plus_imaginary[i] = std::imag(wf_p[i]);
		wf_cross_imaginary[i] = std::imag(wf_c[i]);
	}
	delete [] wf_p;	
	delete [] wf_c;	
	return status;
}


gen_params_base<double>* gen_params_base_py(
	double mass1, 
	double mass2,
	double *spin1,
	double *spin2,
	double Luminosity_Distance,
	double incl_angle,
	double RA,
	double DEC,
	double psi,
	double gmst,
	double f_ref,
	double theta_l,
	double phi_l,
	double theta,
	double phi,
	char * cosmology,
	bool equatorial_orientation,
	bool horizon_coord,
	bool NSflag1,
	bool NSflag2,
	bool dep_postmerger,
	bool shift_time,
	bool shift_phase,
	bool sky_average,
	double LISA_alpha0,
	double LISA_phi0,
	int Nmod_phi,
	int Nmod_sigma,
	int Nmod_beta,
	int Nmod_alpha,
	int *phii,
	int *sigmai,
	int *betai,
	int *alphai,
	double *delta_phi,
	double *delta_sigma,
	double *delta_beta,
	double *delta_alpha,
	int Nmod,
	double *bppe,
	double *betappe
 )
{
	gen_params_base<double> *p = new gen_params_base<double> ;
	p->mass1 = mass1;	
	p->mass2 = mass2;	
	p->spin1[0] = spin1[0];	
	p->spin1[1] = spin1[1];	
	p->spin1[2] = spin1[2];	
	p->spin2[0] = spin2[0];	
	p->spin2[1] = spin2[1];	
	p->spin2[2] = spin2[2];	
	p->Luminosity_Distance = Luminosity_Distance;	
	p->incl_angle = incl_angle;
	p->RA = RA;
	p->DEC = DEC;
	p->psi = psi;
	p->f_ref = f_ref;
	p->theta_l = theta_l;
	p->phi_l = phi_l;
	p->theta = theta;
	p->phi = phi;
	p->cosmology = std::string(cosmology);
	p->equatorial_orientation=equatorial_orientation;
	p->horizon_coord=horizon_coord;
	p->gmst=gmst;
	p->NSflag1=NSflag1;
	p->NSflag2=NSflag2;
	p->dep_postmerger=dep_postmerger;
	p->shift_time=shift_time;
	p->shift_phase=shift_phase;
	p->sky_average=sky_average;
	p->LISA_alpha0=LISA_alpha0;
	p->LISA_phi0=LISA_phi0;
	p->Nmod_phi = Nmod_phi;
	p->Nmod_sigma = Nmod_sigma;
	p->Nmod_beta = Nmod_beta;
	p->Nmod_alpha = Nmod_alpha;
	if(p->Nmod_phi != 0){
		p->phii = new int[p->Nmod_phi];
		p->delta_phi = new double[p->Nmod_phi];
		for(int i = 0 ; i<p->Nmod_phi;i++){
			p->phii[i] = phii[i];
			p->delta_phi[i] = delta_phi[i];
		}
	}
	else{
		p->phii = NULL;
		p->delta_phi = NULL;
	}
	if(p->Nmod_sigma !=0){
		p->sigmai = new int[p->Nmod_sigma];
		p->delta_sigma = new double[p->Nmod_sigma];
		for(int i = 0 ; i<p->Nmod_sigma;i++){
			p->sigmai[i] = sigmai[i];
			p->delta_sigma[i] = delta_sigma[i];
		}
	}
	else{
		p->sigmai = NULL;
		p->delta_sigma = NULL;
	}
	if(p->Nmod_beta !=0){
		p->betai = new int[p->Nmod_beta];
		p->delta_beta = new double[p->Nmod_beta];
		for(int i = 0 ; i<p->Nmod_beta;i++){
			p->betai[i] = betai[i];
			p->delta_beta[i] = delta_beta[i];
		}
	}
	else{
		p->betai = NULL;
		p->delta_beta = NULL;
	}
	if(p->Nmod_alpha!=0){
		p->alphai = new int[p->Nmod_alpha];
		p->delta_alpha = new double[p->Nmod_alpha];
		for(int i = 0 ; i<p->Nmod_alpha;i++){
			p->alphai[i] = alphai[i];
			p->delta_alpha[i] = delta_alpha[i];
		}
	}
	else{
		p->alphai = NULL;
		p->delta_alpha = NULL;
	}
	p->Nmod = Nmod;
	if(Nmod != 0){
		p->bppe = new double[p->Nmod];
		p->betappe = new double[p->Nmod];
		for(int i = 0 ; i<p->Nmod;i++){
			p->bppe[i] = bppe[i];
			p->betappe[i] = betappe[i];
		}
	}
	else{
		p->betappe = NULL;
		p->bppe = NULL;

	}
	

	return p;
}

void gen_params_base_py_destructor(gen_params_base<double> *p)
{
	if(p->betappe){delete [] p->betappe;p->betappe=NULL;}
	//if(p->bppe){delete [] p->bppe;p->bppe=NULL;}
	if(p->delta_phi){delete [] p->delta_phi;p->delta_phi=NULL;}
	//if(p->phii){delete [] p->phii;p->phii=NULL;}
	//if(p->sigmai){delete [] p->sigmai;p->sigmai=NULL;}
	if(p->delta_sigma){delete [] p->delta_sigma;p->delta_sigma=NULL;}
	//if(p->betai){delete [] p->betai;p->betai=NULL;}
	if(p->delta_beta){delete [] p->delta_beta;p->delta_beta=NULL;}
	//if(p->alphai){delete [] p->alphai;p->alphai=NULL;}
	if(p->delta_alpha){delete [] p->delta_alpha;p->delta_alpha=NULL;}
	delete p;
	p=NULL;
}


int get_detector_parameters(char *detector, double *LAT,double *LON, double *location, double *response_tensor)
{
	std::string local_det(detector);
	if(local_det.find("Hanford")!=std::string::npos ||
		local_det.find("hanford")!=std::string::npos){
		*LAT = H_LAT;
		*LON = H_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = H_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=Hanford_D[i][j];	
			}
		}
	}
	else if(local_det.find("Livingston")!=std::string::npos ||
		local_det.find("livingston")!=std::string::npos){
		*LAT = L_LAT;
		*LON = L_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = L_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=Livingston_D[i][j];	
			}
		}
	}
	else if(local_det.find("Virgo")!=std::string::npos ||
		local_det.find("virgo")!=std::string::npos){
		*LAT = V_LAT;
		*LON = V_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = V_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=Virgo_D[i][j];	
			}
		}
	}
	else if(local_det.find("Kagra")!=std::string::npos ||
		local_det.find("kagra")!=std::string::npos){
		*LAT = K_LAT;
		*LON = K_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = K_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=Kagra_D[i][j];	
			}
		}
	}
	else if(local_det.find("Indigo")!=std::string::npos ||
		local_det.find("indigo")!=std::string::npos){
		*LAT = I_LAT;
		*LON = I_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = I_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=Indigo_D[i][j];	
			}
		}
	}
	else if(local_det.find("CosmicExplorer")!=std::string::npos ||
		local_det.find("cosmicexplorer")!=std::string::npos ||
		local_det.find("CE")!=std::string::npos ){
		*LAT = CE_LAT;
		*LON = CE_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = CE_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=CE_D[i][j];	
			}
		}
	}
	else if(local_det.find("Einstein Telescope 1")!=std::string::npos ||
		local_det.find("einstein telescope 1")!=std::string::npos ||
		local_det.find("ET1")!=std::string::npos ){
		*LAT = ET1_LAT;
		*LON = ET1_LONG;
		for(int i= 0 ; i<3 ; i++){
			location[i] = ET1_location[i];
			for(int j= 0 ; j<3 ; j++){
				response_tensor[3*i+j]=ET1_D[i][j];	
			}
		}
	}
	else{
		std::cout<<"Unsupported detector"<<std::endl;
		return -1;
	}
	return 0;
}

int DL_from_Z_py(double z, char * COSMOLOGY, double *out)
{
	*out = DL_from_Z(z,std::string(COSMOLOGY));
	return 0;
}
int calculate_chirpmass_py(double mass1, double mass2,double *out)
{
	*out = calculate_chirpmass(mass1,mass2);
	return 0;
}

int calculate_eta_py(double mass1, double mass2,double *out)
{
	*out = calculate_eta(mass1,mass2);
	return 0;
}
int calculate_mass1_py(double chirpmass, double eta,double *out)
{
	*out = calculate_mass1(chirpmass,eta);
	return 0;
}
int calculate_mass2_py(double chirpmass, double eta,double *out)
{
	*out = calculate_mass2(chirpmass,eta);
	return 0;
}
void populate_noise_py(double *frequencies, char * detector, double *noise_root, int length, double integration_time){
	populate_noise(frequencies, std::string(detector), noise_root, length, integration_time);
	return;
}

void ppE_theory_fisher_transformation_py(double m1,
	double m2,
	double *spin1,
	double *spin2,
	double chip,
	double phip,
	double Luminosity_Distance,
	double phiRef,
	double tc,
	double RA,
	double DEC,
	double phi_l,
	double theta_l,
	double psi,
	double incl_angle,
	double gmst,
	bool reduced_spin,
	bool sky_average ,
	bool NSflag1 ,
	bool NSflag2 ,
	bool horizon_coord ,
	bool equatorial_orientation ,
	int Nmod,
	double *betappe,
	double **original_fisher,
	double **new_fisher,
	char * original_method,
	char * new_method,
	int dimension
	)
{
	gen_params params;
	params.mass1 = m1;
	params.mass2 = m2;
	if(!reduced_spin){
		params.spin1[0] = spin1[0];
		params.spin1[1] = spin1[1];
		params.spin1[2] = spin1[2];
		params.spin2[0] = spin2[0];
		params.spin2[1] = spin2[1];
		params.spin2[2] = spin2[2];
	}
	else{
		params.chip=chip;
		params.chip=phip;
		params.spin1[2]=spin1[2];
		params.spin2[2]=spin2[2];
	}
	params.Luminosity_Distance = Luminosity_Distance;
	params.phiRef = phiRef;
	params.tc = tc;
	params.sky_average = sky_average;
	params.NSflag1 = NSflag1;
	params.NSflag2 = NSflag2;
	params.horizon_coord = horizon_coord;
	params.equatorial_orientation = equatorial_orientation;
	params.RA = RA;
	params.DEC = DEC;
	params.phi_l = phi_l;
	params.theta_l = theta_l;
	params.psi = psi;
	params.incl_angle = incl_angle;
	params.gmst = gmst;
	params.Nmod = Nmod;
	params.betappe = new double[Nmod];
		
	for(int i = 0 ; i<Nmod; i++){
		params.betappe[i]=betappe[i];
	}
	
	ppE_theory_covariance_transformation(std::string(original_method),std::string(new_method),dimension, &params, original_fisher, new_fisher);
	delete [] params.betappe;
}
