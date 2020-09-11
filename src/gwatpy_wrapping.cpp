#include "gwatpy_wrapping.h"
#include "util.h"
#include "fisher.h"
#include "waveform_generator.h"


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


gen_params_base<double>* gen_params_base_py(double mass1, double mass2)
{
	gen_params_base<double> *p = new gen_params_base<double> ;
	p->mass1 = mass1;	
	p->mass2 = mass2;	
	p->Luminosity_Distance = 100;	
	p->spin1[0] = 0;	
	p->spin1[1] = 0;	
	p->spin1[2] = 0;	
	p->spin2[0] = 0;	
	p->spin2[1] = 0;	
	p->spin2[2] = 0;	
	p->incl_angle = M_PI;

	return p;
}

void gen_params_base_py_destructor(gen_params_base<double> *p)
{
	delete p;
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
