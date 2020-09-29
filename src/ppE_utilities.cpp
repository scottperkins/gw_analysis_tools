#include "ppE_utilities.h"
#include "util.h"
#include "D_Z_Config_modified_dispersion.h"
/*! \brief Utilities for ppE parameterizations of gravitational waves
 *
 * Supported theories:
 *
 * dCS -- input beta is \alpha^2 in seconds^4
 * 
 * EdGB -- input beta is \alpha^2 in seconds^4
 * 
 * EdGB_GHOv1 -- EdGB Generic Higher Order version 1-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 *
 * EdGB_GHOv2 -- EdGB Generic Higher Order version 2-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 * 
 * EdGB_GHOv3 -- EdGB Generic Higher Order version 3-- input beta is \alpha^2 in seconds^4 and \gamma, which is dimensionless
 *
 * ExtraDimension -- input beta l^2 in seconds^2
 *
 * TVG -- Time varying G -- input beta is \dot{G}_{z} in Hz
 *
 * DipRad -- Generic Dipole Radiation -- input beta is \delta \dot{E} (dimensionless)
 *
 * NonComm -- Noncommutatitve gravity -- input beta is \Lambda^2 in units of Planck energies (dimensionless), where \sqrt(\Lambda) defines the energy scale of noncommutativity 
 *
 * ModDispersion -- Modified Dispersion -- input beta is A_alpha (eV^(2-alpha))
 */



bool check_mod(std::string generation_method)
{
	if(generation_method.find("ppE") != std::string::npos)
	{
		return true;
	}
	if(generation_method.find("gIMR") != std::string::npos)
	{
		return true;
	}
	if(check_theory_support(generation_method)){
		return true;
	}
	return false;
}
bool check_theory_support(std::string generation_method)
{
	if(generation_method.find("dCS")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB")!=std::string::npos){
		return true;
	}
	if(generation_method.find("EdGB_GHO")!=std::string::npos){
		return true;
	}
	if(generation_method.find("ExtraDimension")!=std::string::npos){
		return true;
	}
	if(generation_method.find("TVG")!=std::string::npos){
		return true;
	}
	if(generation_method.find("DipRad")!=std::string::npos){
		return true;
	}
	if(generation_method.find("NonComm")!=std::string::npos){
		return true;
	}
	if(generation_method.find("ModDispersion")!=std::string::npos){
		return true;
	}
	return false;
} 
template<class T>
void deallocate_mapping(theory_ppE_map<T> *mapping){
	if(mapping->bppe){
		delete [] mapping->bppe;
	}
	if(mapping->beta_fns){
		delete [] mapping->beta_fns;
	}
}
template void deallocate_mapping(theory_ppE_map<adouble> *mapping);
template void deallocate_mapping(theory_ppE_map<double> *mapping);

template<class T>
void assign_mapping(std::string generation_method,theory_ppE_map<T> *mapping, gen_params_base<T> *params_in)
{
	bool ins=false;
	if(generation_method.find("dCS")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		//mapping->beta_fn = new T (**)(source_parameters<T> *)[mapping->Nmod]; 
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->bppe[0] = -1;
		mapping->beta_fns[0] = &dCS_beta ;
		ins = true;
	}
	else if(generation_method.find("EdGB")!= std::string::npos){
		if(generation_method.find("EdGB_GHOv1")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = &EdGB_beta ;
			mapping->beta_fns[1] = &EdGB_GHO_betav1 ;
		}
		if(generation_method.find("EdGB_GHOv2")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = &EdGB_beta ;
			mapping->beta_fns[1] = &EdGB_GHO_betav2 ;
		}
		if(generation_method.find("EdGB_GHOv3")!= std::string::npos){
			mapping->Nmod = 2;
			mapping->bppe = new double[2];
			mapping->bppe[0] =-7;
			mapping->bppe[1] =-5;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = &EdGB_beta ;
			mapping->beta_fns[1] = &EdGB_GHO_betav3 ;
		}
		else{
			mapping->Nmod = 1;
			mapping->bppe = new double[1];
			mapping->bppe[0] =-7;
			mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
			mapping->beta_fns[0] = &EdGB_beta ;
		}
		ins = true;
	}
	else if(generation_method.find("ExtraDimension")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &ExtraDimension_beta ;
		ins = true;
	}
	else if(generation_method.find("TVG")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &TVG_beta ;
		ins = true;
	}
	else if(generation_method.find("DipRad")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-7;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &DipRad_beta ;
		ins = true;
			
	}
	else if(generation_method.find("NonComm")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		mapping->bppe[0] =-1;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &NonComm_beta ;
		ins = true;
			
	}
	else if(generation_method.find("ModDispersion")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new double[1];
		//mapping->bppe[0] =-1;
		mapping->bppe[0] =params_in->bppe[0];
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &ModDispersion_beta ;
		ins = false;
			
	}
	if(generation_method.find("Pv2")!=std::string::npos){
		if(ins){
			mapping->ppE_method = "ppE_IMRPhenomPv2_Inspiral";
		}
		else{
			mapping->ppE_method = "ppE_IMRPhenomPv2_IMR";
		}
	}
	else if(generation_method.find("PhenomD")!=std::string::npos){
		if(ins){
			mapping->ppE_method = "ppE_IMRPhenomD_Inspiral";
		}
		else{
			mapping->ppE_method = "ppE_IMRPhenomD_IMR";
		}

	}
	return ;
}
template void assign_mapping(std::string,theory_ppE_map<double>*,gen_params_base<double>*);
template void assign_mapping(std::string,theory_ppE_map<adouble>*,gen_params_base<adouble>*);

template<class T>
T dCS_beta(source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = (param->betappe[0]);
	T out =  16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) *dCS_phase_factor(param);
	return out;
}
template double dCS_beta(source_parameters<double> *);
template adouble dCS_beta(source_parameters<adouble> *);

template<class T>
T dCS_phase_factor(source_parameters<T> *param)
{
	T g=0;
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T coeff1 = -5./8192.;
	T coeff2 = 15075./114688.;
	T m1 = calculate_mass1(chirpmass,eta);
	T m2 = calculate_mass2(chirpmass,eta);
	T m = m1+m2;
	T chi1 = param->chi_s+param->chi_a;
	T chi2 = param->chi_s-param->chi_a;
	T s1temp = 2.+2.*pow_int(chi1,4) - 2.*sqrt((1.-chi1*chi1)) - chi1*chi1 * ((3. - 2.*sqrt(1.-chi1*chi1)));
	T s2temp = 2.+2.*pow_int(chi2,4) - 2.*sqrt((1.-chi2*chi2)) - chi2*chi2 * ((3. - 2.*sqrt(1.-chi2*chi2)));
	
	T s1  ;
	T s2  ;
	if( fabs(chi1) <DOUBLE_COMP_THRESH){
		s1 = 0;		
	}
	else{
		s1  = s1temp/(2.*chi1*chi1*chi1);
	}
	if( fabs(chi2) <DOUBLE_COMP_THRESH){
		s2 = 0;		
	}
	else{
		s2  = s2temp/(2.*chi2*chi2*chi2);
	}
	//Neutron stars don't source scalar charge
	if(param->NSflag1){s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	g+=coeff1/(pow(eta,14./5.)) * pow((m1*s2 - m2 * s1),2.)/(m*m);
	g+=coeff2/(pow(eta,14./5.)) * (m2*m2* chi1*chi1 - 350./201. * m1*m2*chi1*chi2 + m1*m1 * chi2*chi2)/(m*m);
	return g;
}
template double dCS_phase_factor(source_parameters<double> *);
template adouble dCS_phase_factor(source_parameters<adouble> *);

template<class T>
T EdGB_beta( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * EdGB_phase_factor(param);
} 
template adouble EdGB_beta(source_parameters<adouble> *);
template double EdGB_beta(source_parameters<double> *);


template<class T>
T EdGB_phase_factor( source_parameters<T> *param)
{
 	T M = param->M;	
 	T chirpmass = param->chirpmass;	
 	T eta = param->eta;	
	T m1 = calculate_mass1(chirpmass, eta);
	T m2 = calculate_mass2(chirpmass, eta);
	T chi1 = param->chi_s + param->chi_a;
	T chi2 = param->chi_s - param->chi_a;
	T temp1 = 2.*(sqrt(1.-chi1*chi1) - 1. + chi1*chi1);
	T temp2 = 2.*(sqrt(1.-chi2*chi2) - 1. + chi2*chi2);
	T s1;
        T s2;
	if( fabs(chi1) <DOUBLE_COMP_THRESH){
		s1 = 0;		
	}
	else{
		s1  = temp1/(chi1*chi1);
	}
	if( fabs(chi2) <DOUBLE_COMP_THRESH){
		s2 = 0;		
	}
	else{
		s2  = temp2/(chi2*chi2);
	}
	if(param->NSflag1){debugger_print(__FILE__,__LINE__,"NS 1");s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	return (-5./7168.)* pow_int((m1*m1 * s2 - m2*m2 * s1),2) / (pow_int(M,4) * pow(eta,(18./5)));
} 
template adouble EdGB_phase_factor(source_parameters<adouble> *);
template double EdGB_phase_factor(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav1( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * generic_mod;
} 
template adouble EdGB_GHO_betav1(source_parameters<adouble> *);
template double EdGB_GHO_betav1(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav2( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	T phase_factor = EdGB_phase_factor(param);
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) *phase_factor* generic_mod;
} 
template adouble EdGB_GHO_betav2(source_parameters<adouble> *);
template double EdGB_GHO_betav2(source_parameters<double> *);

template<class T>
T EdGB_GHO_betav3( source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T unredshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	T generic_mod = param->betappe[1];
	T phase_factor = EdGB_phase_factor(param);
	T phase_GR_1PN = (3715./756. + 55. * param->eta / 9.);
	return 16.*M_PI*phase_mod/(pow_int(unredshiftedM,4)) * generic_mod* phase_factor*phase_GR_1PN;
} 
template adouble EdGB_GHO_betav3(source_parameters<adouble> *);
template double EdGB_GHO_betav3(source_parameters<double> *);

template<class T>
T ExtraDimension_beta( source_parameters<T> *param)
{
	T ED_length_sq = param->betappe[0]; //Length in seconds
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T ten_mircometer = 10.e-6 / c; //10 micrometers in seconds
	T m1dot = -2.8e-7 * pow_int(MSOL_SEC * (1+Z)/ param->mass1,2)* ED_length_sq /pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
	T m2dot = -2.8e-7 * pow_int(MSOL_SEC* (1+Z)/ param->mass2,2)* ED_length_sq/pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
	T beta = (m1dot + m2dot) * (25. / 851968.) * 
		( ( 3. - 26.*param->eta + 34. * param->eta*param->eta) / 
		( pow(param->eta, 2./5.) * ( 1- 2*param->eta)) );

	return beta;
} 
template adouble ExtraDimension_beta(source_parameters<adouble> *);
template double ExtraDimension_beta(source_parameters<double> *);

template<class T>
T TVG_beta( source_parameters<T> *param)
{
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T Gdot = param->betappe[0];
	T beta = ( -25. / 65526. ) * ( Gdot * param->chirpmass/(1+Z)) ;
	return beta;
} 
template adouble TVG_beta(source_parameters<adouble> *);
template double TVG_beta(source_parameters<double> *);

template<class T>
T DipRad_beta( source_parameters<T> *param)
{
	T deltaEdot =param->betappe[0];
	T beta = ( -3. / 224. ) * pow(param->eta,2./5.) * deltaEdot;
	return beta;
} 
template adouble DipRad_beta(source_parameters<adouble> *);
template double DipRad_beta(source_parameters<double> *);

template<class T>
T NonComm_beta( source_parameters<T> *param)
{
	T Lambda_sq = param->betappe[0];
	T beta = ( -75. / 256. ) * pow(param->eta,-4./5.) * ( 2.*param->eta - 1.) * Lambda_sq;
	return beta;
} 
template adouble NonComm_beta(source_parameters<adouble> *);
template double NonComm_beta(source_parameters<double> *);

//#############################################################
/* \brief Modified dispersion beta
 *
 * https://arxiv.org/abs/1110.2720
 *
 * b = (alpha*3 - 3)
 *
 * Supported alpha and corresponding b and PN order:
 *
 * 	0 	-3	1
 *
 * 	0.5 	-1.5	1.75
 *
 * 	1.0 	0 	2.5
 *
 * 	1.5 	1.5 	3.25
 *
 * 	2.0 	3	4
 *
 * 	2.5 	4.5	4.75
 *
 * 	3.0 	6	5.5
 *
 * 	3.5 	7.5	6.25
 *
 * 	4.0 	9	7
 */
template<class T>
T ModDispersion_beta( source_parameters<T> *param)
{
	double alpha = (param->bppe[0]+3.)/3.;
	T Z= Z_from_DL(param->DL/MPC_SEC,param->cosmology);
	T Dalpha = DL_from_Z_MD(Z,alpha)*MPC_SEC;
	//T lambdaA = h_planck* pow(param->betappe[0],1./(alpha-2));	
	//T beta = ( pow(M_PI,2.-alpha) / (1.-alpha) ) 
	//	* ( Dalpha / pow(param->betappe[0],2.-alpha) ) 
	//	* ( pow(param->chirpmass,1.-alpha) / pow(1+Z, 1.-alpha) ) ;
	T beta = ( pow(M_PI,2.-alpha) / (1.-alpha) ) 
		* ( Dalpha * param->betappe[0]/pow(h_planck, 2.-alpha) ) 
		* ( pow(param->chirpmass,1.-alpha) / pow(1+Z, 1.-alpha) ) ;
	return beta;
} 
template adouble ModDispersion_beta(source_parameters<adouble> *);
template double ModDispersion_beta(source_parameters<double> *);

/*! \brief Calculates the modified distance given the redshift for modified dispersion relationships
 *
 * Based on Astropy.cosmology calculations -- see python script in the ./data folder of the project -- numerically calculated given astropy.cosmology's definitions (http://docs.astropy.org/en/stable/cosmology/) and used scipy.optimize to fit to a power series, stepping in half powers of Z. These coefficients are then output to a header file (D_Z_config_modified_dispersion.h) which are used here to calculate distance. Custom cosmologies etc can easily be acheived by editing the python script D_Z_config.py, the c++ functions do not need modification. They use whatever data is available in the header file. If the functional form of the fitting function changes, these functions DO need to change.
 *
 * Only Planck15 is implemented.
 *
 * Supported alphas include 0,.5,1,1.5,2,2.5,3,3.5,4
 *
 * See https://arxiv.org/abs/1110.2720 for full definitions
 * 
 */
template <class T>
T DL_from_Z_MD(T Z, double alpha)
{
	int dispersion_index = dispersion_lookup(alpha);
	if (dispersion_index == -1){ std::cout<<"Invalid Dispersion power"<<std::endl;return -1;}
	const double *boundaries = MD_boundaries_Z[dispersion_index];
	int interp_deg = interp_MD_degree[dispersion_index];
	int num_seg = num_MD_segments[dispersion_index];
	T dl;
	for (int i =0; i<num_seg; i++){
		if ( Z<boundaries[i+1]){
			//if(i == 0 ){
			//	debugger_print(__FILE__,__LINE__,"Boundaries: "+std::to_string(boundaries[i])+" "+std::to_string(boundaries[i+1]) +" "+std::to_string(Z));
			//}
			double *coeffs = new double [interp_deg];
			for (int j =0; j<interp_deg;j++)
				coeffs[j]=MD_COEFF_VEC_ZD[dispersion_index][i][j];
			dl =  cosmology_interpolation_function_MD(Z,coeffs, interp_deg);
			delete[] coeffs;
			return dl;
			
		}	
	}
	return -1;
}
template adouble DL_from_Z_MD(adouble, double);
template double DL_from_Z_MD(double, double);

template <class T>
T cosmology_interpolation_function_MD(T x, double *coeffs, int interp_degree)
{
	T result=0;
	for(int i = 0 ; i<interp_degree; i++){
		result+= coeffs[i] * pow(x,-3.5 + i*0.5);
	}
	return result;

}
template adouble cosmology_interpolation_function_MD(adouble, double*,int);
template double cosmology_interpolation_function_MD(double, double*,int);
/*! \brief Helper function for mapping cosmology name to an internal index
 */
int dispersion_lookup(double alpha)
{
	for (int i =0; i<num_MD_alphas; i++){
		if (fabs(alpha- MD_alphas[i]) < DOUBLE_COMP_THRESH){
			return i;
		}
	}
	return -1;
}
//#############################################################
