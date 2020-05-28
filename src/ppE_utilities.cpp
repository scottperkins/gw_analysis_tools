#include "ppE_utilities.h"
#include "util.h"
/*! \brief Utilities for ppE parameterizations of gravitational waves
 *
 * Supported theories:
 *
 * dCS -- input beta is \alpha^2 in seconds^4
 * 
 * EdGB -- input beta is \alpha^2 in seconds^4
 *
 * ExtraDimension -- input beta l^2 in seconds^2
 *
 * TVG -- Time varying G -- input beta is \dot{G}_{z} in Hz
 *
 * DipRad -- Generic Dipole Radiation -- input beta is \delta \dot{E} (dimensionless)
 *
 * NonComm -- Noncommutatitve gravity -- input beta is \Lambda^2 in units of Planck energies (dimensionless), where \sqrt(\Lambda) defines the energy scale of noncommutativity 
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
void assign_mapping(std::string generation_method,theory_ppE_map<T> *mapping)
{
	bool ins=false;
	if(generation_method.find("dCS")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		//mapping->beta_fn = new T (**)(source_parameters<T> *)[mapping->Nmod]; 
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->bppe[0] = -1;
		mapping->beta_fns[0] = &dCS_beta ;
		ins = true;
	}
	else if(generation_method.find("EdGB")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		mapping->bppe[0] =-7;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &EdGB_beta ;
		ins = true;
	}
	else if(generation_method.find("ExtraDimension")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &ExtraDimension_beta ;
		ins = true;
	}
	else if(generation_method.find("TVG")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		mapping->bppe[0] =-13;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &TVG_beta ;
		ins = true;
	}
	else if(generation_method.find("DipRad")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		mapping->bppe[0] =-7;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &DipRad_beta ;
		ins = true;
			
	}
	else if(generation_method.find("NonComm")!= std::string::npos){
		mapping->Nmod = 1;
		mapping->bppe = new int[1];
		mapping->bppe[0] =-1;
		mapping->beta_fns = new beta_fn<T>[mapping->Nmod]; 
		mapping->beta_fns[0] = &NonComm_beta ;
		ins = true;
			
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
template void assign_mapping(std::string,theory_ppE_map<double>*);
template void assign_mapping(std::string,theory_ppE_map<adouble>*);

template<class T>
T dCS_beta(source_parameters<T> *param)
{
 	T M = param->M;	
	T DL = param->DL;
	T Z= Z_from_DL(DL/MPC_SEC,param->cosmology);
	T redshiftedM = M/(1.+Z);
	T phase_mod = (param->betappe[0]);
	T out =  16.*M_PI*phase_mod/(pow_int(redshiftedM,4)) *dCS_phase_factor(param);
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
	T redshiftedM = M/(1.+Z);
	T phase_mod = param->betappe[0];
	return 16.*M_PI*phase_mod/(pow_int(redshiftedM,4)) * EdGB_phase_factor(param);
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
	if(param->NSflag1){s1 =0;}	
	if(param->NSflag2){s2 =0;}	
	return (-5./7168.)* pow_int((m1*m1 * s2 - m2*m2 * s1),2) / (pow_int(M,4) * pow(eta,(18./5)));
} 
template adouble EdGB_phase_factor(source_parameters<adouble> *);
template double EdGB_phase_factor(source_parameters<double> *);

template<class T>
T ExtraDimension_beta( source_parameters<T> *param)
{
	T ED_length_sq = param->betappe[0]; //Length in seconds
	T ten_mircometer = 10.e-6 / c; //10 micrometers in seconds
	T m1dot = -2.8e-7 * pow_int(MSOL_SEC/ param->mass1,2)* ED_length_sq /pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
	T m2dot = -2.8e-7 * pow_int(MSOL_SEC/ param->mass2,2)* ED_length_sq/pow_int( ten_mircometer, 2) *MSOL_SEC/T_year;
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
	T Gdot = param->betappe[0];
	T beta = ( -25. / 65526. ) * ( Gdot * param->chirpmass ) ;
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
