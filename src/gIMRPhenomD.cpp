#include "gIMRPhenomD.h"
#include "util.h"
#include "IMRPhenomD.h"
#include <adolc/adouble.h>

/*! \file 
 *
 * gIMR parameterization for modifications to GR
 *
 * four groups -- phi, beta, sigma, alpha
 *
 * beta can be 2,3
 *
 * alpha can be 2,3,4
 *
 * sigma can be 2,3,4
 *
 * phi can be (-4)-9 (8, 9 correspond to the logarithmic terms at 5 and 6)
 *
 * ***Issues with logarithmic terms at the moment***
 */

/*! \brief Utility to map between the \beta parameter in ppE parameterization to the \delta \phi parameter in gIMR for STATIC PN coefficients only
 */
template<class T>
void gIMRPhenomD<T>::ppE_gIMR_mapping(gen_params_base<T> *parameters, int PN_order,T *beta,int *b  )
{
	source_parameters<T> params;
	//params = params.populate_source_parameters(parameters);
	params.populate_source_parameters(parameters);
	params.phi = parameters->phi;
	params.theta = parameters->theta;
	params.incl_angle = parameters->incl_angle;
	params.f_ref = parameters->f_ref;
	params.phiRef = parameters->phiRef;
	params.cosmology = parameters->cosmology;
	params.shift_time = parameters->shift_time;
	params.sky_average = parameters->sky_average;
	params.shift_phase = parameters->shift_phase;
	params.NSflag1 = parameters->NSflag1;
	params.NSflag2 = parameters->NSflag2;
	params.dep_postmerger = parameters->dep_postmerger;
	
	int id = check_list_id(PN_order, parameters->phii, parameters->Nmod_phi);
	if(id != -1){
		T static_PN_coeffs[11];
		IMRPhenomD<T> model_base;
		model_base.assign_static_pn_phase_coeff(&params, static_PN_coeffs);

		*beta = (3./ 128.) * pow(params.eta, (T)(-PN_order)/5.) * static_PN_coeffs[PN_order] * parameters->delta_phi[id];
		*b = PN_order - 5;
	}
	else{
		std::cout<<"ERROR -- gIMR modification does not exist for PN order: "+std::to_string(PN_order)<<std::endl;
	}
}

//Handles modifications for alpha, beta, and sigma
template<class T>
void gIMRPhenomD<T>::assign_lambda_param(source_parameters<T> *source_param, 
	lambda_parameters<T> *lambda)
{
	IMRPhenomD<T> modelgr;
	modelgr.assign_lambda_param(source_param, lambda);
	//Modify the parameters
	if(source_param->Nmod_beta != 0){
		for(int i = 2 ; i<4; i++){
			int id = check_list_id(i, source_param->betai, source_param->Nmod_beta);
			if(id != -1){
				lambda->beta[i]*=(1. + source_param->delta_beta[id]);	
			}
		}
	}
	if(source_param->Nmod_alpha != 0){
		for(int i = 2 ; i<5; i++){
			int id = check_list_id(i, source_param->alphai, source_param->Nmod_alpha);
			if(id != -1){
				lambda->alpha[i]*=(1. + source_param->delta_alpha[id]);	
			}
		}
	}
	if(source_param->Nmod_sigma != 0){
		for(int i = 2 ; i<5; i++){
			int id = check_list_id(i, source_param->sigmai, source_param->Nmod_sigma);
			if(id != -1){
				lambda->sigma[i]*=(1. + source_param->delta_sigma[id]);	
			}
		}
	}
}

//Covers modifications 0-7
template <class T>
void gIMRPhenomD<T>::assign_static_pn_phase_coeff(source_parameters<T> *source_param, T *coeff)
{
	IMRPhenomD<T> modelgr;
	modelgr.assign_static_pn_phase_coeff(source_param, coeff);
	
	if(source_param->Nmod_phi != 0){
		for(int i = 0 ;i < 8; i++){
			if(i != 6 && i!=5 && i != 1){
				int id = check_list_id(i, source_param->phii, 
					source_param->Nmod_phi);
				if(id != -1){
					coeff[i]*=(1. + source_param->delta_phi[id]);	
				}
			}
			else{
				if(i == 6){
					int id = check_list_id(i, source_param->phii, 
						source_param->Nmod_phi);
					if(id != -1){
						coeff[10]*=(1. + source_param->delta_phi[id]);	
					}
				}
				if(i == 5){
					int id = check_list_id(i, source_param->phii, 
						source_param->Nmod_phi);
					if(id != -1){
						coeff[8]*=(1. + source_param->delta_phi[id]);	
					}
				}
				if(i == 1){
					int id = check_list_id(i, source_param->phii, 
						source_param->Nmod_phi);
					if(id != -1){
						coeff[1]=( source_param->delta_phi[id]);	
					}
				}
					
			}
		}

		int id = check_list_id(8, source_param->phii, 
			source_param->Nmod_phi);
		if(id != -1){
			coeff[9]*=(1. + source_param->delta_phi[id]);	
		}
		id = check_list_id(9, source_param->phii, 
			source_param->Nmod_phi);
		if(id != -1){
			coeff[11]*=(1. + source_param->delta_phi[id]);	
		}
	}
}


template<class T>
T gIMRPhenomD<T>::phase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda, useful_powers<T> *pows)
{
	IMRPhenomD<T> modelgr;
	T phase = modelgr.phase_ins(f, param, pn_coeff, lambda, pows);
	if(param->Nmod_phi != 0 ){
		T pimfcube = pow(M_PI*param->M * f, (T)1./3.);
		for(int i = -4 ;i < 0; i++){
			int id = check_list_id(i, param->phii, 
				param->Nmod_phi);
			if(id != -1){
				phase+=3./(128.*param->eta) *(param->delta_phi[id])*pow_int(pimfcube, i-5);	
			}
		}
	}
	return phase;

}
template<class T>
T gIMRPhenomD<T>::Dphase_ins(T f, source_parameters<T> *param, T *pn_coeff, lambda_parameters<T> *lambda)
{
	IMRPhenomD<T> modelgr;
	T Dphase = modelgr.Dphase_ins(f, param, pn_coeff, lambda);
	if(param->Nmod_phi != 0 ){
		T pimcube = pow(M_PI*param->M , (T)1./3.);
		for(int i = -4 ;i < 0; i++){
			int id = check_list_id(i, param->phii, 
				param->Nmod_phi);
			if(id != -1){
				Dphase+=3./(128.*param->eta) *(param->delta_phi[id])*
					pow_int(pimcube, i-5) * (i-5.)/3. * 
					pow(f, ((i-5.)/3. - 1));	
			}
		}
	}
	return Dphase;

}


template class gIMRPhenomD<double>;
template class gIMRPhenomD<adouble>;
