#ifndef PPE_UTILITIES_H
#define PPE_UTILITIES_H
#include <functional>
#include "util.h"

//typedef T (*beta_fn)(source_parameters<T> *);
template <class T>
using beta_fn = T (*)(source_parameters<T> *);

template<class T>
struct theory_ppE_map{
	double *bppe=NULL;
	int Nmod;
	//T (*beta_fn)(source_parameters<T> *);
	beta_fn<T> *beta_fns;
	std::string ppE_method;
};

bool check_mod(std::string generation_method);
bool check_theory_support(std::string generation_method);

template<class T>
void assign_mapping(std::string generation_method,theory_ppE_map<T> *mapping,gen_params_base<T> *params_in);
template<class T>
void deallocate_mapping(theory_ppE_map<T> *mapping);

template<class T>
T dCS_beta(source_parameters<T> *param);
template<class T>
T dCS_phase_factor(source_parameters<T> *param);

template<class T>
T EdGB_beta( source_parameters<T> *param);

template<class T>
T EdGB_phase_factor( source_parameters<T> *param);

template<class T>
T ExtraDimension_beta( source_parameters<T> *param);

template<class T>
T TVG_beta( source_parameters<T> *param);

template<class T>
T DipRad_beta( source_parameters<T> *param);

template<class T>
T NonComm_beta( source_parameters<T> *param);

template<class T>
T ModDispersion_beta( source_parameters<T> *param);
int dispersion_lookup(double alpha);
template <class T>
T cosmology_interpolation_function_MD(T x, double *coeffs, int interp_degree);
template <class T>
T DL_from_Z_MD(T Z, double alpha);
#endif
