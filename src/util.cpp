#include "util.h"
//#include "general_parameter_structures.h"
#include "noise_util.h"
#include <math.h>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>
#include "adolc/adouble.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

const std::string PROJECT_DIRECTORY="/Users/sperkins/opt/gw_analysis_tools";

/*! \file
 *
 * General utilities that are method independent
 */
//#######################################################################################
//Interpolate Z to DL once per import
const int npts =100000;
double DLvec[npts];
double Zvec[npts];
gsl_interp_accel *Z_DL_accel_ptr ;
gsl_spline *Z_DL_spline_ptr ;
void initiate_LumD_Z_interp()
{
	Z_DL_accel_ptr = gsl_interp_accel_alloc();
	Z_DL_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,npts);
	std::fstream data_table;
	data_table.open(PROJECT_DIRECTORY+"/data/tabulated_LumD_Z.csv",std::ios::in);
	std::vector<std::string> row;
	std::string line, word, temp;
	int i =0,j=0;
	if(data_table){
		while(std::getline(data_table,line)){
			std::stringstream lineStream(line);	
			std::string item;
			while(std::getline(lineStream, item, ','))
			{
				if(i<npts){
				if(j==0){DLvec[i]=std::stod(item);}
				else if(j==1){Zvec[i]=std::stod(item);}
				j++;
				}
			}
			j = 0;
			i ++;
		}
	}
	gsl_spline_init(Z_DL_spline_ptr, DLvec, Zvec, npts);
	data_table.close();
}
void free_LumD_Z_interp()
{
	gsl_interp_accel_free(Z_DL_accel_ptr);
	gsl_spline_free(Z_DL_spline_ptr);
}
//#######################################################################################




adouble Z_from_DL(adouble DL)
{
	//std::fstream data_table;
	//data_table.open(PROJECT_DIRECTORY+"/data/tabulated_LumD_Z.csv",std::ios::in);
	//int npts =100000;
	//double DLvec[npts];
	//double Zvec[npts];
	//std::vector<std::string> row;
	//std::string line, word, temp;
	//int i =0,j=0;
	//if(data_table){
	//	while(std::getline(data_table,line)){
	//		std::stringstream lineStream(line);	
	//		std::string item;
	//		while(std::getline(lineStream, item, ','))
	//		{
	//			if(j==0){DLvec[i]=std::stod(item);}
	//			else if(j==1){Zvec[i]=std::stod(item);}
	//			j++;
	//		}
	//		j = 0;
	//		i ++;
	//	}
	//}
	//gsl_interp_accel *Z_DL_accel_ptr = gsl_interp_accel_alloc();
	//gsl_spline *Z_DL_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,npts);
	//gsl_spline_init(Z_DL_spline_ptr, DLvec, Zvec, npts);
	adouble Z = 0;
	Z = (adouble)gsl_spline_eval(Z_DL_spline_ptr, DL.value(), Z_DL_accel_ptr);
	return Z;
}

double Z_from_DL(double DL)
{
	//std::fstream data_table;
	//data_table.open(PROJECT_DIRECTORY+"/data/tabulated_LumD_Z.csv",std::ios::in);
	//int npts =100000;
	//double DLvec[npts];
	//double Zvec[npts];
	//std::vector<std::string> row;
	//std::string line, word, temp;
	//int i =0,j=0;
	//if(data_table){
	//	while(std::getline(data_table,line)){
	//		std::stringstream lineStream(line);	
	//		std::string item;
	//		while(std::getline(lineStream, item, ','))
	//		{
	//			if(i<npts){
	//			if(j==0){DLvec[i]=std::stod(item);}
	//			else if(j==1){Zvec[i]=std::stod(item);}
	//			j++;
	//			}
	//		}
	//		j = 0;
	//		i ++;
	//	}
	//}
	//gsl_interp_accel *Z_DL_accel_ptr = gsl_interp_accel_alloc();
	//gsl_spline *Z_DL_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,npts);
	//gsl_spline_init(Z_DL_spline_ptr, DLvec, Zvec, npts);
	double Z = 0;
	double DLtemp = DL;
	if(DL>DLvec[npts-1]){
		std::cout<<"WARNING: DL exceeded limit: setting to highest value in table"<<std::endl;
		DLtemp=DLvec[npts-1];
		std::cout<<DL<<std::endl;
		std::cout<<DLtemp<<std::endl;
	}
	
	Z = gsl_spline_eval(Z_DL_spline_ptr, DLtemp, Z_DL_accel_ptr);
	return Z;

}

void printProgress (double percentage)
{
    	int val = (int) (percentage * 100);
    	int lpad = (int) (percentage * PBWIDTH);
    	int rpad = PBWIDTH - lpad;
    	printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    	fflush (stdout);
}

/*! \brief Builds the structure that shuttles source parameters between functions -updated version to incorporate structure argument
 *
 * Populates the structure that is passed to all generation methods - contains all relavent source parameters 
 */
template <class T>
source_parameters<T> source_parameters<T>::populate_source_parameters(
			gen_params *param_in
			) 
{

	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters params;
	params.mass1 = param_in->mass1*MSOL_SEC;
	params.mass2 = param_in->mass2*MSOL_SEC;
	params.spin1x = param_in->spin1[0];
	params.spin2x = param_in->spin2[0];
	params.spin1y = param_in->spin1[1];
	params.spin2y = param_in->spin2[1];
	params.spin1z = param_in->spin1[2];
	params.spin2z = param_in->spin2[2];
	params.chi_s = (1./2)*(params.spin1z+params.spin2z);
	params.chi_a = (1./2)*(params.spin1z-params.spin2z);
	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
	params.chirpmass = calculate_chirpmass(params.mass1,params.mass2);
	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
	params.eta = calculate_eta(params.mass1,params.mass2);	
	params.M = params.mass1 + params.mass2;
	params.chi_eff = (params.mass1*(params.spin1z)+ params.mass2*(params.spin2z))/(params.M);
	params.chi_pn = params.chi_eff - (38*params.eta/113)*(2*params.chi_s);
	params.DL = param_in->Luminosity_Distance*MPC_SEC;
	params.delta_mass = sqrt(1.-4*params.eta);
	params.phic = param_in->phic;
	params.tc = param_in->tc;
	if (param_in->sky_average){
		params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
		params.sky_average=true;
	}
	else{
		params.A0 = sqrt(M_PI*40./192.)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
		params.sky_average=false;
	}
	return params;
}
/*! \brief Builds the structure that shuttles source parameters between functions- outdated in favor of structure argument 
 *
 * Populates the structure that is passed to all generation methods - contains all relavent source parameters 
 */
template <class T>
source_parameters<T> source_parameters<T>::populate_source_parameters_old(
			T mass1, /**< mass of the larger body - in Solar Masses*/ 
			T mass2, /**< mass of the smaller body - in Solar Masses*/
			T Luminosity_Distance,/**< Luminosity Distance in Mpc*/ 
			T *spin1,/** spin vector of the larger body  {sx,sy,sz}*/
			T *spin2, /** spin vector of the smaller body  {sx,sy,sz}*/
			T phi_c,/** coalescence phase*/
			T t_c ,/** coalescence time*/
			bool sky_average
			) 
{

	/* Convert all dimensionful quantities to seconds and build all needed source quantities once*/
	source_parameters params;
	params.mass1 = mass1*MSOL_SEC;
	params.mass2 = mass2*MSOL_SEC;
	params.spin1x = spin1[0];
	params.spin2x = spin2[0];
	params.spin1y = spin1[1];
	params.spin2y = spin2[1];
	params.spin1z = spin1[2];
	params.spin2z = spin2[2];
	params.chi_s = (1./2)*(params.spin1z+params.spin2z);
	params.chi_a = (1./2)*(params.spin1z-params.spin2z);
	//params.chirpmass = (adouble)calculate_chirpmass((double)params.mass1.value(),(double)params.mass2.value());
	params.chirpmass = calculate_chirpmass(params.mass1,params.mass2);
	//params.eta = (adouble)calculate_eta((double)params.mass1.value(),(double)params.mass2.value());	
	params.eta = calculate_eta(params.mass1,params.mass2);	
	params.M = params.mass1 + params.mass2;
	params.chi_eff = (params.mass1*(params.spin1z)+ params.mass2*(params.spin2z))/(params.M);
	params.chi_pn = params.chi_eff - (38*params.eta/113)*(2*params.chi_s);
	params.DL = Luminosity_Distance*MPC_SEC;
	params.delta_mass = sqrt(1.-4*params.eta);
	params.phic = phi_c;
	params.tc = t_c;
	//params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	if(sky_average){
		params.A0 = sqrt(M_PI/30)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	}
	else{
		params.A0 = sqrt(M_PI*40./192.)*params.chirpmass*params.chirpmass/params.DL * pow(M_PI*params.chirpmass,-7./6);
	}
	return params;
}


/*! \brief Calculates the chirp mass from the two component masses
 *
 * The output units are whatever units the input masses are
 */
double calculate_chirpmass(double mass1, double mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}
adouble calculate_chirpmass(adouble mass1, adouble mass2)
{
	return pow(mass1 * mass2,3./5) / pow(mass1 + mass2,1./5);
}

/*!\brief Calculates the symmetric mass ration from the two component masses
 */
double calculate_eta(double mass1, double mass2)
{
	return (mass1 * mass2) / pow(mass1 + mass2 ,2);
}
adouble calculate_eta(adouble mass1, adouble mass2)
{
	return (mass1 * mass2) / pow(mass1 + mass2 ,2);
}

/*! \brief Calculates the larger mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass1(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}

adouble calculate_mass1(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow + sqrt(1.-4*eta)*chirpmass / etapow);
}
/*! \brief Calculates the smaller mass given a chirp mass and symmetric mass ratio
 *
 * Units of the output match the units of the input chirp mass
 */
double calculate_mass2(double chirpmass, double eta)
{
	double etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}
adouble calculate_mass2(adouble chirpmass, adouble eta)
{
	adouble etapow = pow(eta,3./5);
    	return 1./2*(chirpmass / etapow - sqrt(1.-4*eta)*chirpmass / etapow);
}


long factorial(long num)
{
	int prod = 1;
	int step = num;
	while (step>0)
	{
		prod *= step;
		step -=1;
	}
	return prod;
}

double pow_int(double base, int power)
{
	if (power == 0) return 1.;
	double prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++){
		prod = prod * base;
	}
	if (power>0)
		return prod;
	else
		return 1./prod;
}
adouble pow_int(adouble base, int power)
{
	adouble prod = 1;
	int pow = std::abs(power);
	for (int i = 0; i< pow;i++)
		prod = prod * base;
	if (pow>0)
		return prod;
	else
		return 1./prod;
}

double cbrt_internal(double base)
{
	return cbrt(base);
}
adouble cbrt_internal(adouble base)
{
	return pow(base,1./3.);
}

/*! \brief Utility to malloc 2D array
 * 
 */
double** allocate_2D_array( int dim1, int dim2)
{
	double **array = (double **) malloc(sizeof(double*)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (double*)malloc(sizeof(double ) * dim2);
	}
	return array;
}

/*! \brief Utility to free malloc'd 2D array
 * 
 */
void deallocate_2D_array(double **array, int dim1, int dim2)
{
	for(int i =0; i < dim1; i++)
	{
		free(array[i]);
	}
	free(array);
}
/*! \brief Utility to malloc 3D array
 * 
 */
double*** allocate_3D_array( int dim1, int dim2, int dim3)
{
	double ***array = (double ***) malloc(sizeof(double**)*dim1);
	for (int i = 0; i<dim1; i ++)
	{
		array[i] = (double**)malloc(sizeof(double *) * dim2);
		for (int j = 0 ; j< dim2; j++)
		{
			array[i][j] = (double *)malloc(sizeof(double) * dim3);
		}
	}
	return array;
}
/*! \brief Utility to free malloc'd 2D array
 * 
 */
void deallocate_3D_array(double ***array, int dim1, int dim2, int dim3)
{
	for(int i =0; i < dim1; i++)
	{
		for(int j =0 ; j<dim2; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}
	free(array);
}

/*!\brief Utility to read in data
 *
 * Takes filename, and assigns to output[rows][cols]
 *
 * File must be comma separated doubles
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
		double **output, /**<[out] array to store output, dimensions rowsXcols*/
		int rows, /**< first dimension*/
		int cols /**<second dimension*/
		)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word;
	int i=0, j=0;
	double *temp = (double *)malloc(sizeof(double)*rows*cols);
	
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				temp[i]=std::stod(item);	
				i+=1;	
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
	for(i =0; i<rows;i++){
		for(j=0; j<cols;j++)
			output[i][j] = temp[cols*i + j];
	}
	free(temp);
}


/*!\brief Utility to read in data (single dimension vector) 
 *
 * Takes filename, and assigns to output[i*rows + cols]
 *
 * Output vector must be long enough, no check is done for the length
 *
 * File must be comma separated doubles
 */
void read_file(std::string filename, /**< input filename, relative to execution directory*/
	double *output /**<[out] output array, assumed to have the proper length of total items*/
	)
{
	std::fstream file_in;
	file_in.open(filename, std::ios::in);
	std::string line, word, temp;
	int i =0;
	if(file_in){
		while(std::getline(file_in, line)){
			std::stringstream lineStream(line);
			std::string item;
			while(std::getline(lineStream,item, ',')){
				output[i]=std::stod(item);	
				i+=1;
			}	
		}	
	}
	else{std::cout<<"ERROR -- File "<<filename<<" not found"<<std::endl;exit(1);}
}

/*! \brief Utility to write 2D array to file
 *
 * Grid of data, comma separated
 *
 * Grid has rows rows and cols columns
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		double **input, /**< Input 2D array pointer array[rows][cols]*/
		int rows, /**< First dimension of array*/
		int cols /**< second dimension of array*/
		)
{
	
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	for(int i =0; i<rows; i++){
		for(int j=0; j<cols;j++){
			if(j==cols-1)
				out_file<<input[i][j]<<std::endl;
			else
				out_file<<input[i][j]<<" , ";
		}
	}
	out_file.close();
}
/*! \brief Utility to write 1D array to file
 *
 * Single column of data
 */
void write_file(std::string filename, /**<Filename of output file, relative to execution directory*/
		double *input, /**< input 1D array pointer array[length]*/
		int length /**< length of array*/
		)
{
	std::ofstream out_file;
	out_file.open(filename);
	out_file.precision(15);
	for(int j =0; j<length; j++)
		out_file<<input[j]<<std::endl;
	out_file.close();

}



template <class T>
std::complex<T> cpolar(T mag, T phase)
{
	return mag * std::exp(std::complex<T>(0,1) * phase);
}

/*! Shamelessly stolen from LALsuite
 *
 */
template <class T>
std::complex<T> XLALSpinWeightedSphericalHarmonic(
                                   T theta,  /**< polar angle (rad) */
                                   T phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    )
{
  T fac = 0.0;
  std::complex<T> ans = std::complex<T>(0.0,0.0);

  /* sanity checks ... */
  //if ( l < abs(s) ) 
  //{
  //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", __func__, s, l, m );
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  //if ( l < abs(m) ) 
  //{
  //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  if ( s == -2 ) 
  {
    if ( l == 2 ) 
    {
      switch ( m ) 
      {
        case -2:
          fac = sqrt( 5.0 / ( 64.0 * M_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
          break;
        case -1:
          fac = sqrt( 5.0 / ( 16.0 * M_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
          break;

        case 0:
          fac = sqrt( 15.0 / ( 32.0 * M_PI ) ) * sin( theta )*sin( theta );
          break;

        case 1:
          fac = sqrt( 5.0 / ( 16.0 * M_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
          break;

        case 2:
          fac = sqrt( 5.0 / ( 64.0 * M_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      } /*  switch (m) */
    }  /* l==2*/
    else if ( l == 3 ) 
    {
      switch ( m ) 
      {
        case -3:
          fac = sqrt(21.0/(2.0*M_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
          break;
        case -2:
          fac = sqrt(7.0/(4.0*M_PI))*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
          break;
        case -1:
          fac = sqrt(35.0/(2.0*M_PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;
        case 0:
          fac = (sqrt(105.0/(2.0*M_PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
          break;
        case 1:
          fac = -sqrt(35.0/(2.0*M_PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;

        case 2:
          fac = sqrt(7.0/M_PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
          break;

        case 3:
          fac = -sqrt(21.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
          break;

        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }   /* l==3 */
    else if ( l == 4 ) 
    {
      switch ( m ) 
      {
        case -4:
          fac = 3.0*sqrt(7.0/M_PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = 3.0*sqrt(7.0/(2.0*M_PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
          break;

        case -2:
          fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(M_PI));
          break;
        case -1:
          fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*M_PI));
          break;
        case 0:
          fac = (3.0*sqrt(5.0/(2.0*M_PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
          break;
        case 1:
          fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*M_PI));
          break;
        case 2:
          fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(M_PI));
          break;
        case 3:
          fac = -3.0*sqrt(7.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
          break;
        case 4:
          fac = 3.0*sqrt(7.0/M_PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }    /* l==4 */
    else if ( l == 5 ) 
    {
      switch ( m ) 
      {
        case -5:
          fac = sqrt(330.0/M_PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
          break;
        case -4:
          fac = sqrt(33.0/M_PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = (sqrt(33.0/(2.0*M_PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
          break;
        case -2:
          fac = (sqrt(11.0/M_PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
          break;
        case -1:
          fac = (sqrt(77.0/M_PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
          break;
        case 0:
          fac = (sqrt(1155.0/(2.0*M_PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
          break;
        case 1:
          fac = sqrt(77.0/M_PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
          break;
        case 2:
          fac = sqrt(11.0/M_PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
          break;
        case 3:
          fac = -sqrt(33.0/(2.0*M_PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
          break;
        case 4:
          fac = sqrt(33.0/M_PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
          break;
        case 5:
          fac = -sqrt(330.0/M_PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    }  /* l==5 */
    else if ( l == 6 )
    {
      switch ( m )
      {
        case -6:
          fac = (3.*sqrt(715./M_PI)*pow(cos(theta/2.0),4)*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),3)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),7))/2.0;
          break;
        case -4:
          fac = (sqrt(195./(2.0*M_PI))*pow(cos(theta/2.0),2)*(35. + 44.*cos(theta) 
          + 33.*cos(2.*theta))*pow(sin(theta/2.0),6))/8.0;
          break;
        case -3:
          fac = (3.*sqrt(13./M_PI)*cos(theta/2.0)*(98. + 185.*cos(theta) + 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*pow(sin(theta/2.0),5))/32.0;
          break;
        case -2:
          fac = (sqrt(13./M_PI)*(1709. + 3096.*cos(theta) + 2340.*cos(2.*theta) + 1320.*cos(3.*theta) 
          + 495.*cos(4.*theta))*pow(sin(theta/2.0),4))/256.0;
          break;
        case -1:
          fac = (sqrt(65./(2.0*M_PI))*cos(theta/2.0)*(161. + 252.*cos(theta) + 252.*cos(2.*theta) 
          + 132.*cos(3.*theta) + 99.*cos(4.*theta))*pow(sin(theta/2.0),3))/64.0;
          break;
        case 0:
          fac = (sqrt(1365./M_PI)*(35. + 60.*cos(2.*theta) + 33.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (sqrt(65./(2.0*M_PI))*pow(cos(theta/2.0),3)*(161. - 252.*cos(theta) + 252.*cos(2.*theta) 
          - 132.*cos(3.*theta) + 99.*cos(4.*theta))*sin(theta/2.0))/64.0;
          break;
        case 2:
          fac = (sqrt(13./M_PI)*pow(cos(theta/2.0),4)*(1709. - 3096.*cos(theta) + 2340.*cos(2.*theta) 
          - 1320*cos(3*theta) + 495*cos(4*theta)))/256.0;
          break;
        case 3:
          fac = (-3.*sqrt(13./M_PI)*pow(cos(theta/2.0),5)*(-98. + 185.*cos(theta) - 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*sin(theta/2.0))/32.0;
          break;
        case 4:
          fac = (sqrt(195./(2.0*M_PI))*pow(cos(theta/2.0),6)*(35. - 44.*cos(theta) 
          + 33.*cos(2*theta))*pow(sin(theta/2.0),2))/8.0;
          break;
        case 5:
          fac = -(sqrt(2145./M_PI)*pow(cos(theta/2.0),7)*(-1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
          break;
        case 6:
          fac = (3.*sqrt(715./M_PI)*pow(cos(theta/2.0),8)*pow(sin(theta/2.0),4))/2.0;
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==6 */
    else if ( l == 7 )
    {
      switch ( m )
      {
        case -7:
          fac = sqrt(15015./(2.0*M_PI))*pow(cos(theta/2.0),5)*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),4)*(2. + 7.*cos(theta))*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),3)*(93. + 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),2)*(140. + 285.*cos(theta) 
          + 156.*cos(2.*theta) + 91.*cos(3.*theta))*pow(sin(theta/2.0),6))/16.0;
          break;
        case -3:
          fac = (sqrt(15./(2.0*M_PI))*cos(theta/2.0)*(3115. + 5456.*cos(theta) + 4268.*cos(2.*theta) 
          + 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(15./M_PI)*(5220. + 9810.*cos(theta) + 7920.*cos(2.*theta) + 5445.*cos(3.*theta) 
          + 2860.*cos(4.*theta) + 1001.*cos(5.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (3.*sqrt(5./(2.0*M_PI))*cos(theta/2.0)*(1890. + 4130.*cos(theta) + 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) + 1430.*cos(4.*theta) + 1001.*cos(5*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(35./M_PI)*cos(theta)*(109. + 132.*cos(2.*theta) 
          + 143.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (3.*sqrt(5./(2.0*M_PI))*pow(cos(theta/2.0),3)*(-1890. + 4130.*cos(theta) - 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) - 1430.*cos(4.*theta) + 1001.*cos(5.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(15./M_PI)*pow(cos(theta/2.0),4)*(-5220. + 9810.*cos(theta) - 7920.*cos(2.*theta) 
          + 5445.*cos(3.*theta) - 2860.*cos(4.*theta) + 1001.*cos(5.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(15./(2.0*M_PI))*pow(cos(theta/2.0),5)*(3115. - 5456.*cos(theta) + 4268.*cos(2.*theta) 
          - 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*sin(theta/2.0))/128.0;
          break;  
        case 4:
          fac = (sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),6)*(-140. + 285.*cos(theta) - 156.*cos(2*theta) 
          + 91.*cos(3.*theta))*pow(sin(theta/2.0),2))/16.0;
          break;
        case 5:
          fac = -(sqrt(165./(2.0*M_PI))*pow(cos(theta/2.0),7)*(93. - 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = (sqrt(2145./M_PI)*pow(cos(theta/2.0),8)*(-2. + 7.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
          break;
        case 7:
          fac = -(sqrt(15015./(2.0*M_PI))*pow(cos(theta/2.0),9)*pow(sin(theta/2.0),5));
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==7 */
    else if ( l == 8 )
    {
      switch ( m )
      {
        case -8:
          fac = sqrt(34034./M_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),10);
          break;
        case -7:
          fac = sqrt(17017./(2.0*M_PI))*pow(cos(theta/2.0),5)*(1. + 4.*cos(theta))*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = sqrt(255255./M_PI)*pow(cos(theta/2.0),4)*(1. + 2.*cos(theta))
          *sin(M_PI/4.0 - theta/2.0)*sin(M_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),8);
          break;
        case -5:
          fac = (sqrt(12155./(2.0*M_PI))*pow(cos(theta/2.0),3)*(19. + 42.*cos(theta) 
          + 21.*cos(2.*theta) + 14.*cos(3.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(935./(2.0*M_PI))*pow(cos(theta/2.0),2)*(265. + 442.*cos(theta) + 364.*cos(2.*theta) 
          + 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),6))/32.0;
          break;
        case -3:
          fac = (sqrt(561./(2.0*M_PI))*cos(theta/2.0)*(869. + 1660.*cos(theta) + 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) + 455.*cos(4.*theta) + 182.*cos(5.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(17./M_PI)*(7626. + 14454.*cos(theta) + 12375.*cos(2.*theta) + 9295.*cos(3.*theta) 
          + 6006.*cos(4.*theta) + 3003.*cos(5.*theta) + 1001.*cos(6.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (sqrt(595./(2.0*M_PI))*cos(theta/2.0)*(798. + 1386.*cos(theta) + 1386.*cos(2.*theta) 
          + 1001.*cos(3.*theta) + 858.*cos(4.*theta) + 429.*cos(5.*theta) + 286.*cos(6.*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(595./M_PI)*(210. + 385.*cos(2.*theta) + 286.*cos(4.*theta) 
          + 143.*cos(6.*theta))*pow(sin(theta),2))/4096.0;
          break;
        case 1:
          fac = (sqrt(595./(2.0*M_PI))*pow(cos(theta/2.0),3)*(798. - 1386.*cos(theta) + 1386.*cos(2.*theta) 
          - 1001.*cos(3.*theta) + 858.*cos(4.*theta) - 429.*cos(5.*theta) + 286.*cos(6.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(17./M_PI)*pow(cos(theta/2.0),4)*(7626. - 14454.*cos(theta) + 12375.*cos(2.*theta) 
          - 9295.*cos(3.*theta) + 6006.*cos(4.*theta) - 3003.*cos(5.*theta) + 1001.*cos(6.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(561./(2.0*M_PI))*pow(cos(theta/2.0),5)*(-869. + 1660.*cos(theta) - 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) - 455.*cos(4.*theta) + 182.*cos(5.*theta))*sin(theta/2.0))/128.0;
          break;
        case 4:
          fac = (sqrt(935./(2.0*M_PI))*pow(cos(theta/2.0),6)*(265. - 442.*cos(theta) + 364.*cos(2.*theta) 
          - 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),2))/32.0;
          break;
        case 5:
          fac = -(sqrt(12155./(2.0*M_PI))*pow(cos(theta/2.0),7)*(-19. + 42.*cos(theta) - 21.*cos(2.*theta) 
          + 14.*cos(3.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = sqrt(255255./M_PI)*pow(cos(theta/2.0),8)*(-1. + 2.*cos(theta))*sin(M_PI/4.0 - theta/2.0)
          *sin(M_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),4);
          break;
        case 7:
          fac = -(sqrt(17017./(2.0*M_PI))*pow(cos(theta/2.0),9)*(-1. + 4.*cos(theta))*pow(sin(theta/2.0),5));
          break;
        case 8:
          fac = sqrt(34034./M_PI)*pow(cos(theta/2.0),10)*pow(sin(theta/2.0),6);
          break;
        //default:
        //  XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
        //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
        //  break;
      }
    } /* l==8 */
    //else 
    //{
    //  XLALPrintError("XLAL Error - %s: Unsupported mode l=%d (only l in [2,8] implemented)\n", __func__, l);
    //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
    //}
  }
  //else 
  //{
  //  XLALPrintError("XLAL Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", __func__, s);
  //  XLAL_ERROR_VAL(0, XLAL_EINVAL);
  //}
  if (m)
    ans = cpolar((T)(1.0), (T)(m*phi)) * fac;
  else
    ans = fac;
  return ans;
}

template std::complex<double> XLALSpinWeightedSphericalHarmonic<double>(double,double,int,int,int);
template std::complex<adouble> XLALSpinWeightedSphericalHarmonic<adouble>(adouble,adouble,int,int,int);
template std::complex<double> cpolar<double>(double, double);
template std::complex<adouble> cpolar<adouble>(adouble,adouble);
template class source_parameters<double>;
template class source_parameters<adouble>;
