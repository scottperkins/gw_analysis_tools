#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include <iostream>



int AD_v_N(int argc, char *argv[]);
void RT_ERROR_MSG();

int main(int argc, char *argv[])
{
	std::cout<<"TESTING FISHER CALCULATIONS"<<std::endl;
		
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = stoi(argv[1]);	
	if(runtime_opt == 0){
		return AD_v_N(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int AD_v_N(int argc, char *argv[])
{
	
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .1;
	params.spin2[2] = -.1;
	params.chip = .3;
	params.phip = 1.0;
	params.Luminosity_Distance = 400;
	params.phiRef = .0;
	params.RA = 2.;
	params.DEC = -0.9;
	params.f_ref = 20;

	params.NSflag1 = false;
	params.NSflag2 = false;
	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.mass1 = 36;
	params.mass2 = 29;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	//params.theta_l = 1.1;
	//params.phi_l = 2.1;
	params.psi = 1.;
	params.gmst = 2.;

	//double fmin = 3e-2;
	//double fmax = 1e-1;
	//double T = T_year/2;

	//double fmin = 5;
	//double fmax = 2048;
	double fmin = f_0PN(4*T_year,chirpmass);
	double fmax = 1;
	//double T = 32;
	double T = 4*T_year;

	//params.tc = 3.*T/4.;
	params.tc = T;
	//params.tc = 0;

	//int length = T*(fmax-fmin);
	int length = 1000;
	double *frequency = new double[length];
	//int Ndetect = 4;
	int Ndetect = 1;
	double **psd = new double*[Ndetect];
	//std::string SN[4] = {"CE2_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[4] = {"LISA_SADC_CONF","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	//std::string SN[3] = {"Hanford_O1_fitted","Hanford_O1_fitted","Hanford_O1_fitted"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	//for(int i = 0 ; i<length; i++){
	//	frequency[i]=fmin + (double)i /T;
	//}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	//std::string detector = "LISA";
	std::string detector = "Hanford";
	//std::string method = "IMRPhenomD";
	std::string method = "ppE_IMRPhenomPv2_Inspiral";
	//transform_orientation_coords(&params, method, detector);

	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new int[1];
	params.betappe[0] = 0;
	//params.bppe[0] = -1.;
	params.bppe[0] = -13;
	//params.bppe[0] = 1;

	//std::string detectors[4] = {"CE","Hanford","Livingston","Virgo"};
	std::string detectors[4] = {"LISA","Hanford","Livingston","Virgo"};
		
	int dim = 13;
	int dimD = 12;
	int dimDSA = 8;

	double **jac_spins = allocate_2D_array(dim,dim);
	for (int i = 0 ;i<dim; i++){
		for(int j =0 ;j<dim; j++){
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

	double **output_N = allocate_2D_array(dim,dim);
	double **output_N_temp = allocate_2D_array(dim,dim);
	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD3 = allocate_2D_array(dimD,dimD);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD3_temp = allocate_2D_array(dimD,dimD);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD3 = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD2 = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_N[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimD; i++){
		for(int j = 0 ; j<dimD; j++){
			output_AD3[i][j]= 0;
			output_AD3_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
		}
	}

	params.equatorial_orientation = true;
	//params.equatorial_orientation = false;
	params.sky_average = true;
	params.theta_l = 0.12;
	params.phi_l = 0.1;
	params.incl_angle  = 0;
	double snr = calculate_snr(SN[0],"LISA",method, &params, frequency, length, "GAUSSLEG",weights,true);
	std::cout<<snr<<std::endl;
	double SNR_TARGET = 100;
	params.Luminosity_Distance = snr/SNR_TARGET*params.Luminosity_Distance;
	params.sky_average = false;
	params.incl_angle = M_PI-.01;



	//for(int i = 0 ;i < Ndetect; i++){
	//	fisher_numerical(frequency, length, method, detectors[i], output_N_temp, dim, &params, 2, NULL,NULL, psd[i]);
	//	for(int k = 0 ; k<dim; k++){
	//		for(int j = 0 ; j<dim; j++){
	//			output_N[k][j]+= output_N_temp[k][j];
	//		}
	//	}
	//}
	
	//std::cout<<"NUMERICAL:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_N[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	snr = calculate_snr(SN[0],"CE",method, &params, frequency, length, "GAUSLEG",weights,true);
	std::cout<<"SNR: "<<snr<<std::endl;
	//snr = calculate_snr(SN[1],"Hanford",method, &params, frequency, length, "SIMPSONS",NULL,false);
	//std::cout<<"SNR (Hanford): "<<snr<<std::endl;
	
	
	//std::cout<<"AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<COV_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	method = "ppE_IMRPhenomD_Inspiral";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i], output_AD3_temp, dimD, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimD; k++){
			for(int j = 0 ; j<dimD; j++){
				output_AD3[k][j]+= output_AD3_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD3, jac_spins,output_AD3_temp,dimD,dimD,dimD);
	matrix_multiply(jac_spins,output_AD3_temp, output_AD3,dimD,dimD,dimD);
	
	
	//std::cout<<"AD-D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<output_AD3[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_AD3,COV_AD3,dimD);
	std::cout<<"COV AD - D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<COV_AD3[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dimD; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD3[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	method = "ppE_IMRPhenomD_Inspiral";
	params.sky_average = true;
	params.incl_angle = 0;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	
	//std::cout<<"AD-DSA:"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<output_ADSA[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_ADSA,COV_ADSA,dimDSA);
	std::cout<<"COV AD - DSA:"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<COV_ADSA[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_ADSA[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"SNR (SA): "<<sqrt(output_ADSA[0][0])<<std::endl;
	snr = calculate_snr(SN[0],"CE",method, &params, frequency, length, "GAUSSLEG",weights,true);
	std::cout<<"SNR (SA -- CE): "<<snr<<std::endl;
	//snr = calculate_snr(SN[1],"Hanford",method, &params, frequency, length, "SIMPSONS",NULL,false);
	//std::cout<<"SNR (SA -- Hanford): "<<snr<<std::endl;



	//##########################################################################
	//TESTING
	//If you use 13 parameters (phip and phiRef), these two parameters are exactly,
	//linearly correlated. Simply use one. 
	//##########################################################################
	//params.phiRef-=1;
	//params.phip+=1;
	//fisher_autodiff(frequency, length, method, detector, output_AD2, dim, &params, "SIMPSONS",NULL,false, psd,NULL,NULL);
	//
	//std::cout<<"FD AD / AD2:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_AD2[i][j] - output_AD[i][j])*2./(output_AD2[i][j] + output_AD[i][j])<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//##########################################################################
	//##########################################################################
	//std::cout<<"FRACTIONAL DIFF (N-AD)*2/(N+AD):"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_N[i][j] - output_AD[i][j])*2./(fabs(output_N[i][j]) + fabs(output_AD[i][j]))<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	//std::cout<<"DIFF (N-AD):"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<(output_N[i][j] - output_AD[i][j])<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD3,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dim,dim);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD3,dimD,dimD);
	deallocate_2D_array(output_AD3_temp,dimD,dimD);
	deallocate_2D_array(output_N,dim,dim);
	deallocate_2D_array(output_N_temp,dim,dim);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] psd;
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Compare AD to Numerical"<<std::endl;
}


