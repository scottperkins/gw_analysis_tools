#include "gwat/util.h"
#include "gwat/fisher.h"
#include "gwat/detector_util.h"
#include "gwat/waveform_util.h"
#include "gwat/ortho_basis.h"
#include "gwat/pn_waveform_util.h"
#include <iostream>



int AD_v_N(int argc, char *argv[]);
int network_fishers(int argc, char *argv[]);
int dCS_EdGB(int argc, char *argv[]);
int test_jac_transform(int argc, char *argv[]);
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
	if(runtime_opt == 1){
		return network_fishers(argc,argv);
	}
	if(runtime_opt == 2){
		return dCS_EdGB(argc,argv);
	}
	if(runtime_opt == 3){
		return test_jac_transform(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int test_jac_transform(int argc, char *argv[])
{
	std::cout.precision(15);
	gen_params params;	
	params.mass1 = 14.9;
	params.mass2 = 8.3;
	params.spin1[2] = .18* (params.mass1+params.mass2)/params.mass1;
	params.spin2[2] = 0.1 ;
	params.chip = .01;
	params.phip = 1.0;
	params.Luminosity_Distance = 50;
	params.incl_angle = M_PI/3.;

	params.NSflag1 = false;
	params.NSflag2 =false;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	std::cout<<"chirpmass: "<<chirpmass/MSOL_SEC<<std::endl;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;
	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new int[1];
	params.betappe[0] =0;
	//params.betappe[0] =1;
	params.bppe[0] = -13;
	

	double fmin = 5;
	double fmax = 2048;
	double T = 32;

	params.tc = 3.*T/4.;
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 3;
	//int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//std::string SN[3] = {"CE1","AdLIGOMidHigh","AdVIRGOPlus1"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	std::string method = "ExtraDimension_IMRPhenomPv2";
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
	int dim = 13;
	int dimDSA = 8;

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD2 = allocate_2D_array(dim,dim);
	double **output_AD2_T = allocate_2D_array(dim,dim);
	double **output_AD2_temp = allocate_2D_array(dim,dim);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2 = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2_T = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA2_temp = allocate_2D_array(dimDSA,dimDSA);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD2 = allocate_2D_array(dim,dim);
	double **COV_AD2_T = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **COV_ADSA2 = allocate_2D_array(dimDSA,dimDSA);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_AD2[i][j]= 0;
			output_AD2_T[i][j]= 0;
			output_AD2_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
			output_ADSA2[i][j]= 0;
			output_ADSA2_temp[i][j]= 0;
		}
	}



	double snr; 

	std::cout<<"dCS fishers: "<<std::endl;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	//std::cout<<"AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<output_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

	std::cout<<"dCS fishers done "<<std::endl;
	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	//gsl_cholesky_matrix_invert(output_AD,COV_AD,dim);
	//std::cout<<"COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim; j++){
	//		std::cout<<COV_AD[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//std::cout<<"Variances (90%):"<<std::endl;
	//for(int i = 0 ; i<dim; i++){
	//	std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	//}
	//std::cout<<std::endl;
	//std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	//if(params.NSflag2){
	//	std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	//}
	//else{
	//	std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	//}
	//std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD[dim-1][dim-1],1./8.)*3.e5<<std::endl;
	//std::cout<<std::endl;


	double **identity_full = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD,COV_AD,identity_full,dim,dim,dim);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<identity_full[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity_full, dim,dim);




	//double **sub_AD_F = allocate_2D_array(dim-4,dim-4);
	//int ids[4] = {0,1,2,3};
	//rm_fisher_dim(output_AD,dim, sub_AD_F,dim-4,ids);
	//gsl_LU_matrix_invert(sub_AD_F,COV_AD,dim-4);
	//std::cout<<"SUB COV AD:"<<std::endl;
	//for(int i = 0 ; i<dim-4; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dim-4; j++){
	//		std::cout<<sub_AD_F[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//std::cout<<"Variances (90%):"<<std::endl;
	//for(int i = 0 ; i<dim-4; i++){
	//	std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	//}
	//std::cout<<std::endl;
	//std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	//if(params.NSflag2){
	//	std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	//}
	//else{
	//	std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	//}
	//std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
	//std::cout<<std::endl;
	//deallocate_2D_array(sub_AD_F,dim-4,dim-4);
	std::cout<<"ppE fishers: "<<std::endl;
	method = "ppE_IMRPhenomPv2_Inspiral";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD2_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD2[k][j]+= output_AD2_temp[k][j];
			}
		}
	}
	ppE_theory_fisher_transformation("ppE_IMRPhenomPv2_Inspiral","ExtraDimension_IMRPhenomPv2",dim, &params, output_AD2,output_AD2_T);
	//double **derivatives = new double*[1];
	//derivatives[0]=new double[13];
	//ppE_theory_transformation_calculate_derivatives("ppE_IMRPhenomPv2_Inspiral","dCS_IMRPhenomPv2",dim,12, &params, derivatives);
	//for(int i = 0; i<dim; i++){
	//	std::cout<<derivatives[0][i]<<std::endl;
	//}	
	
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<i<<" "<<j<<" "<<(output_AD[i][j]-output_AD2_T[i][j])*2./
				(fabs(output_AD[i][j])+fabs(output_AD2_T[i][j]))
				//<<std::endl;
				<<" "<<output_AD[i][j]<<" "<<output_AD2_T[i][j]<<" "<<output_AD2[i][j]<<std::endl;;
		}
	}
	gsl_LU_matrix_invert(output_AD2,COV_AD2,dim);
	ppE_theory_covariance_transformation("ppE_IMRPhenomPv2_Inspiral","ExtraDimension_IMRPhenomPv2",dim, &params, COV_AD2,COV_AD2_T);
	
	std::cout<<"Covariance difference: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<i<<" "<<j<<" "<<(COV_AD[i][j]-COV_AD2_T[i][j])*2./
				(fabs(COV_AD[i][j])+fabs(COV_AD2_T[i][j]))
				//<<std::endl;
				<<" "<<COV_AD[i][j]<<" "<<COV_AD2_T[i][j]<<" "<<COV_AD2[i][j]<<std::endl;;
		}
	}



	double **identity_full2 = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD2_T,COV_AD2_T,identity_full2,dim,dim,dim);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<identity_full2[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity_full2, dim,dim);





	//matrix_multiply(output_AD2, jac_spins,output_AD2_temp,dimD,dimD,dimD);
	//matrix_multiply(jac_spins,output_AD2_temp, output_AD2,dimD,dimD,dimD);
	
	
	std::cout<<"AD-2:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<output_AD2[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//gsl_LU_matrix_invert(output_AD2,COV_AD2,dimD);
	//std::cout<<"COV AD - D:"<<std::endl;
	//std::cout<<"Variances (90%):"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" "<<1.64*sqrt(COV_AD2[i][i])<<std::endl;
	//}
	//std::cout<<std::endl;
	//std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	//if(params.NSflag2){
	//	std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	//}
	//else{
	//	std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	//}
	//std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD2[dimD-1][dimD-1],1./8.)*3.e5<<std::endl;
	//std::cout<<std::endl;

	params.sky_average = true;
	params.incl_angle = 0;
	method = "ExtraDimension_IMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	method = "ppE_IMRPhenomD_Inspiral";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA2_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA2[k][j]+= output_ADSA2_temp[k][j];
			}
		}
	}
	ppE_theory_fisher_transformation("ppE_IMRPhenomD_Inspiral","ExtraDimension_IMRPhenomD",dimDSA, &params, output_ADSA2,output_ADSA2_T);
	
	//std::cout<<"AD-DSA:"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<output_ADSA[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//gsl_LU_matrix_invert(output_ADSA,COV_ADSA,dimDSA);
	//double **identity = allocate_2D_array(dimDSA,dimDSA);
	//matrix_multiply(output_ADSA,COV_ADSA,identity,dimDSA,dimDSA,dimDSA);
	//std::cout<<"IDENTITY: "<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	for(int j = 0 ; j<dimDSA; j++){
	//		std::cout<<identity[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//
	//deallocate_2D_array(identity, dimDSA,dimDSA);
	//std::cout<<"Variances (90%):"<<std::endl;
	//for(int i = 0 ; i<dimDSA; i++){
	//	std::cout<<i<<" "<<1.64*sqrt(COV_ADSA[i][i])<<std::endl;
	//}
	//std::cout<<std::endl;
	//std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	//if(params.NSflag2){
	//	std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	//}
	//else{
	//	std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	//}
	//std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_ADSA[dimDSA-1][dimDSA-1],1./8.)*3.e5<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<"SNR (SA): "<<sqrt(output_ADSA[0][0])<<std::endl;
	


	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			std::cout<<i<<" "<<j<<" "<<(output_ADSA[i][j]-output_ADSA2_T[i][j])*2./
				(fabs(output_ADSA[i][j])+fabs(output_ADSA2_T[i][j]))
				//<<std::endl;
				<<" "<<output_ADSA[i][j]<<" "<<output_ADSA2_T[i][j]<<" "<<output_ADSA2[i][j]<<std::endl;;
		}
	}

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD2,dim,dim);
	deallocate_2D_array(COV_AD2_T,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(COV_ADSA2,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2_T,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA2_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dim,dim);
	deallocate_2D_array(output_AD2_T,dim,dim);
	deallocate_2D_array(output_AD2_temp,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] params.betappe;
	delete [] params.bppe;
	//delete [] params.phii;
	//delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;

}
int network_fishers(int argc, char *argv[])
{
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .1;
	params.spin2[2] = -.1;
	params.chip = .3;
	params.phip = 1.0;
	params.Luminosity_Distance = 600;
	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = -0.0;
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
	params.psi = 1.;
	params.gmst = 2.;

	//double fmin = 3e-2;
	//double fmax = 1e-1;
	//double T = T_year/2;

	double fmin = 5;
	double fmax = 2048;
	//double fmin = f_0PN(4*T_year,chirpmass);
	//double fmax = 1;
	double T = 32;
	//double T = 4*T_year;

	params.tc = 3.*T/4.;
	//params.tc = T;
	//params.tc = 0;

	//int length = T*(fmax-fmin);
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 2;
	double **psd = new double*[Ndetect];
	std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	std::string method = "ppE_IMRPhenomPv2_Inspiral";
	//transform_orientation_coords(&params, method, detector);

	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new int[1];
	params.betappe[0] = 0;
	//params.bppe[0] = -1.;
	params.bppe[0] = -13;
	//params.bppe[0] = 1;
	params.Nmod_phi = 1;
	params.delta_phi = new double[1];
	params.phii = new int[1];
	params.delta_phi[0] = 0;
	//params.bppe[0] = -1.;
	params.phii[0] = 4;
	//params.bppe[0] = 1;

	//std::string detectors[4] = {"CE","Hanford","Livingston","Virgo"};
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
	int dim = 13;

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

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD3 = allocate_2D_array(dim,dim);
	double **output_AD3_temp = allocate_2D_array(dim,dim);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD3 = allocate_2D_array(dim,dim);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
			output_AD3[i][j]= 0;
			output_AD3_temp[i][j]= 0;
		}
	}

	params.equatorial_orientation = false;
	params.sky_average = false;
	params.incl_angle = M_PI-.01;
	params.incl_angle = M_PI-.01;


	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	//matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	
	
	std::cout<<"AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD3_temp, dim, &params,"GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD3[k][j]+= output_AD3_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD3, jac_spins,output_AD3_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD3_temp, output_AD3,dim,dim,dim);
	
	
	std::cout<<"AD-Network:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD3[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	gsl_LU_matrix_invert(output_AD3,COV_AD3,dim);
	std::cout<<"COV AD - Network:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD3[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD3[i][i])<<std::endl;
	}
	std::cout<<std::endl;


	std::cout<<"FRACTIONAL DIFF (Netowrk-AD)*2/(Network+AD):"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<(output_AD3[i][j] - output_AD[i][j])*2./(fabs(output_AD3[i][j]) + fabs(output_AD[i][j]))<<" ";
		}
		std::cout<<std::endl;
	}


	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD3,dim,dim);
	deallocate_2D_array(output_AD3,dim,dim);
	deallocate_2D_array(output_AD3_temp,dim,dim);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] psd;
	return 0;
}
int dCS_EdGB(int argc, char *argv[])
{
	
	std::cout.precision(15);
	gen_params params;	
	//params.mass1 = 31.5;
	//params.mass2 = 10.3;
	//params.spin1[2] = .38;
	//params.spin2[2] = (.21*(params.mass1+params.mass2)-params.mass1* params.spin1[2])/params.mass2 ;
	//params.chip = .29;
	//params.phip = 1.0;
	//params.Luminosity_Distance = 730;
	//params.incl_angle = .76;
	
	params.mass1 = 14.9;
	params.mass2 = 8.3;
	params.spin1[2] = .18* (params.mass1+params.mass2)/params.mass1;
	params.spin2[2] = 0 ;
	params.chip = .01;
	params.phip = 1.0;
	params.Luminosity_Distance = 50;
	params.incl_angle = M_PI/3.;

	//params.mass1 = 5.0;
	//params.mass2 = 1.4;
	//params.spin1[2] = .801;
	//params.spin2[2] = -.001;//(.21*(params.mass1+params.mass2)-params.mass1* params.spin1[2])/params.mass2 ;
	//params.chip = .001;
	//params.phip = 1.0;
	//params.Luminosity_Distance = 70;
	//params.incl_angle = M_PI/3.;
	

	params.NSflag1 = false;
	params.NSflag2 =false;

	params.phiRef = .0;
	params.RA = 1.;
	params.DEC = 0.6;
	params.f_ref = 20;
	double chirpmass = calculate_chirpmass(params.mass1,params.mass2)*MSOL_SEC;
	//params.spin1[2] = .38;

	params.horizon_coord = false;
	params.shift_time=false;
	params.shift_phase=false;
	params.dep_postmerger=true;
	params.sky_average=false;
	
	params.psi = 1.;
	params.gmst = 2.;
	params.sky_average = false;
	//params.Nmod_phi =1;
	//params.phii = new int[1];
	//params.phii[0]=4;
	//params.delta_phi=new double[1];
	//params.delta_phi[0]=0;
	params.Nmod = 1;
	params.betappe = new double[1];
	params.bppe = new int[1];
	//params.betappe[0] = pow_int(8*1000./(c),4);
	params.betappe[0] =0;
	//params.betappe[0] =1;
	params.bppe[0] = -7;
	

	double fmin = 5;
	double fmax = 2048;
	double T = 32;

	params.tc = 3.*T/4.;
	int length = 1000;
	double *frequency = new double[length];
	int Ndetect = 3;
	//int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGOMidHigh","AdLIGOMidHigh","AdVIRGOPlus1"};
	//std::string SN[3] = {"CE1","AdLIGOMidHigh","AdVIRGOPlus1"};
	
	double *weights = new double[length];
	gauleg(log10(fmin), log10(fmax),frequency,weights,length);
	for(int i = 0 ; i<length; i++){
		frequency[i] = pow(10,frequency[i]);	
	}
	for(int i = 0 ; i<Ndetect; i++){
		psd[i]= new double[length];
		populate_noise(frequency, SN[i],psd[i], length, 48);
		//populate_noise(frequency, "LISA_CONF",psd, length, 12);
		for(int j = 0 ; j<length; j++){
			psd[i][j]*=psd[i][j];	
		}
	}

	std::string method = "dCS_IMRPhenomPv2";
	//std::string method = "ppE_IMRPhenomPv2_IMR";
	//std::string method = "gIMRPhenomPv2";
	//std::string method = "EdGB_IMRPhenomPv2";


	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
		
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

	double **output_AD = allocate_2D_array(dim,dim);
	double **output_AD_temp = allocate_2D_array(dim,dim);
	double **output_AD2 = allocate_2D_array(dimD,dimD);
	double **output_ADSA = allocate_2D_array(dimDSA,dimDSA);
	double **output_ADSA_temp = allocate_2D_array(dimDSA,dimDSA);
	double **output_AD2_temp = allocate_2D_array(dimD,dimD);
	double **COV_AD = allocate_2D_array(dim,dim);
	double **COV_AD2 = allocate_2D_array(dim,dim);
	double **COV_ADSA = allocate_2D_array(dimDSA,dimDSA);
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			output_AD[i][j]= 0;
			output_AD_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimD; i++){
		for(int j = 0 ; j<dimD; j++){
			output_AD2[i][j]= 0;
			output_AD2_temp[i][j]= 0;
		}
	}
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			output_ADSA[i][j]= 0;
			output_ADSA_temp[i][j]= 0;
		}
	}



	double snr; 

	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	std::cout<<"AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<output_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	gsl_LU_matrix_invert(output_AD,COV_AD,dim);
	//gsl_cholesky_matrix_invert(output_AD,COV_AD,dim);
	std::cout<<"COV AD:"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim; j++){
			std::cout<<COV_AD[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dim; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD[dim-1][dim-1],1./8.)*3.e5<<std::endl;
	std::cout<<std::endl;


	double **identity_full = allocate_2D_array(dim,dim);
	matrix_multiply(output_AD,COV_AD,identity_full,dim,dim,dim);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dim; i++){
		for(int j = 0 ; j<dim; j++){
			std::cout<<identity_full[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity_full, dim,dim);




	double **sub_AD_F = allocate_2D_array(dim-4,dim-4);
	int ids[4] = {0,1,2,3};
	rm_fisher_dim(output_AD,dim, sub_AD_F,dim-4,ids);
	gsl_LU_matrix_invert(sub_AD_F,COV_AD,dim-4);
	std::cout<<"SUB COV AD:"<<std::endl;
	for(int i = 0 ; i<dim-4; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dim-4; j++){
			std::cout<<sub_AD_F[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dim-4; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD[dim-5][dim-5],1./8.)*3.e5<<std::endl;
	std::cout<<std::endl;
	deallocate_2D_array(sub_AD_F,dim-4,dim-4);






	method = "dCS_IMRPhenomD";
	//method = "ppE_IMRPhenomD_IMR";
	//method = "gIMRPhenomD";
	//method = "EdGB_IMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD2_temp, dimD, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimD; k++){
			for(int j = 0 ; j<dimD; j++){
				output_AD2[k][j]+= output_AD2_temp[k][j];
			}
		}
	}
	//matrix_multiply(output_AD2, jac_spins,output_AD2_temp,dimD,dimD,dimD);
	//matrix_multiply(jac_spins,output_AD2_temp, output_AD2,dimD,dimD,dimD);
	
	
	std::cout<<"AD-D:"<<std::endl;
	//for(int i = 0 ; i<dimD; i++){
	//	std::cout<<i<<" ";
	//	for(int j = 0 ; j<dimD; j++){
	//		std::cout<<output_AD2[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	gsl_LU_matrix_invert(output_AD2,COV_AD2,dimD);
	std::cout<<"COV AD - D:"<<std::endl;
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dimD; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_AD2[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_AD2[dimD-1][dimD-1],1./8.)*3.e5<<std::endl;
	std::cout<<std::endl;

	params.sky_average = true;
	params.incl_angle = 0;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	
	std::cout<<"AD-DSA:"<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		std::cout<<i<<" ";
		for(int j = 0 ; j<dimDSA; j++){
			std::cout<<output_ADSA[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	gsl_LU_matrix_invert(output_ADSA,COV_ADSA,dimDSA);
	double **identity = allocate_2D_array(dimDSA,dimDSA);
	matrix_multiply(output_ADSA,COV_ADSA,identity,dimDSA,dimDSA,dimDSA);
	std::cout<<"IDENTITY: "<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		for(int j = 0 ; j<dimDSA; j++){
			std::cout<<identity[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(identity, dimDSA,dimDSA);
	std::cout<<"Variances (90%):"<<std::endl;
	for(int i = 0 ; i<dimDSA; i++){
		std::cout<<i<<" "<<1.64*sqrt(COV_ADSA[i][i])<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Variances root delta alpha (90%):"<<std::endl;
	if(params.NSflag2){
		std::cout<<"Cutoff:"<<params.mass1*.5*1.5<<std::endl;
	}
	else{
		std::cout<<"Cutoff:"<<params.mass2*.5*1.5<<std::endl;
	}
	std::cout<<"(delta alpha^2)^(1/4) (KM): "<<1.64*pow(COV_ADSA[dimDSA-1][dimDSA-1],1./8.)*3.e5<<std::endl;
	std::cout<<std::endl;
	std::cout<<"SNR (SA): "<<sqrt(output_ADSA[0][0])<<std::endl;

	deallocate_2D_array(output_AD,dim,dim);
	deallocate_2D_array(output_AD_temp,dim,dim);
	deallocate_2D_array(COV_AD,dim,dim);
	deallocate_2D_array(COV_AD2,dim,dim);
	deallocate_2D_array(COV_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA,dimDSA,dimDSA);
	deallocate_2D_array(output_ADSA_temp,dimDSA,dimDSA);
	deallocate_2D_array(output_AD2,dimD,dimD);
	deallocate_2D_array(output_AD2_temp,dimD,dimD);
	deallocate_2D_array(jac_spins,dim,dim);
	
	delete [] frequency;
	for(int i = 0 ; i<Ndetect; i++){
		delete [] psd[i];
	}
	delete [] params.betappe;
	delete [] params.bppe;
	//delete [] params.phii;
	//delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;
}
int AD_v_N(int argc, char *argv[])
{
	
	std::cout.precision(5);
	gen_params params;	
	params.spin1[2] = .1;
	params.spin2[2] = -.1;
	params.chip = .03;
	params.phip = 1.0;
	params.Luminosity_Distance = 1000;
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

	double fmin = 5;
	double fmax = 2048;
	//double fmin = f_0PN(4*T_year,chirpmass);
	//double fmax = 1;
	double T = 32;
	//double T = 4*T_year;

	params.tc = 3.*T/4.;
	//params.tc = T;
	//params.tc = 0;

	//int length = T*(fmax-fmin);
	int length = 1000;
	double *frequency = new double[length];
	//int Ndetect = 4;
	int Ndetect = 2;
	double **psd = new double*[Ndetect];
	//std::string SN[4] = {"CE2_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	std::string SN[3] = {"AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
	//std::string SN[4] = {"LISA_SADC_CONF","AdLIGODesign_smoothed","AdLIGODesign_smoothed","AdLIGODesign_smoothed"};
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
	//std::string method = "ppE_IMRPhenomPv2_Inspiral";
	//std::string method = "ppE_IMRPhenomD_Inspiral";
	//std::string method = "gIMRPhenomD";
	std::string method = "gIMRPhenomPv2";
	//transform_orientation_coords(&params, method, detector);

	//params.Nmod = 1;
	//params.betappe = new double[1];
	//params.bppe = new int[1];
	//params.betappe[0] = 0;
	//params.bppe[0] = -1.;
	//params.bppe[0] = -1;
	//params.bppe[0] = 1;
	//
	params.Nmod_phi = 1;
	params.delta_phi = new double[1];
	params.phii = new int[1];
	params.delta_phi[0] = 0;
	params.phii[0] = 2;
	std::cout<<"Phase power: "<<params.phii[0]<<std::endl;

	//params.Nmod_alpha = 1;
	//params.delta_alpha = new double[1];
	//params.alphai = new int[1];
	//params.delta_alpha[0] = 0;
	//params.alphai[0] = 2;

	//std::string detectors[4] = {"CE","Hanford","Livingston","Virgo"};
	std::string detectors[3] = {"Hanford","Livingston","Virgo"};
	//std::string detectors[4] = {"LISA","Hanford","Livingston","Virgo"};
		
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

	//params.equatorial_orientation = true;
	params.equatorial_orientation = false;
	params.sky_average = true;
	//params.theta_l = 0.12;
	//params.phi_l = 0.1;
	//params.incl_angle  = 0;
	double snr;
	//snr = calculate_snr(SN[0],"LISA",method, &params, frequency, length, "GAUSSLEG",weights,true);
	//std::cout<<snr<<std::endl;
	//double SNR_TARGET = 100;
	//params.Luminosity_Distance = snr/SNR_TARGET*params.Luminosity_Distance;
	params.sky_average = false;
	params.incl_angle = M_PI-.51;



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
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD_temp, dim, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dim; k++){
			for(int j = 0 ; j<dim; j++){
				output_AD[k][j]+= output_AD_temp[k][j];
			}
		}
	}
	matrix_multiply(output_AD, jac_spins,output_AD_temp,dim,dim,dim);
	matrix_multiply(jac_spins,output_AD_temp, output_AD,dim,dim,dim);
	std::cout<<"SNR: "<<sqrt(output_AD[6][6])<<std::endl;
	snr = calculate_snr(SN[0],"CE",method, &params, frequency, length, "GAUSSLEG",weights,true);
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

	//method = "ppE_IMRPhenomD_Inspiral";
	method = "gIMRPhenomD";
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[0], output_AD3_temp, dimD, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
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

	//method = "ppE_IMRPhenomD_Inspiral";
	//method = "IMRPhenomD";
	//dimDSA = 7;
	params.sky_average = true;
	params.incl_angle = 0;
	for(int i = 0 ;i < Ndetect; i++){
		fisher_autodiff(frequency, length, method, detectors[i],detectors[i], output_ADSA_temp, dimDSA, &params, "GAUSSLEG",weights,true, psd[i],NULL,NULL);
		for(int k = 0 ; k<dimDSA; k++){
			for(int j = 0 ; j<dimDSA; j++){
				output_ADSA[k][j]+= output_ADSA_temp[k][j];
			}
		}
	}
	
	
	std::cout<<"AD-DSA:"<<std::endl;
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
	//delete [] params.betappe;
	//delete [] params.bppe;
	delete [] params.phii;
	delete [] params.delta_phi;
	delete [] psd;
	delete [] weights;
	return 0;
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Compare AD to Numerical"<<std::endl;
	std::cout<<"1 --- Network Fishers"<<std::endl;
	std::cout<<"2 --- dCS or EdGB"<<std::endl;
	std::cout<<"3 --- Check ppE-theory jac transform"<<std::endl;
}


