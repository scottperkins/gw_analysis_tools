#include <iostream>
#include <gwat/ppE_utilities.h>
#include <gwat/D_Z_Config_modified_dispersion.h>
#include <gwat/io_util.h>

int MD_testing(int argc,char*argv[]);
void RT_ERROR_MSG();
int main(int argc, char *argv[])
{
	std::cout<<"TESTING COSMOLOGY CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"Modified Dispersion testing"<<std::endl;
		return MD_testing(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int MD_testing(int argc,char*argv[])
{
	int num_segments=3;
	const double *alphas = MD_alphas;
	int index = 0;
	int samples = 100;
	double** zs=new double*[2];	
	zs[0] = new double[num_segments*samples];
	zs[1] = new double[num_segments*samples];
	for(int k = 0 ; k<num_MD_alphas; k++){
		index = k;
		double alpha = alphas[k];
		const double *boundaries=MD_boundaries_Z[index];
		
		for(int i=0; i<num_segments; i++){
			zs[0][i*samples]=boundaries[i];
			double deltaz = pow(boundaries[i+1]/boundaries[i],1./samples);
			for(int j = 0 ; j<samples-1; j++){
				zs[0][i*samples + j+1] = zs[0][i*samples+j]*deltaz;
			}
			for(int j = 0 ; j<samples; j++){
				zs[1][i*samples+j] = DL_from_Z_MD(zs[0][i*samples+j],alpha);
			}
		}
		write_file("data/MD_DL_from_Z_"+std::to_string(alpha)+".csv",zs,2,num_segments*samples);	
	}
	//Cleanup
	delete [] zs[0];
	delete [] zs[1];
	delete [] zs;
	return 0;

}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Modified Dispersion testing"<<std::endl;
}
