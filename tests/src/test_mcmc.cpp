#include <gwat/waveform_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/io_util.h>
#include <gwat/mcmc_sampler.h>
#include <iostream>


void RT_ERROR_MSG();
int mcmc_standard_test(int argc, char *argv[]);
double log_test (double *c,int dim,void *parameters);
double log_test_prior (double *c,int dim,void *parameters);
void fisher_test(double *c,int dim, double **fisher,  void *parameters);

int main(int argc, char *argv[])
{
	std::cout<<"TESTING MCMC CALCULATIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		std::cout<<"MCMC student t"<<std::endl;
		return mcmc_standard_test(argc,argv);
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int mcmc_standard_test(int argc, char *argv[])
{
	int dimension = 2;
	double initial_pos[2]={1,0.};
	double *seeding_var = NULL;
	int N_steps = 5000;
	int chain_N= 10;
	double ***output;
	output = allocate_3D_array( chain_N, N_steps, dimension );
	//double *initial_pos_ptr = initial_pos;
	int swp_freq = 3;
	//double chain_temps[chain_N] ={1,2,3,10,12};
	double chain_temps[chain_N];
	chain_temps[0] = 1.;
	double c = 1.5;
	for(int i =1; i < chain_N;  i ++)
		chain_temps[i] =  chain_temps[i-1] * c;
	std::string autocorrfile = "";
	std::string chainfile = "data/mcmc_output.csv";
	std::string statfilename = "data/mcmc_statistics.txt";
	std::string checkpointfile = "data/mcmc_checkpoint.csv";
	
	int numThreads = 10;
	bool pool = true;
	bool show_progress = true;
	
	PTMCMC_MH(output, dimension, N_steps, chain_N, initial_pos,seeding_var,chain_temps, swp_freq, log_test_prior, log_test,fisher_test,(void **)NULL,numThreads, pool,show_progress, statfilename,chainfile,autocorrfile, "",checkpointfile );	
	std::cout<<"ENDED"<<std::endl;


	deallocate_3D_array(output, chain_N, N_steps, dimension);
		
	return 0;

}
double log_test (double *c,int dim,void *parameters)
{
	double x = c[0];
	double y = c[1];
	double prefactor = 16./(M_PI*3.);
	double pow1 = -x*x - pow((9+4*x*x +8*y),2);
	double pow2 = -8*x*x -8*pow(y-2,2);
	return log(prefactor*(std::exp(pow1) + .5*std::exp(pow2)));
	//return 2.;
}
double log_test_prior (double *c,int dim,void *parameters)
{
	return 1.;
}

double fisher11(double *c)
{
	double x = c[0];
	double y = c[1];

	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(x,2) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(-2 - 128*pow(x,2) - 16*(9 + 4*pow(x,2) + 8*y)) + 
     exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y),2));
}
double fisher12(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*x + 128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*x*(-2 + y) - 
     16*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*(9 + 4*pow(x,2) + 8*y)*(-2*x - 16*x*(9 + 4*pow(x,2) + 8*y)));

}
double fisher22(double *c)
{
	double x = c[0];
	double y = c[1];
	return 1.6976527263135504*(-8.*exp(-8*pow(x,2) - 8*pow(-2 + y,2)) - 128*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2)) + 
     128.*exp(-8*pow(x,2) - 8*pow(-2 + y,2))*pow(-2 + y,2) + 256*exp(-pow(x,2) - pow(9 + 4*pow(x,2) + 8*y,2))*pow(9 + 4*pow(x,2) + 8*y,2));

}
void fisher_test(double *c,int dim, double **fisher,  void *parameters)
{
	
	fisher[0][0] = fisher11(c);
	fisher[1][1] = fisher22(c);
	fisher[1][0] = fisher12(c);
	fisher[0][1] = fisher12(c);
	return;
} 
	

void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Standard"<<std::endl;
}
