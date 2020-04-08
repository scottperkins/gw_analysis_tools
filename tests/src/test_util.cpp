#include <iostream>
#include <math.h>
#include <gwat/util.h>
#include <gwat/io_util.h>

int test_newton_raphson(int argc, char *argv[]);
void cos_sin(double x, double *func, double *func_prime, void *param);
int HDF5_testing_write(int argc,char *argv[]);
void RT_ERROR_MSG();
int main(int argc, char *argv[])
{
	std::cout<<"TESTING UTILITY FUNCTIONS"<<std::endl;
	if(argc != 2){
		RT_ERROR_MSG();
		return 1;
	}
	
	int runtime_opt = std::stoi(argv[1]);	
	if(runtime_opt == 0){
		return test_newton_raphson(argc,argv);	
	}
	else if(runtime_opt == 1){
		return HDF5_testing_write(argc,argv);	
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int HDF5_testing_write(int argc,char *argv[])
{
	//HDF5wrapper io();
	return 0;
}
int test_newton_raphson(int argc, char *argv[])
{
	std::cout.precision(15);
	double tol = 1e-6;
	double initial = M_PI/3.;
	int max_iter = 100;
	double solution;
	newton_raphson_method_1d(&cos_sin, initial, tol,max_iter, (void *)NULL, &solution);
	std::cout<<"Final Result: "<<solution<<std::endl;
	std::cout<<"Final absolute error: "<<M_PI/2. - solution<<std::endl;
	
	return 0 ;
}
void cos_sin(double x, double *func, double *func_prime, void *param)
{
	*func = cos(x);
	*func_prime= - sin(x);
}
void RT_ERROR_MSG()
{
	std::cout<<"ERROR -- incorrect arguments"<<std::endl;
	std::cout<<"Please supply function option:"<<std::endl;
	std::cout<<"0 --- Test newton raphson"<<std::endl;
	std::cout<<"1 --- Test HDF5 file write"<<std::endl;
}
