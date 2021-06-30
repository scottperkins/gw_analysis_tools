#include <iostream>
#include <math.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>

//#define _LAL
#ifdef _LAL
	#include <lal/TimeDelay.h>
	#include <lal/LALDetectors.h>
#endif

int test_newton_raphson(int argc, char *argv[]);
void cos_sin(double x, double *func, double *func_prime, void *param);
int HDF5_testing_write(int argc,char *argv[]);
int time_delay_testing(int argc,char *argv[]);
int matrix_multiplication(int argc,char *argv[]);
int vector_union_test(int argc,char *argv[]);
int mvn_testing(int argc,char *argv[]);
int GPS_to_GMST(int argc,char *argv[]);
void RT_ERROR_MSG();
int main(int argc, char *argv[])
{
	std::cout<<"TESTING UTILITY FUNCTIONS"<<std::endl;
	if(argc < 2){
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
	else if(runtime_opt == 2){
		return time_delay_testing(argc,argv);	
	}
	else if(runtime_opt == 3){
		return matrix_multiplication(argc,argv);	
	}
	else if(runtime_opt == 4){
		return vector_union_test(argc,argv);	
	}
	else if(runtime_opt == 5){
		return GPS_to_GMST(argc,argv);	
	}
	else if(runtime_opt == 6){
		return mvn_testing(argc,argv);	
	}
	else{
		RT_ERROR_MSG();
		return 1;
	}
}
int mvn_testing(int argc,char *argv[])
{
	int samples = 10000;
	int dim = 5;
	double **cov = new double*[dim];
	double mean[dim] ;
	for(int i = 0 ; i<dim ; i++){
		mean[i] = 2*i;
		cov[i] = new double[dim];
		for(int j = 0 ; j<dim ; j++){
			if(j == i){
				cov[i][j] = 1*(1+j);

			}
			else{
				cov[i][j] = pow_int(-1,(i+j))*.9*(i+j)/2;
			}
		}
	}
	
	double **output = allocate_2D_array(samples, dim );
	mvn_sample(samples, mean, cov, dim , output);
	write_file("data/test_mvn.csv",output, samples, dim);
	write_file("data/test_mvn_cov.csv",cov, dim, dim);
	write_file("data/test_mvn_mean.csv",mean,  dim);
	deallocate_2D_array(output,samples,dim);
	for(int i = 0 ; i<dim ; i++){
		delete [] cov[i];
	}
	delete [] cov;
	return 0;

}
int GPS_to_GMST(int argc,char *argv[])
{
	if (argc != 3){
		std::cout<<"Supply GPS time!"<<std::endl;
		return 1;
	}
	double gps = std::stod(argv[2]);	
		
	std::cout<<"In radian: "<<gps_to_GMST_radian(gps)<<std::endl;
	std::cout<<"In hours: "<<gps_to_GMST(gps)<<std::endl;
	return 0;

}
int vector_union_test(int argc,char *argv[])
{
	std::vector<double> A = {1,3,8,9};
	std::vector<double> B = {3,4,5,10};
	std::vector<double> C;
	vector_union(A, B, &C);
	for(auto i = C.begin(); i!=C.end(); i++){
		std::cout<<*i<<std::endl;
	}
	return 0;
}
int matrix_multiplication(int argc,char *argv[])
{
	int dim1=4;
	int dim2=3;
	int dim3=5;
	double **A =allocate_2D_array(dim1,dim2);
	double **B =allocate_2D_array(dim2,dim3);
	double **C =allocate_2D_array(dim1,dim3);
	for(int i =0 ; i<dim1; i++){
		for(int j = 0 ; j<dim2; j++){
			A[i][j] = (i+1)*(j+2);	
		}
	}
	for(int i =0 ; i<dim2; i++){
		for(int j = 0 ; j<dim3; j++){
			B[i][j] = (i+1)*(j+2);	
		}
	}
	matrix_multiply(A,B,C,dim1,dim2,dim3);
	std::cout<<"A: "<<std::endl;
	for(int i =0 ; i<dim1; i++){
		for(int j = 0 ; j<dim2; j++){
			std::cout<<A[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"B: "<<std::endl;
	for(int i =0 ; i<dim2; i++){
		for(int j = 0 ; j<dim3; j++){
			std::cout<<B[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"C: "<<std::endl;
	for(int i =0 ; i<dim1; i++){
		for(int j = 0 ; j<dim3; j++){
			std::cout<<C[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	
	deallocate_2D_array(A,dim1,dim2);	
	deallocate_2D_array(B,dim2,dim3);	
	deallocate_2D_array(C,dim1,dim3);	
	return 0 ;
}
int HDF5_testing_write(int argc,char *argv[])
{
	//HDF5wrapper io();
	return 0;
}
int time_delay_testing(int argc,char *argv[])
{
#ifdef _LAL
	std::cout.precision(15);
	LALDetector LALD1 = lalCachedDetectors[6];		
	LALDetector LALD2 = lalCachedDetectors[5];		
	double gps_time = 1126259462.4;
	//double gps_time = 1126259462.8984298429;
	const LIGOTimeGPS *gps_timeLAL;
	LIGOTimeGPS gps_timeLAL_temp;
	gps_timeLAL_temp.gpsSeconds= (int)gps_time;
	gps_timeLAL_temp.gpsNanoSeconds= (int)((gps_time - (int)gps_time)*pow_int(10,9));
	gps_timeLAL = &gps_timeLAL_temp;
	double RA = 2.;
	double DEC = -1.4;
	double DTOA = XLALArrivalTimeDiff(LALD1.location,LALD2.location,RA,DEC,gps_timeLAL);	
	std::cout<<"LAL DTOA: "<<DTOA<<std::endl;
	std::cout<<gps_timeLAL->gpsSeconds+gps_timeLAL->gpsNanoSeconds*pow_int(10,-9)<<std::endl;
	std::cout<<gps_timeLAL->gpsNanoSeconds<<std::endl;

	double DTOA_GWAT = DTOA_DETECTOR(RA,DEC,gps_to_GMST_radian(gps_time),"Hanford","Livingston");
	std::cout<<"GWAT DTOA: "<<DTOA_GWAT<<std::endl;
	std::cout<<"fractional diff: "<<(DTOA_GWAT - DTOA)/(DTOA)<<std::endl;

#endif
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
	std::cout<<"2 --- Test Time Delay "<<std::endl;
	std::cout<<"3 --- Matrix Multiplication "<<std::endl;
	std::cout<<"4 --- Vector Union "<<std::endl;
	std::cout<<"5 --- GPS to GMST "<<std::endl;
	std::cout<<"6 --- MVN sampling test "<<std::endl;
}
