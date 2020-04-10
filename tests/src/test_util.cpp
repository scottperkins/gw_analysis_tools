#include <iostream>
#include <math.h>
#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <gwat/io_util.h>

#include <lal/TimeDelay.h>
#include <lal/LALDetectors.h>

int test_newton_raphson(int argc, char *argv[]);
void cos_sin(double x, double *func, double *func_prime, void *param);
int HDF5_testing_write(int argc,char *argv[]);
int time_delay_testing(int argc,char *argv[]);
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
	else if(runtime_opt == 2){
		return time_delay_testing(argc,argv);	
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
int time_delay_testing(int argc,char *argv[])
{
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
}
