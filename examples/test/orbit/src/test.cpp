#include <gwat/util.h>
#include <gwat/detector_util.h>
#include <iostream>
#include <math.h>
#include <fstream>

int main(){
	
	double time_step;
	time_step = 24.0*3600.;
	double time_max;
	time_max = 365*time_step;
	int N_time = int(time_max/time_step);
	
	std::cout << "cp0" << std::endl;
	
	double *t_array = new double[N_time];
	double **p0 = new double*[N_time];
	for(int i=0;i<N_time;i++)
	{
		p0[i]= new double[3];
	}
	
	
	std::cout << "cp1" << std::endl;
	
	for(int i; i<N_time; i++){
		t_array[i] = i*time_step;
	}
	
	
	
	funcp0(t_array, p0, N_time);
	
	std::ofstream out_file;
	out_file.open("p0.txt");
	out_file.precision(15);
	
	std::cout << "cp2" << std::endl;

	for(int i; i<N_time; i++){
		out_file << t_array[i] << " " << p0[i][0] << " " << p0[i][1] << " " << p0[i][2] << std::endl;
	}
	
	out_file.close();
	
	return 0;
	


} 
