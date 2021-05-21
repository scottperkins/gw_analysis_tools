#include <gwat/psd_estimation.h>
#include <gwat/util.h>
#include <gwat/io_util.h>
#include <iostream>
#include <string>

int test_bayesline_LL(int argc, char *argv[]);
int test_bayesline(int argc, char *argv[]);
int RT_ERROR();
int main(int argc, char *argv[]){
	if(argc != 2){

		return RT_ERROR();
	}
	int runtime = std::atoi(argv[1]);
	if(runtime == 0){
		return test_bayesline(argc, argv);
	}
	else if(runtime == 1){
		return test_bayesline_LL(argc, argv);
	}
	else{
		return RT_ERROR();
	}
	return 0;
}

double FFT_data(double *data, int length, double T,double *freqs,std::complex<double> *data_out)
{
	double *window = new double[length];
	double df = 1/T;
	double alpha = 2.*.4/T;
	
	tukey_window(window, length, alpha);
	//for(int i = 0 ; i<length; i++){
	//	window[i] = 1;
	//}

	fftw_outline plan;
	allocate_FFTW_mem_forward(&plan, length);
	for (int i = 0 ; i<length; i++){
		plan.in[i][0] = data[i]*window[i];
		plan.in[i][1] = 0;
	}
	fftw_execute(plan.p);
	for(int i = 0 ; i<length; i++){
		data_out[i] = std::complex<double>(plan.out[i][0],plan.out[i][1]);
		freqs[i] = i*df;
	}
	deallocate_FFTW_mem(&plan);
	delete [] window;
}

int test_bayesline(int argc, char *argv[])
{
	std::string input_file("data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt");
	int length_in = 0;
	count_lines_LOSC_data_file(input_file, &length_in);

	double start_time, duration, fs;
	double *data_in = new double[length_in];
	read_LOSC_data_file(input_file, data_in, &start_time, &duration, &fs);
	double T = 4;
	int init_id = 0;
	int final_id = init_id + int(fs*T);
	int temp_length=final_id - init_id;
	
	int length = temp_length;
	double *freqs= new double[length];
	std::complex<double> *data_stream = new std::complex<double>[length];
	
	double FMIN = 50;
	//double FMAX = 2048;
	double FMAX = 512;
	//double FMAX = 20;
	FFT_data(&(data_in[init_id]), temp_length, T, freqs, data_stream);
	
	length = 0 ; 
	int initial_id = -1;
	for(int i = 0 ; i<temp_length ; i++){
		if(freqs[i]>FMIN && freqs[i]<FMAX){
			if(initial_id <0){
				initial_id = i;
			}
			
			length++;
		}
	}
	std::cout<<"Length: "<<length<<std::endl;
	double **output_RAW = new double*[length];
	for(int i = 0 ; i<length; i++){
		output_RAW[i] = new double[3];
		output_RAW[i][0] = freqs[initial_id + i];
		output_RAW[i][1] = std::real(data_stream[initial_id+i]);
		output_RAW[i][2] = std::imag(data_stream[initial_id+i]);
	}
	write_file("data/bayesline_raw_data.csv",output_RAW,length,3);
	for(int i = 0 ; i<length; i++){
		delete [] output_RAW[i];	
	}

	delete [] output_RAW;
	

	int N_L_MAX = 20;
	int N_L_MIN = 2;
	int N_S_MAX = 100;
	int N_S_MIN = 2;

	std::string chain_allocation("double");
	int samples = 1e4;
	int chain_N = 60;
	double chain_temps[chain_N];
	int ensemble_chain_N = 15;
	int swap_freq = 20;
	int t0 = 2000;
	int nu = 100;
	int max_chunk_size = samples/10;
	int threads = 8;
	bool pool = true;
	bool show_prog = true;
	std::string chain_file("data/bayesline_output.hdf5");
	std::string stat_file("data/bayesline_stat.txt");
	std::string checkpoint_file("data/bayesline_checkpoint.csv");
	PSD_output output;
	output.SN = new double[length];
	output.length = length;
	output.T = T;

	bayesline_psd_estimation(&(data_stream[initial_id]), &(freqs[initial_id]), length, T,N_L_MIN,N_L_MAX,N_S_MIN,N_S_MAX,(double*)NULL, (int*)NULL,(double*)NULL,chain_allocation, samples, chain_N, ensemble_chain_N, chain_temps, swap_freq, t0,nu,max_chunk_size, threads, pool, show_prog, chain_file, stat_file, checkpoint_file, &output);
	write_file("data/bayesline_estimate_psd.csv",output.SN, length);

	delete [] freqs;
	delete [] data_stream;
	delete [] data_in;
	delete [] output.SN;
	return 0 ;
}
int test_bayesline_LL(int argc, char *argv[])
{
	int NS_MAX = 5;
	int NL_MAX = 5;
	int NS_MIN = 2;
	int NL_MIN = 2;
	return 0;
}

int RT_ERROR(){
	std::cout<<"Runtime options: "<<std::endl;
	std::cout<<" 0 -- test bayesline"<<std::endl;
	std::cout<<" 1 -- test bayesline LL"<<std::endl;
}

