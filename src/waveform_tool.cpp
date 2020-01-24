#include <gwat/waveform_util.h>
#include <gwat/io_util.h>
#include <gwat/util.h>
#include <iostream>

int main(int argc, char *argv[])
{
	
	if(argc !=2){
		std::cout<<"ERROR -- A parameter file is required"<<std::endl;
		return 1;
	}
	std::string param_file(argv[1]);

	std::unordered_map<std::string, int> int_dict;
	std::unordered_map<std::string, double> dbl_dict;
	std::unordered_map<std::string, float> flt_dict;
	std::unordered_map<std::string, bool> bool_dict;
	std::unordered_map<std::string, std::string> str_dict;
	std::string input_param_file(argv[1]);
	int status = unpack_input_io_file(input_param_file, &int_dict, &str_dict, &dbl_dict, &flt_dict, &bool_dict);
	gen_params params;
	params.mass1 = dbl_dict["Mass 1"];
	params.mass2 = dbl_dict["Mass 2"];
	params.spin1[0] = dbl_dict["Spin 1 x"];
	params.spin1[1] = dbl_dict["Spin 1 y"];
	params.spin1[2] = dbl_dict["Spin 1 z"];
	params.spin2[0] = dbl_dict["Spin 2 x"];
	params.spin2[1] = dbl_dict["Spin 2 y"];
	params.spin2[2] = dbl_dict["Spin 2 z"];
	params.Luminosity_Distance = dbl_dict["Luminosity Distance"];
	params.phiRef = dbl_dict["phiRef"];
	params.f_ref = dbl_dict["Reference Frequency"];
	params.tc = dbl_dict["tc"];
	params.RA = dbl_dict["Right Ascension"];
	params.DEC = dbl_dict["Declination"];
	params.gmst = dbl_dict["GMST"]*M_PI/12;
	params.NSflag1 = bool_dict["NSflag1"];
	params.NSflag2 = bool_dict["NSflag2"];
	params.equatorial_orientation = bool_dict["Equatorial Orientation"];
	params.sky_average = bool_dict["Sky Average"];
	std::string detector = str_dict["Detector"];
	std::string gen_meth = str_dict["Generation Method"];
	std::string file = str_dict["Output File Path"];

	if(!params.equatorial_orientation){
		params.psi = dbl_dict["Polarization Angle"];
		params.incl_angle = dbl_dict["Inclination Angle"];
	}
	else{
		params.theta_l = dbl_dict["Equatorial polar angle of L"];
		params.phi_l = dbl_dict["Equatorial azimuthal angle of L"];
	}

	double T = dbl_dict["Observation Time"];
	if(detector=="LISA"){
		T *=T_year;	
		std::cout<<1./T<<std::endl;
	}
	double fmin = dbl_dict["Minimum Frequency"];
	double fmax = dbl_dict["Maximum Frequency"];
	double delta_f = 1./T;
	int length = (fmax-fmin)*T;
	std::cout<<length<<std::endl;
	double *freqs = new double[length];
	std::complex<double> *response = new std::complex<double>[length];
	double *times = NULL;
	for(int i = 0 ; i<length; i++){
		freqs[i]=fmin + i * delta_f;
	}

	if(detector=="LISA"){
		times = new double[length];
		time_phase_corrected_autodiff(times, length, freqs, &params, gen_meth, false,NULL);	
		params.LISA_alpha0 =  0;
		params.LISA_phi0 =  0;
		
	}

	fourier_detector_response(freqs, length, response, detector, gen_meth, &params,times);
	double **output = allocate_2D_array(length, 3);
	for(int i = 0 ; i<length; i++){
		output[i][0]=freqs[i];
		output[i][1]=std::real(response[i]);
		output[i][2]=std::imag(response[i]);
	}
	write_file(file, output, length, 3);
	delete [] freqs;
	delete [] response;
	if(detector=="LISA"){
		delete[] times;
	}
	deallocate_2D_array(output,length,3);	
	return 0;
}
