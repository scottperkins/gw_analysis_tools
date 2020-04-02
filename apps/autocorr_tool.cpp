#include <iostream>
#include <string>
#include <gwat/autocorrelation.h>

int main( int argc, char * argv[])
{
	if(argc != 5){
		std::cout<<"ERROR -- wrong number of inputs"<<std::endl;
		std::cout<<"Inputs are:"<<std::endl;
		std::cout<<"Input file"<<std::endl;
		std::cout<<"Output file"<<std::endl;
		std::cout<<"dimension"<<std::endl;
		std::cout<<"Number of Segments"<<std::endl;
	}
	
	std::string input_file(argv[1]);
	std::string output_file(argv[2]);
	int dim = atoi(argv[3]);
	int segs = atoi(argv[4]);

	
	return 0;
}
