#include <error.h>
#include <math.h>
#include "fisher.h"
#include "util.h"
#include "detector_util.h"
#include "waveform_generator.h"
#include "waveform_util.h"

using namespace std; 

/* This is a file to calculate statistical and systematic error. See fisher.cpp (or fisher.h) for a description of the function fisher_numerical (which is what we should use in these calculations). See util.h for a description of the gen_params_base and all of the different parameters it contains (this is what we will use to pass parameters throughout the code). The file tests/src/test_fishers.cpp gives a good example of how to call the fisher function and how to specify the waveform template used with the string "generation_method".  */
