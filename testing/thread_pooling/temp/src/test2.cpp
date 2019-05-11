#include "test.h"
#include <adolc/adouble.h>
#undef data_type
#define data_type adouble
#include "test.cpp"
