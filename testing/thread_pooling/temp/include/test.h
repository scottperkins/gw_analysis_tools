#ifndef TEST_H
#define TEST_H
//#include "types.h"
//#include <adolc/adouble.h>
////#include "wrapper.h"
////#include "wrapper2.h"
//#if type == type_double
////#include "data_type_double.h"
//#define data_type double
//#endif 
////
//#if type == type_adouble
////#include "data_type_adouble.h"
//#define data_type adouble
//#endif
//if (type == "double")
//#define data_type type
#include <adolc/adouble.h>
class test{
public:
#define data_type adouble
data_type c;
data_type func(data_type a, data_type b );
#undef data_type
#define data_type double
data_type func(data_type a, data_type b );
};
#endif 
