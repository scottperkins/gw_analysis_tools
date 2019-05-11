#include "test3.h"
#include <adolc/adouble.h>
template <class data_type>
data_type test3<data_type>::function(data_type a)
{
	return a*a;
}

template class test3<double>;
template class test3<adouble>;
