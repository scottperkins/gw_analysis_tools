#ifndef TEST2_H
#define TEST2_H
#include "adolc/adouble.h"
#include "test.h"
#define data_type adouble

class test2: public test
{
public:
data_type func2(data_type a);
};
#endif
