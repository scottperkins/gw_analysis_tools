#ifndef IMRPHENOMP_H
#define IMRPHENOMP_H
#include "IMRPhenomD.h"
#include "util.h"


template<class T> 
class IMRPhenomPv2: public IMRPhenomD<T>
{
public:

virtual T alpha(T omega, T q,T chi2l, T chi2);

virtual T epsilon(T omega, T q, T chi2l, T chi2);

virtual T d(int l, int mp, int m,T s);
};

#endif
