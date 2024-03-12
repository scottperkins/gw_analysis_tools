#ifndef IMRPHENOMPV3_H
#define IMRPHENOMPV3_H

#include "IMRPhenomD.h"
#include "IMRPhenomP.h"
#include "util.h"
#include "IMRPhenomPv3utils.h"


/*! \file
 * Header file for IMRPhenomPv3 utilities
 */


template<class T>
class IMRPhenomPv3: public IMRPhenomPv2<T>
{
public:

};

template <class T> static void IMRPhenomPv3_Compute_a_b_e(
    T *alpha, T *beta, T *two_epsilon,
    T fHz, const T pi_Msec,
    PhenomPv3Storage<T> *params, sysprecquant<T> *pAngles);

#endif /* IMRPHENOMPV3_H */