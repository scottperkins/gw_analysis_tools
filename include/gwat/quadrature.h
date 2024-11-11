#ifndef QUADRATURE_H
#define QUADRATURE_H


#include "util.h"


//!\class Quadrature
//!\brief Class to evaluate integrals with established spacing and weights.
class Quadrature
{
protected:
    // Length of integrand
    int length = 0;

public:
    ~Quadrature() = default;

    virtual double integrate(double *integrand) = 0;
    virtual int get_length() {return length;}
};

//! \brief Simpson's rule for uniformly-spaced integrals.
//!
//! Quadrature with the classic extended 3-point rule
//! (see, e.g., Numerical Recipes, extended Simpsons rule).
//! For even lengths, the trapezoidal rule is used at the last interval.
class SimpsonsQuad : public Quadrature
{
private:
    double del;         //< Overall factor h/3
    bool evenFlag = false;  //< Track if the length is even
    int SimpsonsEnd;    //< End index of integrand array for Simpson's
    double TrapDel;     //< Overall factor for trapezoid h/2

public:
    SimpsonsQuad(int length, double delta);

    double integrate(double *integrand);
};


#endif // QUADRATURE_H