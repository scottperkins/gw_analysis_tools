#include "quadrature.h"


SimpsonsQuad::SimpsonsQuad(
    int length,
    double delta
)
{
    del = delta/3.;
    this->length = length;

    // Check if there are not enough samples
    if (length < 3)
    {
        std::cerr << "Not enough points for Simpson's rule. "
         << "Result may be inaccurate.\n";
    }

    if ((length % 2) == 0)
    {
        // If length is even, setup trapezoidal rule for the last interval
        evenFlag = true;
        SimpsonsEnd = length-2;
        TrapDel = 0.5*delta;
    }
    else
    {
        // Otherwise the whole integrand can be integrated with Simpson's
        SimpsonsEnd = length-1;
    }
}

double SimpsonsQuad::integrate(double *integrand)
{
    // Integrand array index
    int i = 1;
    // Result of sum
    double integral = 4.*integrand[i++]; // First term of inner sum

    // Inner sum, 2*f(k) + 4*f(k+1) for even k
    while (i < SimpsonsEnd)
    {
        integral += 2.*integrand[i++];
        integral += 4.*integrand[i++];
    }

    // Add in endpoints
    integral += integrand[0] + integrand[SimpsonsEnd];
    // Multiply by overall factor
    integral *= del;

    if (evenFlag)
    {
        // Integrate final interval with trapezoidal rule.
        // Inaccurate compared to Simpson's,
        // but should be enough for small enough spacing.
        integral += TrapDel*(integrand[SimpsonsEnd] + integrand[SimpsonsEnd+1]);
    }

    return integral;
}