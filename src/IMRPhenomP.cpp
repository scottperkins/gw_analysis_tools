#include "IMRPhenomP.h"
#include "IMRPhenomD.h"
#include "util.h"
#include <adolc/adouble.h>
#include <math.h>
#include <algorithm>



template<class T>
T IMRPhenomPv2<T>::alpha(T omega, T q,T chi2l, T chi2){

	T alpha;
//alpha = ((-35.*chi2l*q)/128. - (15.*chi2l*(1 + q)*(1./(1 + q) - q/(1. + q)))/128.)/pow(omega,0.6666666666666666) + 
//   (-0.18229166666666666 - (5.*(1 + q)*(1./(1 + q) - q/(1 + q)))/(64.*q))/omega + 
//   (-1.7952473958333333 - (35.*pow(chi2,2)*pow(q,2))/128. - (515.*q)/(384.*pow(1 + q,2)) - 
//      (175.*(1./(1 + q) - q/(1 + q)))/(256.*(1 + q)) - (4555.*(1 + q)*(1/(1 + q) - q/(1 + q)))/(7168.*q) - 
//      (15.*pow(chi2,2)*q*(1 + q)*(1./(1 + q) - q/(1 + q)))/128. - 
//      (15.*pow(1./(1 + q) - q/(1 + q),2))/(256.*q))/pow(omega,0.3333333333333333) + 
//   pow(omega,0.3333333333333333)*(4.318908476114694 - (35.*chi2l*M_PI*q)/16. + 
//      (475.*pow(chi2,2)*pow(q,2))/6144. - (35.*pow(chi2,4)*pow(q,4))/512. + 
//      (35.*pow(chi2,2)*pow(chi2l,2)*pow(q,4))/128. + (955.*pow(q,2))/(576.*pow(1 + q,4)) + 
//      (39695.*q)/(86016.*pow(1 + q,2)) + (575.*pow(chi2,2)*pow(q,3))/(1536.*pow(1 + q,2)) + 
//      (1645.*pow(chi2l,2)*pow(q,3))/(192.*pow(1 + q,2)) + 
//      (2725.*q*(1/(1 + q) - q/(1 + q)))/(3072.*pow(1 + q,3)) - 
//      (265.*(1/(1 + q) - q/(1 + q)))/(14336.*(1 + q)) + 
//      (145.*pow(chi2,2)*pow(q,2)*(1./(1 + q) - q/(1 + q)))/(512.*(1 + q)) + 
//      (1815.*pow(chi2l,2)*pow(q,2)*(1./(1 + q) - q/(1 + q)))/(256.*(1 + q)) - 
//      (15.*chi2l*M_PI*(1 + q)*(1./(1 + q) - q/(1 + q)))/16. + 
//      (27895885.*(1 + q)*(1./(1 + q) - q/(1 + q)))/(2.1676032e7*q) - 
//      (485.*pow(chi2,2)*q*(1 + q)*(1./(1 + q) - q/(1. + q)))/14336. - 
//      (15.*pow(chi2,4)*pow(q,3)*(1 + q)*(1./(1 + q) - q/(1 + q)))/512. + 
//      (15.*pow(chi2,2)*pow(chi2l,2)*pow(q,3)*(1 + q)*(1./(1 + q) - q/(1 + q)))/128. + 
//      (1615.*pow(1./(1 + q) - q/(1 + q),2))/(28672.*q) + 
//      (15.*pow(chi2,2)*q*pow(1./(1 + q) - q/(1 + q),2))/256. + 
//      (375.*pow(chi2l,2)*q*pow(1./(1 + q) - q/(1 + q),2))/256. + 
//      (35.*pow(1/(1 + q) - q/(1 + q),2))/(256.*pow(1 + q,2)) + 
//      (15.*pow(1./(1 + q) - q/(1 + q),3))/(1024.*q*(1 + q))) - (35.*M_PI*log(omega))/48. + 
//   (2995.*chi2l*q*log(omega))/9216. - (35.*pow(chi2,2)*chi2l*pow(q,3)*log(omega))/384. + 
//   (2545.*chi2l*pow(q,2)*log(omega))/(1152.*pow(1 + q,2)) + 
//   (5.*chi2l*q*(1./(1 + q) - q/(1 + q))*log(omega))/(3.*(1 + q)) + 
//   (2035.*chi2l*(1 + q)*(1/(1 + q) - q/(1 + q))*log(omega))/21504. - 
//   (5.*M_PI*(1 + q)*(1./(1 + q) - q/(1 + q))*log(omega))/(16.*q) - 
//   (5.*pow(chi2,2)*chi2l*pow(q,2)*(1 + q)*(1./(1 + q) - q/(1 + q))*log(omega))/128. + 
//   (5.*chi2l*pow(1./(1 + q) - q/(1 + q),2)*log(omega))/16.;
	alpha = (-5*(338688*pow(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q + 168*pow(chi2,2)*pow(1 + q,2)*(3 + 4*q)))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow(1 + q,4)*(3 + 4*q) - 
             84672*pow(chi2l,2)*q*pow(1 + q,2)*(3 + 4*q)*
              (75 + q*(113 + 6*pow(chi2,2)*q*pow(1 + q,2))) + 
             q*(-271003598 + 127008*pow(chi2,4)*pow(q,2)*pow(1 + q,4)*(3 + 4*q) - 
                q*(327403764 + q*(181442579 + 39432548*q)) - 
                1512*pow(chi2,2)*pow(1 + q,2)*(213 + q*(1802 + q*(2909 + 956*q)))))) + 
       1008*omega*pow(1 + q,2)*(1344*M_PI*pow(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(-5253 + q*(-18854 + q*(-18197 - 2972*q + 168*pow(chi2,2)*pow(1 + q,2)*(3 + 4*q)))))*
        log(omega)))/(6.5028096e7*omega*q*pow(1 + q,4));
	return alpha;
}


template<class T>
T IMRPhenomPv2<T>::epsilon(T omega, T q, T chi2l, T chi2)
{
	T epsilon;
//epsilon = ((-35.*chi2l*q)/128. - (15.*chi2l*(1 + q)*(1./(1 + q) - q/(1 + q)))/128.)/pow(omega,0.6666666666666666) + 
//   (-0.18229166666666666 - (5.*(1 + q)*(1./(1 + q) - q/(1 + q)))/(64.*q))/omega + 
//   (-1.7952473958333333 - (515.*q)/(384.*pow(1 + q,2)) - (175.*(1./(1 + q) - q/(1 + q)))/(256.*(1 + q)) - 
//      (4555.*(1 + q)*(1/(1 + q) - q/(1 + q)))/(7168.*q) - (15.*pow(1./(1 + q) - q/(1 + q),2))/(256.*q))/
//    pow(omega,0.3333333333333333) + pow(omega,0.3333333333333333)*
//    (4.318908476114694 - (35.*chi2l*M_PI*q)/16. + (955.*pow(q,2))/(576.*pow(1 + q,4)) + 
//      (39695.*q)/(86016.*pow(1. + q,2)) + (1645.*pow(chi2l,2)*pow(q,3))/(192.*pow(1 + q,2)) + 
//      (2725.*q*(1./(1 + q) - q/(1 + q)))/(3072.*pow(1 + q,3)) - 
//      (265.*(1./(1 + q) - q/(1 + q)))/(14336.*(1 + q)) + 
//      (1815.*pow(chi2l,2)*pow(q,2)*(1./(1 + q) - q/(1 + q)))/(256.*(1 + q)) - 
//      (15.*chi2l*M_PI*(1 + q)*(1./(1 + q) - q/(1 + q)))/16. + 
//      (27895885.*(1 + q)*(1./(1 + q) - q/(1 + q)))/(2.1676032e7*q) + 
//      (1615.*pow(1./(1 + q) - q/(1 + q),2))/(28672.*q) + 
//      (375.*pow(chi2l,2)*q*pow(1./(1 + q) - q/(1 + q),2))/256. + 
//      (35.*pow(1./(1 + q) - q/(1 + q),2))/(256.*pow(1 + q,2)) + 
//      (15.*pow(1./(1 + q) - q/(1 + q),3))/(1024.*q*(1 + q))) - (35.*M_PI*log(omega))/48. + 
//   (2995.*chi2l*q*log(omega))/9216. + (2545.*chi2l*pow(q,2)*log(omega))/(1152.*pow(1 + q,2)) + 
//   (5.*chi2l*q*(1./(1 + q) - q/(1 + q))*log(omega))/(3.*(1 + q)) + 
//   (2035.*chi2l*(1 + q)*(1./(1 + q) - q/(1 + q))*log(omega))/21504. - 
//   (5.*M_PI*(1 + q)*(1./(1 + q) - q/(1 + q))*log(omega))/(16.*q) + 
//   (5.*chi2l*pow(1./(1 + q) - q/(1 + q),2)*log(omega))/16.;
	epsilon = (-5*(338688*pow(1 + q,4)*(3 + 4*q) + 508032*chi2l*pow(omega,0.3333333333333333)*q*pow(1 + q,4)*
        (3 + 4*q) + 3024*pow(omega,0.6666666666666666)*pow(1 + q,2)*
        (2985 + q*(12890 + q*(15789 + 4988*q))) + 
       pow(omega,1.3333333333333333)*(-17660607 + 
          q*(-107348840 + 4064256*chi2l*M_PI*pow(1 + q,4)*(3 + 4*q) - 
             84672*pow(chi2l,2)*q*pow(1 + q,2)*(3 + 4*q)*(75 + 113*q) - 
             q*(271003598 + q*(327403764 + q*(181442579 + 39432548*q))))) - 
       1008*omega*pow(1 + q,2)*(-1344*M_PI*pow(1 + q,2)*(3 + 4*q) + 
          chi2l*q*(5253 + q*(18854 + q*(18197 + 2972*q))))*log(omega)))/(6.5028096e7*omega*q*pow(1 + q,4));
	return epsilon;
}

template<class T>
T IMRPhenomPv2<T>::d(int l, int m_prime, int m,T s)
{
	T sqrt2 = sqrt(2);
	/*cos(beta/2)*/
	T cb = (1./sqrt2)*sqrt(1 + 1./sqrt(1+s*s));
	/*sin(beta/2)*/
	T sb = (1./sqrt2) * sqrt( 1. - 1./sqrt(1. + s*s ));
	
	T overall_factor = sqrt(factorial(l+m) * factorial(l-m) * factorial(l +m_prime) *factorial(l-m_prime));
	
	/*Limits of the sum (determined by keeping all factorials positive)*/
	int kmax;	
	int kmax_options[3];
	if (l == std::abs(m))
		kmax = 0;
	else
	{
		if(m < m_prime) return 0;
		else
		{
			kmax_options[0] = m-m_prime;	
			kmax_options[1] = l+m;
			kmax_options[2] = l-m;
			kmax = *std::min_element(kmax_options, kmax_options+3);
		}
	}

	/*Compute the rotation matrix*/
	int k = 0;
	T sum = 0;
	while(k<=kmax)
	{
		sum=sum+ pow(-1.,k+m_prime-m)/ ( factorial(l+m-k)*factorial(l-m-k)*
					factorial(k)*factorial(m_prime-m+k) ) *
					pow(cb,2*l-2*k-m_prime+m)*pow(sb,2*k+m_prime-m);
		k++;
	}
	return sum;
	
	
}


template class IMRPhenomPv2<double>;
template class IMRPhenomPv2<adouble>;
