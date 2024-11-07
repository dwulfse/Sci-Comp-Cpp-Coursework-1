//////////////////////////////////////////////////////////////
//Module that implements a double precision vector and a 
//number of vector operations
//////////////////////////////////////////////////////////////

#include "legendre.hpp"

//////////////////////////////////////////////////////////////
double ComputeLegendre(double x, int p)
//Computes the value of the pth order Legendre polynomial at
//a given point x
{
	if (p==0)
	{
		return 1;
	}
	else if (p==1)
	{
		return x;
	}
	else
	{
        double val2 = 1;
        double val1 = x;
        double val;
        for (int k=2;k<=p;k++)
        {
            val = ((2.0*double(k)-1.0)*x*val1-double(k-1)*val2)/double(k);

            val2 = val1;
            val1 = val;
        }

        return val;
	}
}

double computeLegendreDerivative(double x, int p)
// Computes the value of the derivative of the pth order Legendre
// polynomial at a given point x
{
	// (1-x^2)P'_p(x) = p(xP_p(x)-P_{p-1}(x))
	if (p==0)
	{
		return 0;
	}
	else if (p==1)
	{
		return 1;
	}
	else
	{
		return p*(x*ComputeLegendre(x,p)-ComputeLegendre(x,p-1))/(1-x*x);
	}
}

double findRootNewtonRaphson(double (*f)(double, int), double (*df)(double, int), double x0, double tol, int p)
{
	double x_old = x0;
	double x = x0;

	do 
	{
		x_old = x;
		x = x_old - f(x_old,p)/df(x_old,p);
	} while (abs(x - x_old) > tol);

	return x;
}
