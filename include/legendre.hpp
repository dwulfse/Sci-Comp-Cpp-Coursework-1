#ifndef legendre_header
#define legendre_header

//////////////////////////////////////////////////////////////
//Module that computes the Legendre polynomials and finds
//their roots
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

double ComputeLegendre(double x, int p);
double computeLegendreDerivative(double x, int p);
double findRootNewtonRaphson(double (*f)(double, int), double (*df)(double, int), double x0, double tol, int p);

#endif