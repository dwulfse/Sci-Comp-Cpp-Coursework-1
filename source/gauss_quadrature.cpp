#include <vector>

#include "legendre.hpp"

struct GaussQuadrature {
	std::vector<double> x, w;
};

GaussQuadrature generateGaussQuadrature(int n)
{
	GaussQuadrature gauss_quadrature;
	const double PI = 3.14159265358979323846;
	const double tol = 1e-10;
	double initial_guess;
	std::vector<double> x(n), w(n);


	for (int i=0; i<n; i++)
	{
		initial_guess = cos(PI * (i - 1.0/2) / n);
		x[i] = findRootNewtonRaphson(ComputeLegendre, computeLegendreDerivative, initial_guess, tol, n);
		w[i] = 2.0 / ((1 - x[i]*x[i]) * pow(computeLegendreDerivative(x[i], n), 2));
	}

	gauss_quadrature.x = x;
	gauss_quadrature.w = w;

	return gauss_quadrature;
}