#include <iostream>
#include <cmath>

#include "gauss_quadrature.hpp"

int main()
{
	const double pi = 3.14159265358979323846;
	const double e = 2.71828182845904523536;
	int n;

	std::cout << "Enter a number of Gauss Quadrature points: ";
	std::cin >> n;

	GaussQuadrature gauss_quadrature = generateGaussQuadrature(n);

	auto integrand_1 = [e](double x) { return pow(e, x); };
	auto integrand_2 = [pi](double x) { return cos(pi * x / 2); };

	double integral_1 = 0.0;
	double integral_2 = 0.0;

	for (int i = 0; i < n; i++)
	{
		integral_1 += gauss_quadrature.w[i] * integrand_1(gauss_quadrature.x[i]);
		integral_2 += gauss_quadrature.w[i] * integrand_2(gauss_quadrature.x[i]);
	}

	std::cout << "Integral of e^x from -1 to 1: " << integral_1 << std::endl;
	std::cout << "Integral of cos(pi*x/2) from -1 to 1: " << integral_2 << std::endl;

	return 0.0;
}