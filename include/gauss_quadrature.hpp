#ifndef gauss_quadratue_header
#define gauss_quadratue_header

#include <vector>

struct GaussQuadrature {
	std::vector<double> x, w;
};

GaussQuadrature generateGaussQuadrature(int n);

#endif