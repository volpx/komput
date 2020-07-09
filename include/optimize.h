#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include <cstdint>
#include <iostream>
#include <vector>

extern "C"
{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
}

int fit_to_exp(double lambda[3],
			   std::vector<double> &datax,
			   std::vector<double> &datay);

#endif // __OPTIMIZE_H__
