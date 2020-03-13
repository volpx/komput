#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <cstdint>
#include <vector>
#include <fstream>
#include <functional>
#include <cmath>

/* Plain implementation of the predictor corrector integrator
 * 
 */
double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h);

/* Find zero with Newton Raphson
 * F is the function that this method wants to find the zero of
 * x0 is the starting value for the search
 * epsilon precision i want for the found zero
 * h is the mesh size for the derivation
 */
double findzero_newton_raphson_xeps(std::function<double(double)> F,double x0,double epsilon,double h,
                                    const double xmin=-INFINITY,const double xmax=INFINITY);
double findzero_newton_raphson_yeps(std::function<double(double)> F,double x0,double epsilon,double h);

/* Find zero with bisection
 * F function for the zero
 * xmin and xmax are the boundary
 * epsilon is the accuracy of the zero
 * NOTE: must change sign between xmin and xmax otherwise returns the middle
 */
double findzero_bisection_xeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon);
double findzero_bisection_yeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon);

/* Derivator with 5 points rule
 * of function F
 * around point x0
 * mesh size of h
 */
double derive_5points(std::function<double(double)> F,const double x0,const double h);

/* Integrate with simpson cubic rule
 * the function F
 * from xmin to xmax
 * with M steps
 */
double integrator_simpson_cubic(std::function<double(double)> F,const double xmin,const double xmax,const uint32_t M);

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, std::ofstream& file);
void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, const std::string& fname);

#endif //__FUNCTIONS_H__