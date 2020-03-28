#ifndef __VECTOR_HELP_H__
#define __VECTOR_HELP_H__

#include <cstdint>
#include <vector>
#include <functional>

// vector manipulation
void map(std::vector<double>& vec,std::function<double(uint32_t)> f);
void arange(std::vector<double>& vec, const double start, const double step);
void linespace(std::vector<double>& vec, const double xmin, const double xmax);
void fill(std::vector<double>& vec,double val);

#endif // __VECTOR_HELP_H__