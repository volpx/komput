#ifndef __VECTOR_HELP_H__
#define __VECTOR_HELP_H__

#include <cstdint>
#include <vector>
#include <functional>

// vector manipulation
void map(std::vector<double>& vec,std::function<double(uint32_t)> f);
void arange(std::vector<double>& vec, const double start, const double step);
void linspace(std::vector<double>& vec, const double xmin, const double xmax);
void fill(std::vector<double>& vec,double val);

uint32_t ind_min(const std::vector<double>& vec,
				std::function<double(double)> map=nullptr);
uint32_t ind_max(const std::vector<double>& vec,
				std::function<double(double)> map=nullptr);
#endif // __VECTOR_HELP_H__
