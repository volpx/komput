#ifndef __VECTOR_HELP_H__
#define __VECTOR_HELP_H__

#include <cstdint>
#include <vector>
#include <functional>

// vector manipulation
void map(std::vector<double> &vec, std::function<double(uint32_t)> f);
void arange(std::vector<double> &vec, const double start, const double step);
void linspace(std::vector<double> &vec, const double xmin, const double xmax);
void fill(std::vector<double> &vec, double val);

uint32_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr);
uint32_t ind_max(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr);

template <typename T>
T mean(const std::vector<T> &x)
{
	T r = 0;
	for (T v : x)
	{
		r += v;
	}
	return r / x.size();
}

template <typename T>
T variance(const std::vector<T> &x)
{
	T s1 = 0;
	T s2 = 0;
	for (T v : x)
	{
		s1 += v;
		s2 += v * v;
	}
	return s2 / x.size() - s1 * s1 / x.size() / x.size();
}
#endif // __VECTOR_HELP_H__
