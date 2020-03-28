#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <cstdint>
#include <cmath>
#include <vector>
#include <fstream>
#include <boost/format.hpp>


int number_of_significant_digits(double x,double corr);

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, std::ofstream& file,
    const int digits=5);
void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, const std::string& fname,
    const int digits=5);

#endif // __FUNCTIONS_H__