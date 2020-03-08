#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <cstdint>
#include <vector>
#include <fstream>
#include <functional>

double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h);
void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, std::ofstream& file);
void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, const std::string& fname);

#endif //__FUNCTIONS_H__