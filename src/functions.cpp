#include "functions.h"
#include <iostream>

double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h){
	double k1{ F(t,y) };
	double y_pred{ y+h*k1 };
	double k2{ F(t+h,y_pred) };
	return y + h*0.5*(k1+k2);
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, std::ofstream& file){
	uint64_t cols=data.size();

	file << "#";
	// Send column names to the stream
    for(uint64_t j = 0; j < cols; ++j)
    {
        file << data.at(j).first;
        if(j != cols - 1)
			file << ","; // No comma at end of line
    }
    file << "\n";

	for(uint64_t i{0};i<data[0].second.size(); ++i){
		for(uint64_t j{0};j<cols; ++j){
			file<<data[j].second[i];
			if(j != cols - 1)
				file << ","; // No comma at end of line
		}
		file<<'\n';
	}

	return;
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, const std::string& fname){
	std::ofstream file(fname);
	tocsv(data,file);
	file.close();
}