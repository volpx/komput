#include "functions.h"

int number_of_significant_digits(double x,double corr){
	return (int) std::log10(x/corr);
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, std::ofstream& file, const int digits){
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
			file<< boost::format((boost::format("%%.%de")%digits).str()) % data[j].second[i];
			if(j != cols - 1)
				file << ","; // No comma at end of line
		}
		file<<'\n';
	}

	return;
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>>& data, const std::string& fname, const int digits){
	std::ofstream file(fname);
	tocsv(data,file,digits);
	file.close();
}