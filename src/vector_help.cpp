#include "vector_help.h"

void map(std::vector<double>& vec, const std::function<double(uint32_t)> f){
	const uint32_t M{ static_cast<uint32_t>(vec.size()) };
	for(uint32_t i{0};i<M;++i){
		vec[i]=f(i);
	}
}
void arange(std::vector<double>& vec, const double start, const double step){
	std::function<double(uint32_t)> f{
		[&] (uint32_t i)->double { return start+i*step; }
	};
	map(vec,f);
}
void linespace(std::vector<double>& vec, const double xmin, const double xmax){
	const double h{(xmax-xmin)/(vec.size()-1)};
	arange(vec,xmin,h);
}
void fill(std::vector<double>& vec,const double val){
	std::function<double(uint32_t)> f{
		[&] (uint32_t i)->double { return val; }
	};
	map(vec,f);
}