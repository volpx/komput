#include "functions.h"

#include <cstdint>
#include <iostream>
// #include <fstream>
// #include <vector>
// #include <cmath>
// #include <boost/format.hpp>

double F(double x){
	return (x-0.2132123)*(x-1.1)*(x+1.1)/(x-0.5)/(x+0.5);
}
int main(){
	std::cout<<findzero_secants_xdigits(F,0,0.51,2,-10,10)<<std::endl;
	return 0;
}
