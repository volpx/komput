#include "functions.h"

#include <iostream>


double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h){
	double k1{ F(t,y) };
	double y_pred{ y+h*k1 };
	double k2{ F(t+h,y_pred) };
	return y + h*0.5*(k1+k2);
}

double findzero_newton_raphson_xeps(std::function<double(double)> F,const double x0,const double epsilon,const double h,
									const double xmin,const double xmax){
	double x{x0},corr{F(x)/derive_5points(F,x,h)};

	while( std::abs(corr) > epsilon ){

		#ifdef DEBUG
		std::cout<<"cond= "<< (std::abs(F(x)) - epsilon) <<std::endl;
		#endif

		#ifdef DEBUG
		std::cout<<"NR:x= "<<x<<std::endl;
		#endif 

		#ifdef DEBUG
		std::cout<<"NR:corr= "<<F(x)/derive_5points(F,x,h)<<std::endl;
		#endif

		corr=F(x)/derive_5points(F,x,h);
		x=x-corr;
		x=(x>=xmax)?xmax:x;
		x=(x<=xmin)?xmin:x;
	}
	return x;
}

double findzero_newton_raphson_yeps(std::function<double(double)> F,const double x0,const double epsilon,const double h){
	double x{x0};

	while( std::abs(F(x)) > epsilon ){

		#ifdef DEBUG
		std::cout<<"cond= "<< (std::abs(F(x)) - epsilon) <<std::endl;
		#endif

		#ifdef DEBUG
		std::cout<<"NR:x= "<<x<<std::endl;
		#endif 

		#ifdef DEBUG
		std::cout<<"NR:corr= "<<F(x)/derive_5points(F,x,h)<<std::endl;
		#endif


		x=x-F(x)/derive_5points(F,x,h);
	}
	return x;
}

double findzero_bisection_yeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon){
	// must change sign between min and max
	double xmiddle{(xmax+xmin)/2};
	if(F(xmin)*F(xmax) > 0)
		return xmiddle;

	while( std::abs(F(xmiddle)) > epsilon ){
		if (F(xmin)*F(xmiddle) < 0){
			// zero here
			xmax=xmiddle;
		}
		else{
			xmin=xmiddle;
		}
		xmiddle=(xmax+xmin)/2;
	}
	return xmiddle;
}

double findzero_bisection_xeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon){
	// must change sign between min and max
	double xmiddle{(xmax+xmin)/2};
	if(F(xmin)*F(xmax) > 0)
		return xmiddle;

	while( std::abs(xmax-xmin) > epsilon ){
		if (F(xmin)*F(xmiddle) < 0){
			// zero here
			xmax=xmiddle;
		}
		else{
			xmin=xmiddle;
		}
		xmiddle=(xmax+xmin)/2;
	}
	return xmiddle;
}

double findzero_secants_xeps(std::function<double(double)> F, double x0,double x1,const double epsilon){
	double x2{x1-F(x1)*(x1-x0)/(F(x1)-F(x0))}; 

	while(std::abs(x2-x1) > epsilon){
		x0=x1; 
		x1=x2;
		x2=x1-F(x1)*(x1-x0)/(F(x1)-F(x0)); 
	}
	return (x2+x1)/2;
}
/*
x0=-.5
x1=0
x2=sensato
x1=sensato
x0=0

*/

double derive_5points(std::function<double(double)> F,const double x0,const double h){
	// 5 points derivative
	return ( F(x0-2*h)-8*F(x0-h)+8*F(x0+h)-F(x0+2*h) ) /(12*h);
}

double integrator_simpson_cubic(std::function<double(double)> F,const double xmin,const double xmax,const uint32_t M){
	std::vector<double> f(M+1);
	// std::vector<double> x(M);
	double h{(xmax-xmin)/M};

	for (uint32_t i{0}; i<=M; ++i){
		f[i]=F(xmin+i*h);
	}

	double result{0};
	result+=f[0]+f[M];
	for (uint32_t i{1}; i<=M/2-1; ++i){
		result+=4*f[2*i-1]+2*f[2*i];
	}
	result+=4*f[M-1];

	return h/3*result;
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