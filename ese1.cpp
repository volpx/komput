#include "functions.h"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


double V(double x){
	// adimensional version
	return 4*(std::pow(1.0/x,12)-std::pow(1.0/x,6));
}

double F(double E_tilde,double gamma, uint8_t n){
	// adimensional version
	// const double d{1e-3};
	// const double eps{1e-6};
	// const double h{1e-8};
	const uint32_t M{(uint32_t)1e5};

	//double xmin { std::pow(2,1.0/6) };
	// la funzione era invertibile
	// double xin { 
	// 	findzero_newton_raphson([&] (double x) -> double { return V(x)-E_tilde; } , xmin-d, eps, h) 
	// };
	// double xout { 
	// 	findzero_newton_raphson([&] (double x) -> double { return V(x)-E_tilde; } , xmin+d, eps, h) 
	// };
	double xin{ std::pow( -2/E_tilde + 2*std::sqrt( (1/E_tilde+1)/E_tilde) ,1.0/6) };
	double xout{ std::pow( -2/E_tilde - 2*std::sqrt( (1/E_tilde+1)/E_tilde) ,1.0/6) };
	
	#ifdef DEBUG
	std::cout<<xmin<<" "<<xin<<" "<<xout<<std::endl;
	#endif

	double I{
		integrator_simpson_cubic([&] (double x) -> double { return std::sqrt(std::abs(E_tilde-V(x))) ; }, xin, xout, M)
	};

	return gamma*I-2*M_PI*(n+0.5);
}


int main(){
	// #ifdef DEBUG
	// 	std::cout<<std::pow(1/1.122,12)<<std::endl;
	// #else

	double gamma{21.7};
	double autov_E_tilde;

	for (int n{1};n<=4;++n){
		double E_tilde_0{-0.5};
		double epsilon{1e-6};
		double h_diff{1e-8};

		auto G{
			[&] (double E_tilde)->double { return F(E_tilde,gamma,n); }
		};
		autov_E_tilde=findzero_newton_raphson(G,E_tilde_0,epsilon,h_diff);
		std::cout<<"Autovalore n="<<n<<", E_tilde="<<autov_E_tilde<<std::endl;
		
	}
	return 0;
	// #endif
}

double V_d(double r){
	const double sigma{1},epsilon{1};
	return 4*epsilon*(std::pow(sigma/r,12)-std::pow(sigma/r,6));
}

double S_d(double E){
	const double sigma{1};
	const double m={1};
	const uint32_t M{(uint32_t)1e3};

	double rmin { std::pow(2,1.0/6)*sigma };
	double rin { 
		findzero_newton_raphson([&] (double r) { return V(r)-E; } , rmin-0.0001, 0.0001, 0.0001) 
	};
	double rout { 
		findzero_newton_raphson([&] (double r) { return V(r)-E; } , rmin+0.0001, 0.0001, 0.0001) 
	};

	return integrator_simpson_cubic( [&] (double r) { return std::abs(std::sqrt(2*m*(E-V(r)))); } ,
		rin,rout,M); 
}

double F_d(double E){
	const uint16_t n{1};

	return S_d(E) - 2*M_PI*1*(n+1.0/2);
}