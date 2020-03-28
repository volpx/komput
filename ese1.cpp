#include "integrate.h"
#include "findzero.h"
#include "functions.h"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>


double V(double x){
	// adimensional version
	return 4*(std::pow(1.0/x,12)-std::pow(1.0/x,6));
}

double F(double E_tilde,double gamma, uint8_t n){
	// adimensional version
	// const double d{1e-3};
	// const double eps{1e-6};
	// const double h{1e-8};
	const uint32_t M{(uint32_t)1e6};

	//double xmin { std::pow(2,1.0/6) };
	// la funzione era invertibile
	// double xin { 
	// 	findzero_newton_raphson([&] (double x) -> double { return V(x)-E_tilde; } , xmin-d, eps, h) 
	// };
	// double xout { 
	// 	findzero_newton_raphson([&] (double x) -> double { return V(x)-E_tilde; } , xmin+d, eps, h) 
	// };
	double xin{ std::pow( -2/E_tilde - 2*std::sqrt( (1/E_tilde+1)/E_tilde) ,1.0/6) };
	double xout{ std::pow( -2/E_tilde + 2*std::sqrt( (1/E_tilde+1)/E_tilde) ,1.0/6) };

	double I{
		(E_tilde>=0) ? 0.8413092631952725567050114474301765:
			integrator_simpson_cubic([&] (double x) -> double { return std::sqrt(std::abs(E_tilde-V(x))) ; }, xin, xout, M)
	};

	return gamma*I-M_PI*(n+0.5);
}

int main(){
	
	// double x{},corr{};
	// std::cin>>x;
	// std::cin>>corr;
	// std::cout<<"Sig dig: "<<number_of_significant_digits(x,corr)<<std::endl;
	// return 0;

	const double gamma{150};
	const int n_max{39};

	std::vector<double> aut(n_max+1);

	std::cout<<"Gamma= "<<gamma<<std::endl;
	
	double E_tilde_0{-0.5};
	// double epsilon{1e-10};
	// double h_diff{1e-10};
	// double hw{1};
	double autov_E_tilde{E_tilde_0};
	double autov_E_tilde_harmonic{};

	// std::cout<<"E_tilde_0= "<<E_tilde_0<<std::endl;
	std::cout<<"n_max= "<<n_max<<std::endl;
	// std::cout<<"h_diff= "<<h_diff<<std::endl;


	// double E;
	// for (int i{0};i<101;++i){
	// 	E=-1+1.0/100*i;
	// 	std::cout<<"E= "<<E<<",F(E)= "<<F(E,150,25)<<std::endl;
	// }

	// int n{0};
	// while(autov_E_tilde<0){
	for (int n{0};n<=n_max;++n){

		std::function<double(double)> G{
			[&] (double E_tilde)->double { return F(E_tilde,gamma,n); }
		};
		// autov_E_tilde=findzero_newton_raphson_xeps(G,autov_E_tilde,epsilon,h_diff,-INFINITY,0-epsilon);
		// autov_E_tilde=findzero_secants_xeps(G,autov_E_tilde,0,1e-10);
		autov_E_tilde=findzero_secants_xdigits(G,autov_E_tilde,0,7,-1,0);
		// autov_E_tilde=findzero_bisection_xeps(G,autov_E_tilde,0-epsilon,epsilon);

		autov_E_tilde_harmonic=1/gamma*std::sqrt(24*13/std::pow(2,7.0/3)-12*7/std::pow(2,4.0/3))*(n+0.5)-1;

		std::cout<<"Autovalore n="<<n<<boost::format(", E_tilde= %.7e, E_tilde_h= %.7e") % autov_E_tilde % autov_E_tilde_harmonic <<std::endl;

		aut[n]=autov_E_tilde;
		// ++n;
	}


	tocsv({{"autov",aut}},"output_data/autv_150.csv",7);
	return 0;
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
		findzero_newton_raphson_xeps([&] (double r) { return V(r)-E; } , rmin-0.0001, 0.0001, 0.0001) 
	};
	double rout { 
		findzero_newton_raphson_xeps([&] (double r) { return V(r)-E; } , rmin+0.0001, 0.0001, 0.0001) 
	};

	return integrator_simpson_cubic( [&] (double r) { return std::abs(std::sqrt(2*m*(E-V(r)))); } ,
		rin,rout,M); 
}

double F_d(double E){
	const uint16_t n{1};

	return S_d(E) - 2*M_PI*1*(n+1.0/2);
}