#include "vector_help.h"
#include "integrate.h"
#include "findzero.h"
#include "ode_solvers.h"
#include "functions.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

// dtheta_t/dx = phi(x,theta,phi)
double f1(double x,double theta,double phi){
	return phi;
}

// dphi/dx = f2(x,theta,phi)
double f2_tmp(double x,double theta,double phi,
				double n){
	return -2/x*phi-std::pow(theta,n);
}


// modified
// dtheta_t/dx = phi(x,theta,m)
double f1_m(double x,double theta,double m){
	return -m/x/x;
}

// dm/dx = f2(x,theta,m)
double f2_tmp_m(double x,double theta,double m,
				double n){
	return std::pow(x,2)*std::pow(theta,n);
}

// TOV
double tov_f1(double x, double theta,double m,
	double n,double gamma){
		return -1.0/(2*(n+1)*gamma)*m/(x*x)/(1-m/x)*
			(1+gamma*theta)*(1+gamma*std::pow(x,3)*std::pow(theta,n+1)/m);
}
double tov_f2(double x, double theta,double m,
	double n){
		return x*x*std::pow(theta,n);
}

#define PART 1
int main(int argc, char const *argv[]){
	#if PART == 0
	// mathematical portion

	double theta_0{1};
	double phi_0{0};
	// start of integrating interval
	double x_0{0+1e-5};

	// end of the interval (exluded)
	double x_end{7};
	uint32_t M{1000};
	double h{(x_end-x_0)/M};

	std::vector<double> n(7);
	linespace(n,1.5,3);

	std::vector<double> theta(M);
	std::vector<double> phi(M);
	std::vector<double> x(M);
	theta[0]=theta_0;
	phi[0]=phi_0;
	arange(x,x_0,h);

	std::cout<<boost::format("Initial configuration\n\
		\r\tM: %d\n\tintegration step: %f\n\ttheta_0: %f\n\tphi_0: %f\n") %M % h % theta_0 % phi_0
		<<std::endl;

	for(uint32_t j=0;j<static_cast<uint32_t>(n.size());++j){
	// for(uint32_t j=0;j<;++j){
		std::cout<<"Integrating n: "<<n[j]<<std::endl;

		std::function<double(double,double,double)> f2{
			[&] (double x,double theta,double phi)->double { return f2_tmp(x,theta,phi,n[j]); }
		};


		for(uint32_t i=1;i<M;++i){
			RK_2_step(f1,f2,x[i],theta[i-1],phi[i-1],h,&(theta[i]),&(phi[i]));
		}

		tocsv(
			{{"x",x},
			{"theta",theta},
			{"phi",phi}},
			(boost::format("output_data/neutron_n_%.2f.csv") % n[j]).str());
	}

	n.clear();
	n.push_back(1.0);
	n.push_back(4.5);
	x_end=40;
	h=(x_end-x_0)/M;
	arange(x,x_0,h);

	for(uint32_t j=0;j<static_cast<uint32_t>(n.size());++j){
		std::cout<<"Integrating n: "<<n[j]<<std::endl;

		std::function<double(double,double,double)> f2{
			[&] (double x,double theta,double phi)->double { return f2_tmp(x,theta,phi,n[j]); }
		};


		for(uint32_t i=1;i<M;++i){
			RK_2_step(f1,f2,x[i],theta[i-1],phi[i-1],h,&(theta[i]),&(phi[i]));
		}

		tocsv(
			{{"x",x},
			{"theta",theta}},
			(boost::format("output_data/neutron_n_%.2f.csv") % n[j]).str() );
	}
	#elif PART==1
	// mathematical portion modified version

	double theta_0{1};
	double m_0{0};
	// start of integrating interval
	double x_0{0+1e-5};

	// end of the interval (exluded)
	double x_end{7};
	uint32_t M{1000};
	double h{(x_end-x_0)/M};

	std::vector<double> n(7);
	linespace(n,1.5,3);

	std::vector<double> theta(M);
	std::vector<double> m(M);
	std::vector<double> x(M);
	theta[0]=theta_0;
	m[0]=m_0;
	arange(x,x_0,h);

	std::cout<<boost::format("Initial configuration\n\
		\r\tM: %d\n\tintegration step: %f\n\ttheta_0: %f\n\tm_0: %f\n") %M % h % theta_0 % m_0
		<<std::endl;

	for(uint32_t j=0;j<static_cast<uint32_t>(n.size());++j){
	// for(uint32_t j=0;j<;++j){
		std::cout<<"Integrating n: "<<n[j]<<std::endl;

		std::function<double(double,double,double)> f2_m{
			[&] (double x,double theta,double m)->double { return f2_tmp_m(x,theta,m,n[j]); }
		};


		for(uint32_t i=1;i<M;++i){
			RK_2_step(f1_m,f2_m,x[i],theta[i-1],m[i-1],h,&(theta[i]),&(m[i]));
		}

		tocsv(
			{{"x",x},
			{"theta",theta},
			{"m",m}},
			(boost::format("output_data/neutron_mod_n_%.2f.csv") % n[j]).str());
	}

	n.clear();
	n.push_back(1.0);
	n.push_back(4.5);
	x_end=40;
	h=(x_end-x_0)/M;
	arange(x,x_0,h);

	for(uint32_t j=0;j<static_cast<uint32_t>(n.size());++j){
	// for(uint32_t j=0;j<;++j){
		std::cout<<"Integrating n: "<<n[j]<<std::endl;

		std::function<double(double,double,double)> f2_m{
			[&] (double x,double theta,double m)->double { return f2_tmp_m(x,theta,m,n[j]); }
		};


		for(uint32_t i=1;i<M;++i){
			RK_2_step(f1_m,f2_m,x[i],theta[i-1],m[i-1],h,&(theta[i]),&(m[i]));
		}

		tocsv(
			{{"x",x},
			{"theta",theta},
			{"m",m}},
			(boost::format("output_data/neutron_mod_n_%.2f.csv") % n[j]).str());
	}
	#elif PART==2
	// TOV
	double theta_0{1};
	double m_0{1e-10};
	double k3{41.00421351298523};
	double rho_c{1.288960613915068e+18};
	double c{299792458};
	double gamma;
	// start of integrating interval
	double x_0{0+1e-5};

	// end of the interval (exluded)
	double x_end{0.2};
	uint32_t M{1000};
	double h{(x_end-x_0)/M};

	// std::vector<double> n(7);
	// linespace(n,2.5,3.5);
	double n=3;

	std::vector<double> theta(M);
	std::vector<double> m(M);
	std::vector<double> x(M);
	theta[0]=theta_0;
	m[0]=m_0;
	arange(x,x_0,h);

	std::cout<<boost::format("Initial configuration\n\
		\r\tM: %d\n\tintegration step: %f\n\ttheta_0: %f\n\tm_0: %f\n") %M % h % theta_0 % m_0
		<<std::endl;

	gamma=k3*std::pow(rho_c,1/n)/c/c;
	std::cout<<"Integrating n: "<<n<<"\tgamma: "<<gamma<<std::endl;


	std::function<double(double,double,double)> f1_wrap{
		[&] (double x,double theta,double m)->double { return tov_f1(x,theta,m,n,gamma); }
	};

	std::function<double(double,double,double)> f2_wrap{
		[&] (double x,double theta,double m)->double { return tov_f2(x,theta,m,n); }
	};


	for(uint32_t i=1;i<M;++i){
		RK_2_step(f1_wrap,f2_wrap,x[i],theta[i-1],m[i-1],h,&(theta[i]),&(m[i]));
	}

	tocsv(
		{{"x",x},
		{"theta",theta},
		{"m",m}},
		(boost::format("output_data/neutron_TOV_n_%.2f.csv") % n).str());

	#elif PART==3
	// Physics for Lane-Emden
	double n=3.0/2;
	int M=100;
	// taken from the file
	// double theta_0{2.39488e-3};
	double phi_0{-2.05235e-1};
	double x_0{3.57};
	std::vector<double> mass(M); // M0
	linespace(mass,0.3,4);

	double h_bar=1.05457148e-34; // J*s
	double m_n=1.67492749804e-27; // Kg
	double k_1_5=std::pow(h_bar,2)*std::pow(3*M_PI*M_PI,2.0/3)/(5*std::pow(m_n,8.0/3)); //
	double G=6.67430e-11; // Jm/Kg^2

	std::cout<<"Setted n: "<<n<<"\nphi_0: "<<phi_0<<"\nk3/2: "<<k_1_5<<std::endl;

	std::function<double(double)> radius_from_mass_index{
		[&] (uint32_t i)->double {
			return std::pow( -4*M_PI*std::pow(5.0*k_1_5/8/M_PI/G,3.0)*std::pow(x_0,5)/(1.989e30*mass[i])*phi_0 ,1.0/3);
		}
	};

	std::vector<double> radius(M); // meters
	map(radius,radius_from_mass_index);

	tocsv(
		{{"radius",radius},
		{"mass",mass}},
		"output_data/mass_radius_1_5.csv" );

	#endif
	return 0;
}
