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


int main(int argc, char const *argv[]){
	uint32_t Nx=200; //1000
	uint32_t Nt=2000; //41000

	double r=0.1;
	double sigma=0.4;
	double E=10;
	double T=0.25;

	double dt_max=1/(sigma*sigma*(Nx-1)+0.5*r);
	double dt=T/Nt;
	double dx=3*E/(Nx-1);
	double alpha=sigma*sigma*dt;
	double beta=r*dt;

	std::cout<<
		"dt_max: "<<dt_max<<'\n'<<
		"dt: "<<dt<<'\n'<<
		"alpha: "<<alpha<<'\n'<<
		"beta: "<<beta<<'\n'<<std::endl;

	std::vector<double> x(Nx);
	linspace(x,0,3*E);

	std::vector<double> u_x1(Nx);
	std::vector<double> u_x2(Nx);
	std::vector<double> *u_x_p0=&u_x1;
	std::vector<double> *u_x_p1=&u_x2;
	std::vector<double> *tmp;

	// CI: u(0,x)
	auto g_i=[&] (uint32_t i)->double {
		double res;
		return ((res=2*std::exp(r/sigma/sigma)*std::sinh(x[i]/2))>0)?res:0;
	};
	map(*u_x_p0,g_i);
	// CI: u(t,0)=0
	(*u_x_p0)[0]=0;
	(*u_x_p1)[0]=0;
	// CI: u(t,L)=3*E
	(*u_x_p0)[Nx-1]=3*E;
	(*u_x_p1)[Nx-1]=3*E;

	std::ofstream file("output_data/res.dat");
	file<<"0";
	for (uint32_t n=0;n<Nx;++n){
		file<<','<<x[n];
	}
	file<<'\n';

	// Evolution
	for (uint32_t m=1;m<Nt;++m){
		file<<m*dt;
		file<<','<<(*u_x_p1)[0];
		for (uint32_t n=1;n<Nx-1;++n){
			// (*u_x_p1)[n]=
			// 	(*u_x_p0)[n+1]*(n*n*alpha/2+n*beta)+
			// 	(*u_x_p0)[n]*(1-n*n*alpha-beta*(1+n))+
			// 	(*u_x_p0)[n-1]*(n*n*alpha/2) ;
			(*u_x_p1)[n]=
				(*u_x_p0)[n]+
				dt/(dx*dx)*((*u_x_p0)[n+1]-2*(*u_x_p0)[n]+(*u_x_p0)[n-1]);
			file<<','<<(*u_x_p1)[n];
		}
		file<<','<<(*u_x_p1)[Nx-1];
		file<<'\n';

		// swap
		tmp=u_x_p0;
		u_x_p0=u_x_p1;
		u_x_p1=tmp;
	}
	file.close();

	return 0;
}

/*
u=u(t,x)
CI:
u(0,x)=g(x)
u(t,0)=q(t)
u(t,L)=q(t)

// Finale
u(T,x)=0
*/
