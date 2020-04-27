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
	uint32_t Nx=200; //200,1000
	uint32_t Nt=2000; //2000,41000
	// uint32_t Nx=1000; //200,1000
	// uint32_t Nt=41000; //2000,41000

	double r=0.1;
	double sigma=0.4;
	double E=10;
	double T=0.25;
	double T0=0;
	double xf=3*E;

	double dt_max=1/(sigma*sigma*(Nx-1)+r/2);
	double dt=(T-T0)/(Nt-1);
	double dx=xf/(Nx-1);
	double a=2*r/sigma/sigma-1;
	double b=-2*r/sigma/sigma;
	double alpha=sigma*sigma*dt;
	double beta=r*dt;

	std::cout<<
		"dt_max: "<<dt_max<<'\n'<<
		"dt: "<<dt<<'\n'<<
		"dx: "<<dx<<'\n'<<
		"alpha: "<<alpha<<'\n'<<
		"beta: "<<beta<<'\n'<<
		std::endl;

	// x=S
	std::vector<double> x(Nx);
	linspace(x,0,xf);

	// u=V
	std::vector<double> u_x1(Nx);
	std::vector<double> u_x2(Nx);
	std::vector<double> *u_x_p0=&u_x1;
	std::vector<double> *u_x_p1=&u_x2;
	std::vector<double> *tmp;

	// CI: u(0,x)
	auto g_i=[&] (uint32_t i)->double {
		double res;
		res=x[i]-E;
		// gently re-add E at boundary by exponentiation
		// res+=E*std::exp((x[i]-xf)/(xf/4));
		return (res>0)?res:0;
	};
	map(*u_x_p0,g_i);
	// CI: u(t,0)=0
	(*u_x_p0)[0]=0;
	(*u_x_p1)[0]=0;
	// CI: u(t,L)=2*E+E*(1-exp(-t/tau))
	auto q_i=[&] (uint32_t m)->double {
		// return 2*E+E*(1-std::exp(-(T0+m*dt)/((T-T0)/5)));
		return 2*E;
	};
	(*u_x_p0)[Nx-1]=q_i(0);
	(*u_x_p1)[Nx-1]=q_i(0);

	#ifdef EXPLICIT
	// Setup output stream
	std::ofstream file("output_data/res.dat");
	file<<"0";
	for (uint32_t n=0;n<Nx;++n){
		file<<','<<x[n];
	}
	file<<'\n';

	// Write initial condition
	file<<T0;
	for (uint32_t n=0;n<Nx;++n){
		file<<','<<(*u_x_p0)[n];
	}
	file<<'\n';

	// Evolution
	for (uint32_t m=1;m<Nt;++m){
		// Bound condition
		(*u_x_p1)[Nx-1]=q_i(m);

		// Write time and first bound
		file<<m*dt+T0;
		file<<','<<(*u_x_p1)[0];

		for (uint32_t n=1;n<Nx-1;++n){
			(*u_x_p1)[n]=
				(*u_x_p0)[n-1]*0.5*(n*n*alpha-n*beta)+
				(*u_x_p0)[n]*(1-n*n*alpha-beta)+
				(*u_x_p0)[n+1]*0.5*(n*n*alpha+n*beta);

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
	#endif

	// Full implicit method: LU decomposition
	// A * v_{m+1} = b_m

	auto d_n = [&] (uint32_t n)->double{
		return 1+beta+alpha*n*n;
	};
	auto u_n = [&] (uint32_t n)->double{
		return -1.0/2*(beta*(n-1)+alpha*(n-1)*(n-1));
	};
	auto l_n = [&] (uint32_t n)->double{
		return 1.0/2*(beta*(n+1)+alpha*(n+1)*(n+1));
	};

	std::vector<double> aa(Nx-1);
	std::vector<double> bb(Nx-1);
	std::vector<double> cc(Nx-1);

	std::vector<double> xx(Nx-1);

	// eq. 57
	aa[0]=d_n(0+0+1);
	cc[0]=u_n(0+1+1);
	for(uint32_t n=1;n<Nx-1;++n){
		aa[n]=d_n(n+1)-l_n(n)/aa[n-1]*u_n(n-1+2);
		bb[n]=l_n(n)/aa[n-1];
		cc[n]=u_n(n+2);
	}

	// Setup output stream
	std::ofstream file("output_data/res_LU.dat");
	file<<"0";
	for (uint32_t n=0;n<Nx;++n){
		file<<','<<x[n];
	}
	file<<'\n';

	// Write initial condition
	file<<T0;
	for (uint32_t n=0;n<Nx;++n){
		file<<','<<(*u_x_p0)[n];
	}
	file<<'\n';

	// Evolution
	for(uint32_t m=1;m<Nt;++m){
		// Eq. 44 gina dura
		(*u_x_p0)[1]-=l_n(0)*0;
		(*u_x_p0)[Nx-2]-=u_n(Nx-1)*q_i(m+1);

		// Eq. 61
		xx[0]=(*u_x_p0)[1];
		for(uint32_t n=1;n<Nx-2;++n){
			xx[n]=(*u_x_p0)[n+1]-bb[n]*xx[n-1];
		}

		// eq. 64
		(*u_x_p1)[Nx-2]=xx[Nx-2]/aa[Nx-2];
		for(uint32_t n=1;n<Nx-2;++n){
			(*u_x_p1)[Nx-n-2]=xx[Nx-n-2]/aa[Nx-n-2]
				-cc[Nx-n-2]/aa[Nx-n-2]*(*u_x_p1)[Nx-n-1];
		}


		// TODO: Save u_x_p1
		// Write time and first bound
		file<<m*dt+T0;
		for(uint32_t n=0;n<Nx;++n){
			file<<','<<(*u_x_p1)[n];
		}
		file<<'\n';

		// swap
		tmp=u_x_p0;
		u_x_p0=u_x_p1;
		u_x_p1=tmp;
	}
	file.close();

	return 0;
}
