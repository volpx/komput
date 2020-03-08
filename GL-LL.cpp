#include "functions.h"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

double F_L(double t,double N){
	const double al{0.16},Nl_inf{4724};
	return N*al*(1-N/Nl_inf);
}
double F_G(double t,double N){
	const double ag{0.16},Ng_inf{4724};
	return N*ag*std::log(Ng_inf/N);
}

int main(){
	double h{1.0/24}; // days
	int M{60*24};

	double N0{2};
	double t0{0};

	std::vector<double> Nl(M);
	std::vector<double> Ng(M);
	std::vector<double> Nl_dumb(M);
	std::vector<double> Ng_dumb(M);
	std::vector<double> t(M);

	Ng[0]=Nl[0]=Ng_dumb[0]=Nl_dumb[0]=N0;
	t[0]=t0;

	for (int i=0;i<M-1;++i){
		// y(t+h)=y[i+1]

		t[i+1]=t[i]+h;

		Ng_dumb[i+1]=Ng_dumb[i]+h*F_G(t[i],Ng_dumb[i]);
		Nl_dumb[i+1]=Nl_dumb[i]+h*F_L(t[i],Nl_dumb[i]);
		Ng[i+1]=predictor_corrector(F_G,Ng[i],t[i],h);
		Nl[i+1]=predictor_corrector(F_L,Nl[i],t[i],h);
	}

	tocsv(
	   {{"time",t},
		{"Ng_dumb",Ng_dumb},
		{"Nl_dumb",Nl_dumb},
		{"Ng",Ng},{"Nl",Nl}},
		"output_data/GL-LL.csv");

	return 0;
}