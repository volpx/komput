
#include "simulation.h"
#include "differentiate.h"
#include "functions.h"
#include "vector_help.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>

// Interaction potential
double V_LJ(double x)
{
	return 4 * (std::pow(1.0 / x, 12) - std::pow(1.0 / x, 6));
}

// Derivative of the potential (can be omitted)
double dV_LJ(double x)
{
	return -4 * (12.0 / std::pow(x, 13) - 6.0 / std::pow(x, 7));
}

double Vtot(
	const std::vector<Vec3D> &pos,
	const std::function<double(double)> &Vr,
	const double L)
{
	// Number of particles
	size_t N{pos.size()};

	double d;
	Vec3D alias;

	// Cumulative value to return
	double Vreturn{0};

	for (size_t i{0}; i < N - 1; i++)
	{
		// Acceleration on i
		// caused by all other particles
		for (size_t j{i + 1}; j < N; j++)
		{
			// Copy the particle j
			alias = pos[j];
			// Check if there is interaction
			// and also compute the correct alias
			if ((d = compute_alias(pos[i], alias, L)) > 0)
			{
				// add potential
				Vreturn += Vr(d);
			}
		}
	}
	return Vreturn;
}

int main(int argc, char const *argv[])
{
	// Number of cells on side
	constexpr uint32_t n = 5;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n;
	// Number of evo steps
	constexpr uint32_t M = 100;
	// Number of volumes
	constexpr uint32_t Nv = 10;

	// Physics dimensional quantities
	constexpr double e_r = 1.6e-19;	  // C
	constexpr double kB_r = 1.38e-23; // J*K^-1
	constexpr double T_r = 119.76;	  // K

	constexpr double sigma_r = 3.405e-10; // m
	constexpr double eps_r = kB_r * T_r;  // J

	constexpr double rho_star_r = 1.2;
	constexpr double rho_r = rho_star_r / (sigma_r * sigma_r * sigma_r);

	// Problem adimensional quantities
	constexpr double Temp{kB_r / eps_r * T_r}; // Temperature
	constexpr double beta{1 / (kB_r * T_r) * eps_r};
	constexpr double rho{rho_star_r};

	constexpr double V0{N / rho};

	// Print info data
	std::cout
		<< "\nM:\t" << M
		<< "\nN:\t" << N
		// << "\nKrho:\t" << rho_r / rho
		// << "\nKl:\t" << dl_r / dl
		// << "\nKm:\t" << m_r / m
		<< "\nKe:\t" << eps_r / 1
		<< "\nT:\t" << Temp
		<< "\nM:\t" << M
		// << "\nV_LJ(L/2):\t" << V_LJ(L / 2)
		<< std::endl;

	double Delta;
	double L;
	double A;
	double nu;
	double Prob;
	double sum, sum2;
	double Vp0, Vp1;
	double d;
	double F;
	double P;
	bool cond;
	Vec3D alias;

	std::vector<Vec3D> pos_a(N);
	std::vector<Vec3D> pos_b(N);

	std::vector<Vec3D> *pos0 = &pos_a;
	std::vector<Vec3D> *pos1 = &pos_b;
	std::vector<Vec3D> *tmp;

	std::vector<double> V(Nv);
	// linspace(V, 0.01, 1.2);
	// V[0] = V0;
	auto fill_V_i = [&V0, &Nv](uint32_t i) -> double {
		return V0 * std::pow(1.5, (static_cast<double>(i) - Nv / 2) * 2 / Nv);
	};
	map(V, fill_V_i);

	std::ofstream PV_file("output_data/PV.dat");
	PV_file << "#P V\n";

	std::cout << "Start" << std::endl;
	for (uint32_t v = 0; v < V.size(); v++)
	{
		L = std::pow(V[v], 1.0 / 3);
		Delta = L / n / 200;

		// Initial uniform position
		init_lattice(*pos0, L, n, 1);
		apply_periodic_bounds(*pos0, L);

		Vp0 = Vtot(*pos0, V_LJ, L);

		sum2 = sum = 0;
		for (uint32_t i = 0; i < M; i++)
		{
			// New state
			for (uint32_t j = 0; j < N; j++)
			{
				(*pos1)[j].x = (*pos0)[j].x + Delta * (randu() - 0.5);
				(*pos1)[j].y = (*pos0)[j].y + Delta * (randu() - 0.5);
				(*pos1)[j].z = (*pos0)[j].z + Delta * (randu() - 0.5);
			}
			apply_periodic_bounds(*pos1, L);

			Vp1 = Vtot(*pos1, V_LJ, L);
			A = std::exp(-beta * (Vp1 - Vp0));
			Prob = std::min(1.0, A);

			// Decide if accept the new values
			nu = randu();
			std::cout << "Probability:\t" << Prob << "\tNu:\t" << nu << "\n";
			if (nu < Prob)
			{
				// accept the new configuration
				tmp = pos0;
				pos0 = pos1;
				pos1 = tmp;

				Vp0 = Vp1;
			}
			// else don't accept new configuration

			// Pos0 is now the position to make stats

			// Calculate the quantity wanted
			for (uint32_t ii{0}; ii < N; ii++)
			{
				for (uint32_t kk{0}; kk < N; kk++)
				{
					if (kk != ii)
					{
						// Copy the particle j
						alias = (*pos0)[kk];
						// Check if there is interaction
						// and also compute the correct alias
						if ((d = compute_alias((*pos0)[ii], alias, L)) > 0)
						{
							// add potential
							F = 0.5 * dV_LJ(d) * d;
							sum += F;
							sum2 += F * F;
						}
					}
				}
			}
		}
		//std::cout << "Progress:" << v * 100 / V.size() << std::endl;
		sum /= M;
		sum2 /= M;
		P = N / V[v] - 1.0 / 3 / V[v] * sum;

		PV_file << P << " " << V[v] << "\n";
	}

	PV_file.close();
	std::cout << "End" << std::endl;
	return 0;
}
