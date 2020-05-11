#include "simulation.h"
#include "differentiate.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

double V_LJ(const double x)
{
	return 4 * (std::pow(1.0 / x, 12) - std::pow(1.0 / x, 6));
}
double V_LJ_1(const Vec3D &a, const Vec3D &b)
{
	double d2 = (a - b).norm2();
	return 4 * (std::pow(1 / d2, 6) - std::pow(1 / d2, 3));
}

double V_LJ_diff(double x)
{
	return -4 * (12 / x * std::pow(1.0 / x, 12) - 6 / x * std::pow(1.0 / x, 6));
}

double V(const std::vector<Vec3D> &pos)
{
	double res{0};
	size_t N{pos.size()};
	// TODO: check ranges
	for (size_t i{0}; i < N - 1; i++)
	{
		for (size_t j{i + 1}; j < N; j++)
		{
			res += V_LJ_1(pos[i], pos[j]);
		}
	}
	return res;
}

/* Calculate the vector of forces F
  V the potential between two particles
  pos the vector of positions
*/
void calculate_F(
	std::vector<Vec3D> &F,
	std::function<double(Vec3D, Vec3D)> Vij,
	std::vector<Vec3D> &pos)
{
	double der;
	size_t N = F.size();
	// Loop on every particle to calculate the force on each
	for (size_t k{0}; k < N; k++)
	{
		// Start with 0 to add single contributes
		F[k] = 0;
		// TODO: could this for be halved?
		// Loop on all other particles
		for (size_t i{0}; i < N; i++)
		{
			// Skip the self-case
			if (i != k)
			{
				// x=r_ij-r_ij0 ==> dx=dr
				// move the potential with the original position in the origin
				// this way the derivative can be computed in 0 on variable x
				// which directly varies the radial distance
				auto v = [&](double x) -> double {
					return V_LJ((pos[k] - pos[i]).norm() + x);
				};
				// Compute the derivative on function v around 0
				der = derive_5points(v, 0, 1e-10);
				// add the contribute to the force using the difference versor
				F[k] += -der * (pos[k] - pos[i]).normalize();
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	uint32_t N = 125;
	uint32_t n = 5; // number of celles on side
	uint32_t M = 1000;
	{	// Physics quantities
		// double m = 2.2 * std::pow(10.0, -22);	 //g
		// double k_B = 1.38 * std::pow(10.0, -23); //J/K
		// double T = 300;							 //K
		// double rho = std::pow(10.0, 27);
		// double K = std::sqrt(k_B * T / m);
		// double dl = std::pow(rho, -1.0 / 3);
		// double L = n * dl;
		// double v_std = std::sqrt(k_B * T / m);
	}

	std::vector<Vec3D> pos(N);
	std::vector<Vec3D> vel(N);
	std::vector<Vec3D> F(N);

	printf("Start\n");

	// Initialize with uniform lattice conditions
	printf("Init lattice\n");
	init_lattice(pos, 1, n, 1);

	// Initialize velocities as gaussian on the components
	printf("Init velocities\n");
	for (uint32_t i = 0; i < N; i++)
	{
		init_distribute_maxwell_boltzmann(vel[i], 1);
	}

	// Velocity verlet
	printf("Start verlet\n");
	for (uint32_t i = 0; i < M; i++)
	{
		// TODO: verlet and save datafiles
		apply_periodic_bounds_unitary(pos);
	}

	printf("End\n");
	// Pause
	// getchar();

	return 0;
}
