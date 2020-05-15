#include "simulation.h"
#include "differentiate.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define SAVE_POS

double V_LJ(double x);
void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f);
// Compute the acceleration on particle i
double compute_acc_return_V(
	std::vector<Vec3D> &acc,
	const std::function<double(double)> &Vr,
	const std::vector<Vec3D> &pos,
	const double m,
	const double L);
int tmp_main(int argc, char const *argv[])
{
	// check acceleration
	uint32_t N = 2;
	std::vector<Vec3D> pos(N);

	pos[0] = Vec3D{0.0857488, 6.95682, 4.79487};
	pos[1] = Vec3D{12.7761, 6.95682, 4.79487};

	double sigma_r = 3.94e-10;				 // m
	double rho_r = 1e27;					 // #*m^-3
	double dl_r = std::pow(rho_r, -1.0 / 3); // m
	double dl{dl_r / sigma_r};				 // Cell length
	double L{5 * dl};

	apply_periodic_bounds(pos, L);

	std::cout << pos[0] << "\n"
			  << pos[1] << std::endl;

	return 0;
}

int main(int argc, char const *argv[])
{
	uint32_t n = 5; // number of celles on side
	uint32_t N = n * n * n;
	uint32_t M = 100000;

	// Physics quantities
	double e_r = 1.6e-19;						 // C
	double kB_r = 1.38e-23;						 // J*K^-1
	double T_r = 300;							 // K
	double m_r = 2.2e-25;						 // Kg
	double rho_r = 1e27;						 // #*m^-3
	double dl_r = std::pow(rho_r, -1.0 / 3);	 // m
	double dt_r = 1e-12;						 // s
	double vstd_r = std::sqrt(kB_r * T_r / m_r); // m*s^-1
	double sigma_r = 3.94e-10;					 // m
	double eps_r = 0.02 * e_r;					 // J

	// Modifications
	// T_r = 10;
	dt_r /= 100;
	// M = 50000;
	// rho_r *= 0.5 / 0.061163;
	// m_r *= 1 / 8.17488;
	// dt_r *= 0.01 / 0.875202;
	// dl_r *= 6.29960525 / 12.6904;
	// vstd_r *= 0;

	// Problem quantities
	double dl{dl_r / sigma_r};									  // Cell length
	double L{n * dl};											  // Box length
	double dt{std::sqrt(eps_r / m_r / sigma_r / sigma_r) * dt_r}; // Time step
	double rho{std::pow(sigma_r, 3) * rho_r};					  // Density
	double m{rho * std::pow(dl, 3)};							  // Mass of the particles
	double vstd{vstd_r * std::sqrt(m_r / eps_r)};				  // std velocity
	double Temp{kB_r / eps_r * T_r};							  // Temperature

	// Print info data
	std::cout << "\nM:\t" << M << "\nN:\t" << N << "\nL:\t" << L
			  << "\ndl:\t" << dl
			  << "\nrho:\t" << rho << "\nm:\t" << m << "\nT:\t" << Temp
			  << "\ndt:\t" << dt << "\nM:\t" << M << "\nTend:\t" << dt * M
			  << "\nV_LJ(L/2):\t" << V_LJ(L / 2) << "\nvstd:\t" << vstd
			  << std::endl;

	// Evolution variables
	std::vector<Vec3D> pos_a(N);
	std::vector<Vec3D> vel_a(N);
	std::vector<Vec3D> acc_a(N);
	std::vector<Vec3D> pos_b(N);
	std::vector<Vec3D> vel_b(N);
	std::vector<Vec3D> acc_b(N);

	std::vector<Vec3D> *pos0 = &pos_a;
	std::vector<Vec3D> *pos1 = &pos_b;
	std::vector<Vec3D> *vel0 = &vel_a;
	std::vector<Vec3D> *vel1 = &vel_b;
	std::vector<Vec3D> *acc0 = &acc_a;
	std::vector<Vec3D> *acc1 = &acc_b;
	std::vector<Vec3D> *tmp;

	double E_tot{0};
	double T{0};
	double V{0};

#ifdef SAVE_POS
	std::ofstream pos_file("output_data/positions.dat");
	pos_file << "#n x y z\n";
#endif

	std::ofstream energy_file("output_data/energy.dat");
	energy_file << "#m t T V E\n";

	// Initialize with uniform lattice conditions
	printf("Init lattice\n");
	init_lattice(*pos0, L, n, 1);
	apply_periodic_bounds(*pos0, L);

#ifdef SAVE_POS
	save_positions(*pos0, pos_file);
	pos_file << "\n\n";
#endif

	V = compute_acc_return_V(*acc0, V_LJ, *pos0, m, L);

	printf("V0:\t%f\n", V);

	// Initialize velocities as gaussian on the components
	printf("Init velocities\n");
	T = 0;
	for (uint32_t i = 0; i < N; i++)
	{
		init_distribute_maxwell_boltzmann((*vel0)[i], vstd);
		T += 0.5 * m * (*vel0)[i].norm2();
	}
	E_tot = T + V;
	printf("Initial energy: %f\n", E_tot);

	energy_file << 0 << " " << 0 << " " << T << " " << V << " " << E_tot << "\n";

	// Velocity verlet
	printf("Start velocity verlet\n");
	for (uint32_t i = 1; i < M; i++)
	{
		// Velocity - Verlet
		// Compute all new positions
		for (uint32_t j = 0; j < N; j++)
		{
			(*pos1)[j] = (*pos0)[j] + dt * (*vel0)[j] + 0.5 * dt * dt * (*acc0)[j];
		}
		// Periodic conditions
		apply_periodic_bounds(*pos1, L);

		// Compute potential and new accelerations
		V = compute_acc_return_V(*acc1, V_LJ, *pos1, m, L);

		// Compute all new velocities
		T = 0;
		for (uint32_t j = 0; j < N; j++)
		{
			(*vel1)[j] = (*vel0)[j] + 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);
			T += 0.5 * m * (*vel1)[j].norm2();
		}
		E_tot = T + V;

		// Output
#ifdef SAVE_POS
		save_positions(*pos0, pos_file);
		pos_file << "\n\n";
#endif
		energy_file << i << " " << i * dt << " " << T << " " << V << " " << E_tot << "\n";
		if ((i % 1000) == 0)
		{
			printf("Step: %7d \t Percentage: %3d \tEnergy: %f\n", i, (i * 100) / M, E_tot);
		}

		// Swap pointers
		tmp = pos1;
		pos1 = pos0;
		pos0 = tmp;

		tmp = vel1;
		vel1 = vel0;
		vel0 = tmp;

		tmp = acc1;
		acc1 = acc0;
		acc0 = tmp;
	}

	printf("Final energy: %f\n", E_tot);
	printf("End\n");
	// Pause
	// getchar();

	return 0;
}

double V_LJ(double x)
{
	return 4 * (std::pow(1.0 / x, 12) - std::pow(1.0 / x, 6));
}

// Compute the acceleration on particle i
double compute_acc_return_V(
	std::vector<Vec3D> &acc,
	const std::function<double(double)> &Vr,
	const std::vector<Vec3D> &pos,
	const double m,
	const double L)
{
	size_t N{pos.size()};
	double d;
	Vec3D tmp_acc;
	double Vreturn{0};
	Vec3D alias;

	for (size_t i{0}; i < N; i++)
	{
		acc[i].clear();
		// Acceleration of i
		for (size_t j{0}; j < N; j++)
		{
			// Skip the self case
			if (i != j)
			{
				// Iterate on all possible alias of j
				for (uint8_t a{0}; a < 27; a++)
				{
					alias = pos[j] + L * aliaser[a];
					d = (pos[i] - alias).norm();
					// Check interaction
					if (d < L / 2)
					{
						// Compute the derivative on function v around 0
						// Compute the acc contribute for the pair ij
						// add the contribute to the force using the difference versor
						acc[i] += -derive_5points(Vr, d, 1e-10) / m / d * (pos[i] - alias);
						// add potential on half of the cases
						if (i < j)
							Vreturn += Vr(d);
					}
					// else ignore interaction
				}
			}
		}
	}
	return Vreturn;
}

void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f)
{
	size_t N{pos.size()};
	for (size_t i{0}; i < N; i++)
	{
		f << i << " " << pos[i].x << " " << pos[i].y << " " << pos[i].z << "\n";
	}
}

/* junk functions ****************************************************/
// double V_LJ_diff(double x)
// {
// 	return -4 * (12 / x * std::pow(1.0 / x, 12) - 6 / x * std::pow(1.0 / x, 6));
// }
