#include "vec3d.h"
#include "simulation.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>

// Interaction potential
double V_LJ(double x)
{
	return 4 * (std::pow(1. / x, 12) - std::pow(1. / x, 6));
}
// Derivative of the potential
double dV_LJ(double x)
{
	return -4 * (12. / std::pow(x, 13) - 6. / std::pow(x, 7));
}
double ddV_LJ(double x)
{
	return 4 * (12. * 13. / std::pow(x, 14) - 6. * 7. / std::pow(x, 8));
}

int main()
{
	// Number of cells on side
	constexpr uint32_t n = 5;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n;
	// Number of time steps
	constexpr uint32_t M = 100'000;

	// Physics constants
	constexpr double kB_r = 1.3806503e-23;	  // J
	constexpr double hbar_r = 1.05457148e-34; // J s

	// Argon data
	constexpr double sigma_r = 3.398e-10;			// m
	constexpr double T_r = 122;						// K
	constexpr double eps_r = T_r * kB_r;			// J
	constexpr double m_r = 39.95 * 1e-3 / 6.022e23; // Kg
	constexpr double dt_r = 1e-14;					// s
	constexpr double rho = 0.65;					//
	constexpr double lambda = 0.86e-3;				//

	std::cout
		<< "Scale constants:\n"
		<< "\tKl: " << sigma_r << " m\n"
		<< "\tKe: " << eps_r << " J\n"
		<< "\tKm: " << m_r << " Kg\n"
		<< "\tKrho: " << 1. / (sigma_r * sigma_r * sigma_r) << " m^-3\n"
		<< "\tKmrho: " << m_r / (sigma_r * sigma_r * sigma_r) << " Kg m^-3\n"
		<< "\tKT: " << eps_r / kB_r << " K\n"
		<< "\tKP: " << eps_r / (sigma_r * sigma_r * sigma_r) << " Pa\n"
		<< "\tKt: " << std::sqrt(m_r * sigma_r * sigma_r / eps_r) << " s\n"
		<< "\tKh: " << std::sqrt(eps_r * m_r * sigma_r) << " J s\n"
		<< std::endl;

	const double vstd_r = std::sqrt(2 * kB_r * T_r / m_r); // m*s^-1

	// Problem constants
	// constexpr double m = 1.;
	// constexpr double eps = 1.;
	// constexpr double sigma = 1.;
	const double dt = dt_r / std::sqrt(m_r * sigma_r * sigma_r / eps_r);
	const double L = n * std::pow(1. / rho, 1. / 3);
	const double vstd =
		vstd_r / (sigma_r / std::sqrt(m_r * sigma_r * sigma_r / eps_r));

	// Info data
	std::cout
		<< "\nM:\t" << M
		<< "\nN:\t" << N
		<< "\nL:\t" << L
		<< "\nm:\t" << 1.
		<< "\neps:\t" << 1.
		<< "\nsigma:\t" << 1.
		<< "\nrho:\t" << rho
		<< "\ndt:\t" << dt
		<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
		<< "\nvstd:\t" << vstd
		<< std::endl;

	// Evolution variables
	std::vector<Vec3D> pos_a(N);
	std::vector<Vec3D> vel_a(N);
	std::vector<Vec3D> acc_a(N);
	std::vector<Vec3D> pos_b(N);
	std::vector<Vec3D> vel_b(N);
	std::vector<Vec3D> acc_b(N);

	// Proper pointers to handle them
	std::vector<Vec3D> *pos0 = &pos_a;
	std::vector<Vec3D> *pos1 = &pos_b;
	std::vector<Vec3D> *vel0 = &vel_a;
	std::vector<Vec3D> *vel1 = &vel_b;
	std::vector<Vec3D> *acc0 = &acc_a;
	std::vector<Vec3D> *acc1 = &acc_b;
	std::vector<Vec3D> *tmp;

	// Wanted quantities
	std::vector<double> Bs(M);

	// Boilerplate variables
	typedef struct
	{
		std::vector<Vec3D> **pos = nullptr;
		std::vector<Vec3D> **acc = nullptr;
		double (*Vr)(double) = nullptr;
		double (*dVr)(double) = nullptr;
		double (*ddVr)(double) = nullptr;
		double Vtot;
		double Q;
	} Args1;
	Args1 *args1;
	args1->pos = &pos1;
	args1->acc = &acc1;
	args1->Vr = V_LJ;
	args1->dVr = dV_LJ;
	args1->ddVr = ddV_LJ;
	int (*interaction)(
		const size_t i, const size_t j,
		const Vec3D &alias, const double d,
		Args1 *args1){
		[](const size_t i, const size_t j,
		   const Vec3D &alias, const double d,
		   Args1 *args1) -> int {
			// Obviously update accelerations
			(**args1->acc)[i] +=
				-args1->dVr(d) / d * ((**args1->pos)[i] - alias);
			// Addup the potential
			if (i < j)
				args1->Vtot += args1->Vr(d);
			// Compute the quantity
			args1->Q += args1->ddVr(d) + args1->dVr(d) * 2 / d;
			return 0;
		}};

	typedef struct
	{
		std::vector<Vec3D> **acc = nullptr;
	} Args2;
	Args2 *args2;
	args2->acc = &acc1;
	int (*start_interaction)(
		const size_t i,
		Args2 *args2){
		[](const size_t i,
		   Args2 *args2) -> int {
			(**args2->acc)[i].clear();
			return 0;
		}};

	// STARTING
	// Initialize positions with uniform lattice conditions
	std::cout
		<< "Init lattice" << std::endl;
	init_lattice(*pos0, L, n, 1);
	apply_periodic_bounds(*pos0, L);

	// Compute Potential and acceleration in one shot
	new_acceleration(
		*acc0, *pos0, V_LJ, 1, L, dV_LJ);

	// Initialize velocities as gaussian on the components
	std::cout << "Init velocities" << std::endl;
	for (uint32_t i = 0; i < N; i++)
	{
		init_distribute_maxwell_boltzmann((*vel0)[i], vstd);
	}

	// Main loop of the evolution
	std::cout << "Start velocity verlet" << std::endl;
	for (uint32_t i = 1; i < M; i++)
	{
		// Print progress
		if ((i % 1000) == 0)
		{
			std::cout
				<< "Step: " << i
				<< "\t Perc: " << (i * 100) / M
				<< "\tQ: " << 1.;
		}

		// Velocity - Verlet
		// Compute all new positions
		for (uint32_t j = 0; j < N; j++)
		{
			// Position step
			(*pos1)[j] = (*pos0)[j] +
						 dt * (*vel0)[j] +
						 0.5 * dt * dt * (*acc0)[j];

			// Periodic condition
			(*pos1)[j].x -= L * std::floor((*pos1)[j].x / L);
			(*pos1)[j].y -= L * std::floor((*pos1)[j].y / L);
			(*pos1)[j].z -= L * std::floor((*pos1)[j].z / L);
		}

		// Compute new accelerations and quantities
		args1->Q = 0;
		args1->Vtot = 0;
		do_interactions(*pos1,
						interaction, args1,
						start_interaction, args2,
						L);

		// Compute all new velocities
		for (uint32_t j = 0; j < N; j++)
		{
			// Velocity step
			(*vel1)[j] = (*vel0)[j] +
						 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);
		}

		// Quantities analysis
		args1->Q /= N;
		Bs[i] = args1->Q / 48;

		// Swap pointers for next step
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

	return 0;
}
