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

#undef SAVE_POS

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

// Save positions to file
void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f);

int main(int argc, char const *argv[])
{
	// Number of cells on side
	uint32_t n = 5;
	// Number of particles in the box
	uint32_t N = n * n * n;
	// Number of time steps
	uint32_t M = 100000;

	// Physics dimensional quantities
	double e_r = 1.6e-19;	   // C
	double kB_r = 1.38e-23;	   // J*K^-1
	double T_r = 300;		   // K
	double m_r = 2.2e-25;	   // Kg
	double rho_r = 1e27;	   // #*m^-3
	double sigma_r = 3.94e-10; // m
	double dt_r = 1e-12;	   // s
	double eps_r = 0.02 * e_r; // J
	// rho_r *= 1 / 0.061163;
	double dl_r = std::pow(rho_r, -1.0 / 3);	 // m
	double vstd_r = std::sqrt(kB_r * T_r / m_r); // m*s^-1

	// Modifications
	dt_r /= 100;
	// T_r = 10;
	// m_r *= 1 / 8.17488;
	// dt_r *= 0.01 / 0.875202;
	// dl_r *= 6.29960525 / 12.6904;
	// vstd_r *= 0;

	// Problem adimensional quantities
	double dl{dl_r / sigma_r};									  // Cell length
	double L{n * dl};											  // Box length
	double dt{std::sqrt(eps_r / m_r / sigma_r / sigma_r) * dt_r}; // Time step
	double rho{std::pow(sigma_r, 3) * rho_r};					  // Density
	double m{rho * std::pow(dl, 3)};							  // Mass of the particles
	double vstd{vstd_r * std::sqrt(m_r / eps_r)};				  // std velocity
	double Temp{kB_r / eps_r * T_r};							  // Temperature

	// Autocorrelation velocity variables
	// Starting index
	uint32_t cvv_istart{1000};
	// Final correlation length
	uint32_t cvv_length{500};
	// Instantiate a AutoCorr class
	AutoCorr cvv(cvv_istart, M, cvv_length, N);

	// Print info data
	std::cout
		<< "\nM:\t" << M
		<< "\nN:\t" << N
		<< "\nL:\t" << L
		<< "\ndl:\t" << dl
		<< "\nrho:\t" << rho
		<< "\nm:\t" << m
		<< "\nT:\t" << Temp
		<< "\ndt:\t" << dt
		<< "\nM:\t" << M
		<< "\nTend:\t" << dt * M
		<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
		<< "\nvstd:\t" << vstd
		<< "\ncvv_istart\t" << cvv_istart
		<< "\ncvv_means_number:\t" << cvv.get_means_number()
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

	// Energy variables
	double E_tot{0};
	double T{0};
	double V{0};

#ifdef SAVE_POS
	std::ofstream pos_file("output_data/positions.dat");
	pos_file << "#n x y z\n";
#endif

	// Prepare energy output
	std::ofstream energy_file("output_data/energy.dat");
	energy_file << "#m t T V E\n";

	// Initialize positions with uniform lattice conditions
	std::cout << "Init lattice" << std::endl;
	init_lattice(*pos0, L, n, 1);
	apply_periodic_bounds(*pos0, L);

#ifdef SAVE_POS
	save_positions(*pos0, pos_file);
	pos_file << "\n\n";
#endif

	// Compute Potential and acceleration in one shot
	V = new_acceleration(
		*acc0, *pos0, V_LJ, m, L, dV_LJ);
	std::cout << "V0: " << V << std::endl;

	// Initialize velocities as gaussian on the components
	std::cout << "Init velocities" << std::endl;
	T = 0;
	for (uint32_t i = 0; i < N; i++)
	{
		init_distribute_maxwell_boltzmann((*vel0)[i], vstd);
		// Add contribute to kinetic energy
		T += 0.5 * m * (*vel0)[i].norm2();
	}
	E_tot = T + V;
	std::cout << "Initial energy: " << E_tot << std::endl;

	// Save energy
	energy_file << 0 << " " << 0 << " "
				<< T << " " << V << " " << E_tot << "\n";

	// Main loop of the evolution
	std::cout << "Start velocity verlet" << std::endl;
	for (uint32_t i = 1; i < M; i++)
	{
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

		// Compute potential and new accelerations
		V = new_acceleration(
			*acc1, *pos1, V_LJ, m, L, dV_LJ);

		// Compute all new velocities and add-up the kinetic energy
		T = 0;
		for (uint32_t j = 0; j < N; j++)
		{
			// Velocity step
			(*vel1)[j] = (*vel0)[j] +
						 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);

			T += 0.5 * m * (*vel1)[j].norm2();
		}
		E_tot = T + V;

		// CVV add the data for autocorrelation
		cvv.add_data(*vel1, i);

		// Output
#ifdef SAVE_POS
		save_positions(*pos0, pos_file);
		pos_file << "\n\n";
#endif
		// Save energy
		energy_file << i << " " << i * dt << " "
					<< T << " " << V << " " << E_tot << "\n";

		// Print progress
		if ((i % 1000) == 0)
		{
			printf("Step: %7d \t Percentage: %3d \tEnergy: %f\n",
				   i, (i * 100) / M, E_tot);
		}

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

	// Close the energy file
	energy_file.close();

	// Write out the cvv
	std::ofstream cvv_file("output_data/cvv.dat");
	cvv_file << "#t cvv\n";
	for (uint32_t i{0}; i < cvv_length; i++)
	{
		cvv_file << dt * i << " " << cvv[i] << "\n";
	}
	cvv_file.close();

	std::cout << "Final energy: " << E_tot << std::endl;
	std::cout << "End" << std::endl;

	return 0;
}

// Save positions in file for eventual plot
void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f)
{
	size_t N{pos.size()};
	for (size_t i{0}; i < N; i++)
	{
		f << i
		  << " " << pos[i].x
		  << " " << pos[i].y
		  << " " << pos[i].z << "\n";
	}
}
