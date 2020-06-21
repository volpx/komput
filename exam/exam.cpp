#include "vec3d.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>

// Interaction potential
double V_LJ(double x)
{
	return 4 * (std::pow(1.0 / x, 12) - std::pow(1.0 / x, 6));
}
// Derivative of the potential
double dV_LJ(double x)
{
	return -4 * (12.0 / std::pow(x, 13) - 6.0 / std::pow(x, 7));
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
	constexpr double kB_r = 1;	 // J
	constexpr double hbar_r = 1; // J s

	// Argon
	constexpr double sigma_r = 3.398e-10; // m
	constexpr double T_r = 122;			  // K
	constexpr double eps_r = T_r * kB_r;  // J
	constexpr double m_r = 1;			  // Kg
	constexpr double lambda = 0.86e-3;	  //

	std::cout
		<< "Scale constants:\n"
		<< "\tKl: " << sigma_r << "m\n"
		<< "\tKe: " << eps_r << "J\n"
		<< "\tKrho: " << 1 / (sigma_r * sigma_r * sigma_r) << "m^-3\n"
		<< "\tKmrho: " << m_r / (sigma_r * sigma_r * sigma_r) << "Kg m^-3\n"
		<< "\tKmrho: " << m_r / (sigma_r * sigma_r * sigma_r) << "Kg m^-3\n"
		<< "\tKT: " << eps_r / kB_r << "K\n"
		<< "\tKP: " << eps_r / (sigma_r * sigma_r * sigma_r) << "Pa\n"
		<< "\tKt: " << std::sqrt(m_r * sigma_r * sigma_r / eps_r) << "s\n"
		<< "\tKh: " << std::sqrt(eps_r * m_r * sigma_r) << "J s\n"
		<< std::endl;

	const double vstd_r = std::sqrt(2 * kB_r * T_r / m_r); // m*s^-1

	// Problem constants
	constexpr double L = 1;
	const double vstd =
		vstd_r / (sigma_r / std::sqrt(m_r * sigma_r * sigma_r / eps_r));

	// Info data
	std::cout
		<< "\nM:\t" << M
		<< "\nN:\t" << N
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

	return 0;
}
