#include "simulation.h"
#include "differentiate.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

double Vf(const std::vector<Vec3D> &pos);
double V_LJ(double x);
void compute_acc_k(
	Vec3D &acc,
	const std::function<double(double)> Vr,
	const std::vector<Vec3D> &pos,
	double m,
	size_t k,
	double L);
void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f);

#undef SAVE_POS

int main(int argc, char const *argv[])
{
	uint32_t n = 5; // number of celles on side
	uint32_t N = n * n * n;
	uint32_t M = 1000000;

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
	dt_r /= 1000;
	// T_r = 10;
	M = 50000;

	// Problem quantities
	double dl{dl_r / sigma_r};									  // Cell length
	double L{n * dl};											  // Box length
	double dt{std::sqrt(eps_r / m_r / sigma_r / sigma_r) * dt_r}; // Time step
	double rho{std::pow(sigma_r, 3) * rho_r};					  // Density
	double m{rho * std::pow(dl, 3)};							  // Mass of the particles
	double vstd{vstd_r * std::sqrt(m_r / eps_r)};				  // std velocity

	// Print info data
	printf("dl:\t%f\nL:\t%f\ndt:\t%f\nrho:\t%f\nm:\t%f\nvstd:\t%f\n", dl, L, dt, rho, m, vstd);

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

	printf("Start\n");

	// Initialize with uniform lattice conditions
	printf("Init lattice\n");
	init_lattice(*pos0, dl, n, 1);
	apply_periodic_bounds(*pos0, L);

#ifdef SAVE_POS
	save_positions(*pos0, pos_file);
	pos_file << "\n\n";
#endif

	V = Vf(*pos0);

	// Initialize velocities as gaussian on the components
	printf("Init velocities\n");
	for (uint32_t i = 0; i < N; i++)
	{
		compute_acc_k((*acc0)[i], V_LJ, *pos0, m, i, L);
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

#ifdef SAVE_POS
		save_positions(*pos0, pos_file);
		pos_file << "\n\n";
#endif

		V = Vf(*pos1);
		// Compute all new velocities
		T = 0;
		for (uint32_t j = 0; j < N; j++)
		{
			compute_acc_k((*acc1)[j], V_LJ, *pos1, m, j, L);
			(*vel1)[j] = (*vel0)[j] + 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);
			T += 0.5 * m * (*vel1)[j].norm2();
		}
		E_tot = T + V;

		// TODO: write out pos1,vel1,energy
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
	if (x < 1e-10)
	{
		x = 1e-10;
	}
	return 4 * (std::pow(1.0 / x, 12) - std::pow(1.0 / x, 6));
}
double V_LJ_1(const Vec3D &a, const Vec3D &b)
{
	double d2 = (a - b).norm2();
	return 4 * (std::pow(1.0 / d2, 6) - std::pow(1.0 / d2, 3));
}

double Vf(const std::vector<Vec3D> &pos)
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
  Vij the potential between two particles
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

// Compute the acceleration on particle i
void compute_acc_k(
	Vec3D &acc,
	const std::function<double(double)> Vr,
	const std::vector<Vec3D> &pos,
	double m,
	size_t k,
	double L)
{
	// Reset previous acceleration
	acc = 0 * acc;
	size_t N{pos.size()};
	double d;
	// x=r_ij-r_ij0 ==> dx=dr
	// move the potential with the original position in the origin
	// this way the derivative can be computed in 0 on variable x
	// which directly varies the radial distance
	auto v = [&](double x) -> double {
		return Vr(d + x);
	};

	for (size_t i{0}; i < N; i++)
	{
		// Skip the self-case
		if (i != k)
		{
			// TODO: check this when I'm more active
			// Check  eventual alias
			// The partile _can_ interact only one time as real or aliased
			Vec3D alias = pos[i];
			d = (pos[k] - alias).norm();

			if (d > L / 2)
			{
				// calculate alias
				alias.x += (pos[k].x < L / 2) ? ((alias.x < L / 2) ? 0 : -L) : ((alias.x < L / 2) ? L : 0);
				alias.y += (pos[k].y < L / 2) ? ((alias.y < L / 2) ? 0 : -L) : ((alias.y < L / 2) ? L : 0);
				alias.z += (pos[k].z < L / 2) ? ((alias.z < L / 2) ? 0 : -L) : ((alias.z < L / 2) ? L : 0);
				d = (pos[k] - alias).norm();
			}

			// Check interaction limit
			if (d < L / 2)
			{
				// Compute the derivative on function v around 0
				double der = derive_5points(v, 0, 1e-10);
				// add the contribute to the force using the difference versor
				acc += (-der / d) * (pos[k] - alias);
			}
			// else no interaction:
			// the case of unperfect packing of spheres
		}
	}
}

void save_positions(const std::vector<Vec3D> &pos, std::ofstream &f)
{
	size_t N{pos.size()};
	for (size_t i{0}; i < N; i++)
	{
		f << i << " " << pos[i].x << " " << pos[i].y << " " << pos[i].z << "\n";
	}
}

double V_LJ_diff(double x)
{
	return -4 * (12 / x * std::pow(1.0 / x, 12) - 6 / x * std::pow(1.0 / x, 6));
}
