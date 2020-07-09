#include "vec3d.h"
#include "simulation.h"
#include "optimize.h"
#include "uniconst.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>

#define LJ_CORRECTION_OFFSET
// #define DEBUG_SAVE

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
	// Number of threads
	constexpr int threads_n = 2;
	// Number of temps
	constexpr int Ts_n = 10;
	// Number of rhos
	constexpr int rhos_n = 50;
	// Number of cells on side
	constexpr uint32_t n = 5;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n;
	// Number of time steps
	constexpr uint32_t M = 10'000;
	// Length of correlation
	constexpr uint32_t i_ac_length = 500;

	// Problem constants
	// constexpr double m = 1.;
	// constexpr double eps = 1.;
	// constexpr double sigma = 1.;
	// Evolution timestep,
	// kept it fixed but could be a function of L
	constexpr double dt = 0.0045;

	// vary T from 0.8 to 1.3
	std::vector<double> Ts(Ts_n);
	linspace(Ts, 0.8, 1.3);
	// vary rho from 0.65 to 0.9
	std::vector<double> rhos(rhos_n);
	linspace(rhos, 0.65, 0.9);

	// Create jobs
	struct Job
	{
		const int id;
		const double Temp;
		const double rho;
	};
	std::vector<Job> jobs;
	jobs.reserve(Ts_n * rhos_n);
	for (int i = 0; i < Ts_n; i++)
	{
		for (int j = 0; j < rhos_n; j++)
		{
			jobs.emplace_back(Job{i * rhos_n + j, Ts[i], rhos[j]});
		}
	}

	// Interaction description
	struct IntCall
	{
		std::vector<Vec3D> **pos = nullptr;
		std::vector<Vec3D> **acc = nullptr;
		double (*Vr)(double) = nullptr;
		double (*dVr)(double) = nullptr;
		double (*ddVr)(double) = nullptr;
		double V;
		double Q;
		const double V_offset;

		void operator()(const size_t i, const size_t j,
						const Vec3D &alias, const double d)
		{
			// Obviously update accelerations
			(**acc)[i] +=
				-dVr(d) / d * ((**pos)[i] - alias);
			// Addup the potential
			if (i < j)
				V += (Vr(d) + V_offset);
			// Compute the quantity
			Q += ddVr(d) + dVr(d) * 2 / d;
		}
	};
	struct IntSetup
	{
		std::vector<Vec3D> **acc = nullptr;
		void operator()(const size_t i)
		{
			// Clear the acceleration before start adding contributions
			(**acc)[i].clear();
		}
	};

	// Open main file of output
	std::ofstream B_file("data/B.dat",
						 std::ofstream::out | std::ofstream::app);
	B_file << "#T rho B stdB\n";

#ifdef DEBUG_SAVE
	std::ofstream data_file("data/data_evolution.dat");
	data_file << "#i,t,T,V,E,Q\n";
#endif

	std::vector<double> tac(i_ac_length);
	map(tac, [dt](uint32_t i) -> double {
		return i * dt;
	});

	// In case of multi threading
	// from here on the variables are thread specific

	// Evolution status variables
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

	// Wanted quantities accumulator
	std::vector<double> Qs(M);
	std::vector<double> Qac(i_ac_length);

	// job loops
	for (auto job : jobs)
	{
		// Extract job data
		std::cout << "Job id: " << job.id
				  << " Progress jobs: " << job.id * 100 / jobs.size()
				  << std::endl;
		const double Temp = job.Temp;
		const double rho = job.rho;

		// Derivatives of job data
		const double vstd = std::sqrt(2 * Temp);
		const double L = n * std::pow(1. / rho, 1. / 3);

#ifdef LJ_CORRECTION_OFFSET
		const double V_offset = -V_LJ(L / 2);
#endif

		// Info data
		std::cout
			<< "\nM:\t" << M
			<< "\nN:\t" << N
			<< "\nL:\t" << L
			<< "\nm:\t" << 1.
			<< "\neps:\t" << 1.
			<< "\nsigma:\t" << 1.
			<< "\nrho:\t" << rho
			<< "\nT:\t" << Temp
			<< "\ndt:\t" << dt
			<< "\nvstd:\t" << vstd
			<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
#ifdef LJ_CORRECTION_OFFSET
			<< "\nLJ offest enabled"
#endif
			<< '\n'
			<< std::endl;
#ifdef DEBUG_SAVE
		data_file << "# rho: " << rho << "\n";
#endif

		IntCall interaction{
			&pos1, &acc1, V_LJ, dV_LJ, ddV_LJ, V_offset};
		IntSetup start_interaction{&acc1};

		double E, T;
		double Qsum = 0, Q2sum = 0;
		double &V = interaction.V;
		double &Q = interaction.Q;

		// STARTING

		// Initialize positions with uniform lattice conditions
		std::cout << "Init lattice" << std::endl;
		init_lattice(*pos0, L, n, 1);
		apply_periodic_bounds(*pos0, L);

		// Compute new accelerations and quantities
		interaction.pos = &pos0;
		interaction.acc = &acc0;
		start_interaction.acc = &acc0;
		V = 0;
		Q = 0;
		do_interactions(*pos0,
						&interaction,
						&start_interaction,
						L);
		interaction.pos = &pos1;
		interaction.acc = &acc1;
		start_interaction.acc = &acc1;

		// Initialize velocities as gaussian on the components
		T = 0;
		std::cout << "Init velocities" << std::endl;
		for (uint32_t i = 0; i < N; i++)
		{
			init_distribute_maxwell_boltzmann((*vel0)[i], vstd);
			T += 0.5 * (*vel0)[i].norm2();
		}

		// Quantities analysis
		E = T + V;
		Q /= N;
		Qs[0] = Q;

		// Output
		std::cout
			<< "First data:\n"
			<< "Step: " << 0
			<< "\tP: " << (0 * 100) / M
			<< "\tt: " << 0
			<< "\tE: " << E
			<< "\tT: " << T
			<< "\tV: " << V
			<< "\tQ: " << Q
			<< std::endl;

#ifdef DEBUG_SAVE
		data_file
			<< 0 << ' '
			<< 0 << ' '
			<< T << ' '
			<< V << ' '
			<< E << ' '
			<< Q << '\n';
#endif

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

			// Compute new accelerations and quantities
			Q = 0;
			V = 0;
			do_interactions(*pos1,
							&interaction,
							&start_interaction,
							L);

			// Compute all new velocities
			T = 0;
			for (uint32_t j = 0; j < N; j++)
			{
				// Velocity step
				(*vel1)[j] = (*vel0)[j] +
							 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);
				T += 0.5 * (*vel1)[j].norm2();
			}

			// Quantities analysis
			E = T + V;
			Q /= N;
			Qs[i] = Q;

			// Print progress
			if ((i % 1000) == 0)
			{
				std::cout
					<< "Step: " << i
					<< "\tP: " << (i * 100) / M
					<< "\tt: " << i * dt
					<< "\tE: " << E
					<< "\tT: " << T
					<< "\tV: " << V
					<< "\tQ: " << Q
					<< std::endl;
			}

			// Output
#ifdef DEBUG_SAVE
			data_file
				<< i << ' '
				<< i * dt << ' '
				<< T << ' '
				<< V << ' '
				<< E << ' '
				<< Q << '\n';
#endif

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

		// Post analysis
		std::cout << "Post analysis" << std::endl;

		// B autocorrelation
		autocorrelation(Qac, Qs);
#ifdef DEBUG_SAVE
		std::ofstream Qac_file("data/Qautocorr.dat");
		Qac_file << "#i,t,Qac\n";
		for (uint32_t i = 0; i < Qac.size(); ++i)
		{
			Qac_file << i << ' ' << tac[i] << ' ' << Qac[i] << '\n';
		}
		Qac_file.close();
#endif

		// Correlation factor calculation with fit
		double par[3];
		fit_to_exp(par, tac, Qac);
		// Number of dependent points
		double tau = (1. / par[1]) / dt;

		double meanQ = mean(Qs);
		double DmeanQ = tau / (M - 1) * variance(Qs);

		std::cout
			<< "T: " << Temp
			<< "rho: " << rho
			<< " B: " << meanQ / 48
			<< " +/- " << std::sqrt(DmeanQ) / 48
			<< " tau: " << tau
			<< std::endl;

		B_file
			<< Temp << ' '
			<< rho << ' '
			<< meanQ / 48 << ' '
			<< std::sqrt(DmeanQ) / 48
			<< '\n';

#ifdef DEBUG_SAVE
		data_file << "\n\n\n";
#endif
	}
	B_file.close();

#ifdef DEBUG_SAVE
	data_file.close();
#endif

	return 0;
}
