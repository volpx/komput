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
	// Number of temps
	constexpr int Ts_n = 6;
	// Number of rhos
	constexpr int rhos_n = 10;
	// Number of jobs
	constexpr int jobs_n = rhos_n * Ts_n;
	// Number of cells on side
	constexpr uint32_t n = 3;
	// Bravais lattice type
	constexpr uint32_t q = 4;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n * q;
	// Number of time steps
	constexpr uint32_t M = 10000;

	// Other lengths and simulation partitioning
	// Length of correlation
	constexpr uint32_t ac_length_n = 100;
	// Thermalization length
	constexpr uint32_t thermalization_n = 500;
	constexpr uint32_t Mt = M - thermalization_n;

	// Problem constants
	// constexpr double m = 1.;
	// constexpr double eps = 1.;
	// constexpr double sigma = 1.;
	// Evolution timestep,
	// kept it fixed but could be a function of L
	constexpr double dt = 0.0045;
	constexpr double Ms = 10;
	constexpr double gs = 3 * N + 1;

	// Job subdivision
	// vary T from 0.8 to 1.3
	std::vector<double> Ts(Ts_n);
	linspace(Ts, 0.8, 1.3);
	// Ts[0] = 1.;
	// vary rho from 0.65 to 0.9
	std::vector<double> rhos(rhos_n);
	linspace(rhos, 0.55, 1.);
	// fill(rhos, 0.65);

	// Create jobs
	struct Job
	{
		const int id;
		const double Temp;
		const double rho;
	};
	std::vector<Job> jobs;
	jobs.reserve(jobs_n);
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
		double (*const Vr)(double) = nullptr;
		double (*const dVr)(double) = nullptr;
		double (*const ddVr)(double) = nullptr;
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
	// std::ofstream B_file("data/B.dat",
	// 					 std::ofstream::out | std::ofstream::app);
	std::ofstream B_file{"data/B.dat"};
	B_file << "#T rho B stdB\n";

#ifdef DEBUG_SAVE
	std::ofstream data_file("data/data_evolution.dat");
	data_file << "#i,t,T,V,E,Q\n";
#endif

	std::vector<double> tac(ac_length_n);
	map(tac, [](uint32_t i) -> double {
		return i;
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
	std::vector<Vec3D> *tmp = nullptr;

	// Wanted quantities accumulator
	std::vector<double> Qs(Mt);
	std::vector<double> Qac(ac_length_n);
	std::vector<double> Ks(M);

	// job loops
	for (const auto &job : jobs)
	{
		// Seed for the job
		srand(job.id + 420);

		// Extract job data
		std::cout << "\n\nJob id: " << job.id
				  << " Progress jobs: " << job.id * 100 / jobs.size()
				  << std::endl;
		const double Temp = job.Temp;
		const double rho = job.rho;

		// Derivatives of job data
		const double vstd = std::sqrt(2. * Temp);
		const double L = n * std::pow(1. / (rho / q), 1. / 3);

#ifdef LJ_CORRECTION_OFFSET
		const double V_offset = -V_LJ(L / 2);
#else
		constexpr double V_offset = 0;
#endif

		// Info data
		std::cout
			<< "\nM:\t" << M
			<< "\nMt:\t" << Mt
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

		IntCall interaction{
			&pos1, &acc1, V_LJ, dV_LJ, ddV_LJ, V_offset};
		IntSetup start_interaction{&acc1};

		double E, T, Temp_ist;
		double s0, s1, vs0, vs1, as0, as1;
		double &V = interaction.V;
		double &Q = interaction.Q;

		// STARTING

		// Initialize positions with uniform lattice conditions
		std::cout << "Init lattice" << std::endl;
		init_lattice(*pos0, L, n, q);
		apply_periodic_bounds(*pos0, L);
		s0 = 1;

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
		as0 = 0;

		// Initialize velocities as gaussian on the components
		T = 0;
		std::cout << "Init velocities" << std::endl;
		for (uint32_t i = 0; i < N; i++)
		{
			init_distribute_maxwell_boltzmann((*vel0)[i], vstd);
			T += 0.5 * (*vel0)[i].norm2();
		}
		vs0 = 0;

		// Quantities analysis
		E = T + V;
		Q /= N;
		Ks[0] = T;

		// Thermalizaton quantities
		if (0 >= thermalization_n)
		{
			Qs[0] = Q;
		}

		// Output
		std::cout
			<< "First data:"
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
			std::cout << "aa" << vs0 << " " << vs1 << std::endl;
			// Velocity - Verlet
			// Compute all new positions
			// for (uint32_t j = 0; j < N; j++)
			// {
			// 	// Position step
			// 	(*pos1)[j] = (*pos0)[j] +
			// 				 dt * (*vel0)[j] +
			// 				 0.5 * dt * dt * (*acc0)[j];

			// 	// Periodic condition
			// 	(*pos1)[j].x -= L * std::floor((*pos1)[j].x / L);
			// 	(*pos1)[j].y -= L * std::floor((*pos1)[j].y / L);
			// 	(*pos1)[j].z -= L * std::floor((*pos1)[j].z / L);
			// }
			s1 = s0 + dt * vs0 + 0.5 * dt * dt * as0;
			std::cout << "a" << vs0 << " " << vs1 << std::endl;

			// Compute new accelerations and quantities
			// Q = 0;
			// V = 0;
			// do_interactions(*pos1,
			// 				&interaction,
			// 				&start_interaction,
			// 				L);
			std::cout << "b" << vs0 << " " << vs1 << std::endl;

			// Compute all new velocities
			T = 0;
			// for (uint32_t j = 0; j < N; j++)
			// {
			// 	// Velocity step
			// 	(*vel1)[j] = (*vel0)[j] +
			// 				 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);
			// 	(*vel1)[j] /= (s1 * s1);
			// 	T += 0.5 * (*vel1)[j].norm2();
			// }
			std::cout << "c" << vs0 << " " << vs1 << std::endl;
			as1 = -(T * (-2 / s1) + gs * Temp / s1) / Ms;
			std::cout << "d" << vs0 << " " << vs1 << std::endl;
			vs1 += +0.5 * dt * (as0 + as1);
			std::cout << "e" << vs0 << " " << vs1 << std::endl;

			// Quantities analysis
			// 			E = T + V;
			// 			Q /= N;
			// 			Ks[i] = T;
			// 			Temp_ist = T * 2 / (3 * N - 3);

			// 			// Thermalization quantities
			// 			if (i >= thermalization_n)
			// 			{
			// 				Qs[i - thermalization_n] = Q;
			// 			}

			// 			// Print progress
			// 			if ((i % 1000) == 0)
			// 			{
			// 				std::cout
			// 					<< "Step: " << i
			// 					<< "\tP: " << (i * 100) / M
			// 					<< "\tt: " << i * dt
			// 					<< "\tE: " << E
			// 					<< "\tT: " << T
			// 					<< "\tV: " << V
			// 					<< "\tQ: " << Q
			// 					<< "\tTemp_ist: " << Temp_ist
			// 					<< "\t" << vs0
			// 					<< "\t" << vs1
			// 					<< std::endl;
			// 			}

			// 			// Output
			// #ifdef DEBUG_SAVE
			// 			data_file
			// 				<< i << ' '
			// 				<< i * dt << ' '
			// 				<< T << ' '
			// 				<< V << ' '
			// 				<< E << ' '
			// 				<< Q << '\n';
			// #endif

			// 			// Swap pointers for next step
			// 			tmp = pos1;
			// 			pos1 = pos0;
			// 			pos0 = tmp;

			// 			tmp = vel1;
			// 			vel1 = vel0;
			// 			vel0 = tmp;

			// 			tmp = acc1;
			// 			acc1 = acc0;
			// 			acc0 = tmp;

			s0 = s1;
			std::cout << "f" << vs0 << " " << vs1 << std::endl;
			vs0 = vs1;
			std::cout << "g" << vs0 << " " << vs1 << std::endl;
			as0 = as1;
			getchar();
		} // End main VV loop

		// Post analysis
		std::cout << "Post analysis" << std::endl;

		// B autocorrelation
		autocorrelation(Qac, Qs);
		// Correlation factor calculation with fit
		double par[3] = {1., 0.1, 0.};
		fit_to_exp(par, tac, Qac);
		// Number of dependent points
		double tau = 1. / par[1];

#ifdef DEBUG_SAVE
		std::ofstream Qac_file("data/Qautocorr.dat");
		Qac_file << "#"
				 << " par0:" << par[0]
				 << " par1:" << par[1]
				 << " par2:" << par[2]
				 << " tau:" << tau
				 << "\n";
		Qac_file << "#i,t,Qac\n";
		for (uint32_t i = 0; i < Qac.size(); ++i)
		{
			Qac_file << i << ' ' << tac[i] << ' ' << Qac[i] << '\n';
		}
		Qac_file.close();
#endif

		double meanQ = mean(Qs);
		double DmeanQ = tau / Qs.size() * variance(Qs, 1);
		// DmeanQ *= 25; // Magic factor

		std::cout
			<< "T0: " << Temp
			<< " rho: " << rho
			<< " B: " << meanQ / 48
			<< " +/- " << std::sqrt(DmeanQ) / 48
			<< " tau: " << tau
			<< " Temp_mean: " << mean(Ks) * 2 / (3 * N - 3)
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
