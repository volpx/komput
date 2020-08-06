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
	constexpr uint32_t Ts_n = 10;
	// Number of rhos
	constexpr uint32_t rhos_n = 50;
	// constexpr uint32_t Ts_n = 1;
	// constexpr uint32_t rhos_n = 1;
	// Number of jobs
	constexpr uint32_t jobs_n = rhos_n * Ts_n;
	// Number of cells on side
	constexpr uint32_t n = 3;
	// Bravais lattice type
	constexpr uint32_t q = 4;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n * q;
	// Degrees of freedom
	constexpr double gs = 3 * N + 1;
	// Number of time steps
	constexpr uint32_t M = 3000;

	// Other lengths and simulation partitioning
	// Length of correlation
	constexpr uint32_t ac_length_n = 200;
	// Thermalization length
	constexpr uint32_t thermalization_n = 500;
	constexpr uint32_t Mt = M - thermalization_n;

	// Problem constants
	// constexpr double m = 1.;
	// constexpr double eps = 1.;
	// constexpr double sigma = 1.;

	// Evolution timestep,
	// kept it fixed but could be a function of L
	constexpr double dt = 0.004;

	// Job subdivision
	// vary T from 0.8 to 1.3
	std::vector<double> Ts(Ts_n);
	linspace(Ts, 0.8, 1.5);
	// Ts[0] = 1.4;
	// vary rho from 0.65 to 0.9
	std::vector<double> rhos(rhos_n);
	linspace(rhos, 0.6, 0.9);
	// linspace(rhos, 0.87, 0.895);
	// fill(rhos, 0.65);
	// rhos[0] = 1;

	// Create jobs
	struct Job
	{
		const uint32_t id;
		const double T;
		const double rho;
	};
	std::vector<Job> jobs;
	jobs.reserve(jobs_n);
	for (uint32_t i = 0; i < Ts_n; i++)
	{
		for (uint32_t j = 0; j < rhos_n; j++)
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
		const double V_offset;
		double V;
		double Q;

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
	std::ofstream B_file{"data/B.dat"};
	B_file << "#T rho B stdB\n";

#ifdef DEBUG_SAVE
	std::ofstream data_file("data/data_evolution.dat");
	data_file << "#i,t,T,V,E,Q,Temp_ist,s,vs,as\n";
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

	// special support
	double E, K, K_half, T_ist;
	double s0 = 1, s1 = 1, vs0 = 0, vs1 = 0, as0 = 0, as1 = 0;

	// job loops
	for (const auto &job : jobs)
	{
		// Seed for the job
		srand(job.id + 420);

		// Extract job data
		std::cout << "\n\nJob id: " << job.id
				  << " Progress jobs: " << job.id * 100 / jobs.size()
				  << std::endl;
		const double T = job.T;
		const double rho = job.rho;

		// Derivatives of job data
		const double vstd = std::sqrt(T);
		const double L = n * std::pow(1. / (rho / q), 1. / 3);
		// const double Ms = 0.01 * gs * T;
		// const double Ms = 0.005 * gs * T;
		const double Ms = 1.5;

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
			<< "\nT:\t" << T
			<< "\nMs:\t" << Ms
			<< "\ndt:\t" << dt
			<< "\nvstd:\t" << vstd
			<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
#ifdef LJ_CORRECTION_OFFSET
			<< "\nLJ offest enabled"
#endif
			<< '\n'
			<< std::endl;

		IntCall interaction{
			&pos1, &acc1, V_LJ, dV_LJ, ddV_LJ, V_offset, 0, 0};
		IntSetup start_interaction{&acc1};
		double &V = interaction.V;
		double &Q = interaction.Q;

		// STARTING
		s1 = s0 = 1;
		vs0 = vs1 = 0;
		as1 = as0 = 0;

		// Initialize positions with uniform lattice conditions
		std::cout << "Init lattice" << std::endl;
		init_lattice(*pos0, L, n, q);
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
		K = 0;
		K_half = 0;
		std::cout << "Init velocities" << std::endl;
		init_distribute_maxwell_boltzmann((*vel0), vstd);
		for (uint32_t j = 0; j < N; j++)
		{
			K += (*vel0)[j].norm2();
			K_half += ((*vel0)[j] + 0.5 * dt * (*acc0)[j]).norm2();
		}
		K *= 0.5 / (s1 * s1);
		K_half *= 0.5 / (s1 * s1);
		as0 = 1.0 / (Ms * s0) * (2 * K - gs * T);
		as1 = 1.0 / (Ms * s1) * (2 * K_half - gs * T);
		vs1 = 0;
		vs0 = 0.5 * dt * (as0 + as1);

		// Quantities analysis
		E = K + V;
		E += 0.5 * Ms * vs1 * vs1 + gs * T * std::log(s1);
		Q /= N;
		Ks[0] = K;
		T_ist = K * 2 / gs;

		// Thermalizaton quantities
		if (0 >= thermalization_n)
		{
			Qs[0] = Q;
		}

		// Output
		std::cout
			<< "First data:"
			<< "\tE: " << E
			<< "\tK: " << K
			<< "\tV: " << V
			<< "\tQ: " << Q
			<< "\tT_ist: " << T_ist
			<< std::endl;

#ifdef DEBUG_SAVE
		data_file
			<< 0 << ' '
			<< 0 << ' '
			<< K << ' '
			<< V << ' '
			<< E << ' '
			<< Q << ' '
			<< T_ist << ' '
			<< s1 << ' '
			<< vs1 << ' '
			<< as1 << '\n';
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
							 dt / (s0 * s0) * (*vel0)[j] +
							 0.5 * dt * dt / (s0 * s0) * (*acc0)[j];

				// Periodic condition
				(*pos1)[j].x -= L * std::floor((*pos1)[j].x / L);
				(*pos1)[j].y -= L * std::floor((*pos1)[j].y / L);
				(*pos1)[j].z -= L * std::floor((*pos1)[j].z / L);
			}
			as0 = 1.0 / (Ms * s0) * (2 * K - gs * T);
			s1 = s0 + dt * vs0 + 0.5 * dt * dt * as0;

			// Compute new accelerations and quantities
			Q = 0;
			V = 0;
			do_interactions(*pos1,
							&interaction,
							&start_interaction,
							L);

			// Compute all new velocities
			K = 0;
			K_half = 0;
			for (uint32_t j = 0; j < N; j++)
			{
				// Velocity step
				(*vel1)[j] = (*vel0)[j] + 0.5 * dt * ((*acc0)[j] + (*acc1)[j]);

				K += (*vel1)[j].norm2();
				K_half += ((*vel0)[j] + 0.5 * dt * (*acc0)[j]).norm2();
			}
			K *= 0.5 / (s1 * s1);
			K_half *= 0.5 / (s1 * s1);
			as1 = 1.0 / (Ms * s1) * (2 * K_half - gs * T);
			vs1 = vs0 + 0.5 * dt * (as0 + as1);

			// Quantities analysis
			E = K + V;
			E += 0.5 * Ms * vs1 * vs1 + gs * T * std::log(s1);
			Q /= N;
			Ks[i] = K;
			T_ist = K * 2 / gs;

			// Thermalization quantities
			if (i >= thermalization_n)
			{
				Qs[i - thermalization_n] = Q;
			}

			// Print progress
			if ((i % 1000) == 0)
			{
				std::cout
					<< "Step: " << i
					<< "\tP: " << (i * 100) / M
					<< "\tt: " << i * dt
					<< "\tE: " << E
					<< "\tK: " << K
					<< "\tV: " << V
					<< "\tQ: " << Q
					<< "\tT_ist: " << T_ist
					<< std::endl;
			}

			// Output
#ifdef DEBUG_SAVE
			data_file
				<< i << ' '
				<< i * dt << ' '
				<< K << ' '
				<< V << ' '
				<< E << ' '
				<< Q << ' '
				<< T_ist << ' '
				<< s1 << ' '
				<< vs1 << ' '
				<< as1 << '\n';
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

			s0 = s1;
			vs0 = vs1;
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

		std::cout
			<< "T0: " << T
			<< " T_mean: " << mean(Ks) * 2 / gs
			<< " rho: " << rho
			<< " B: " << meanQ / 48
			<< " +/- " << std::sqrt(DmeanQ) / 48
			<< " tau: " << tau
			<< std::endl;

		B_file
			<< T << ' '
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
