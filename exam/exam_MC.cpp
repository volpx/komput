#include "vec3d.h"
#include "simulation.h"
#include "optimize.h"
#include "uniconst.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>

// #define LJ_CORRECTION_OFFSET
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
	constexpr uint32_t ac_length_n = 500;
	constexpr uint32_t Delta_correction_n = 100;

	// Problem constants
	// constexpr double m = 1.;
	// constexpr double eps = 1.;
	// constexpr double sigma = 1.;

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
	// jobs.reserve(1);
	// jobs.emplace_back(Job{0, 1.2, 0.7});
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
		double (*const Vr)(double) = nullptr;
		double (*const dVr)(double) = nullptr;
		double (*const ddVr)(double) = nullptr;
		const double V_offset;
		double V;
		double Q;

		void operator()(const size_t i, const size_t j,
						const Vec3D &alias, const double d)
		{
			// Addup the potential
			if (i < j)
				V += (Vr(d) + V_offset);
			// Compute the quantity
			Q += ddVr(d) + dVr(d) * 2 / d;
		}
	};
	struct IntSetup
	{
		void operator()(const size_t i)
		{
		}
	};

	// Open main file of output
	// std::ofstream B_file("data/B.dat",
	// 					 std::ofstream::out | std::ofstream::app);
	std::ofstream B_file{"data/B_MC.dat"};
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
	std::vector<Vec3D> pos_b(N);

	// Proper pointers to handle them
	std::vector<Vec3D> *pos0 = &pos_a;
	std::vector<Vec3D> *pos1 = &pos_b;
	std::vector<Vec3D> *tmp = nullptr;

	// Wanted quantities accumulator
	std::vector<double> Qs_all(M);
	std::vector<double> Qs(M);
	std::vector<double> Qac(ac_length_n);
	std::vector<double> As(Delta_correction_n);

	// job loop
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
		const double beta = 1. / Temp;
		const double L = n * std::pow(1. / (rho / q), 1. / 3);
		double Delta = L / n / 100;
		double lambda = 0.05;

#ifdef LJ_CORRECTION_OFFSET
		const double V_offset = -V_LJ(L / 2);
#else
		constexpr double V_offset = 0;
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
			<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
#ifdef LJ_CORRECTION_OFFSET
			<< "\nLJ offest enabled"
#endif
			<< '\n'
			<< std::endl;

		IntCall interaction{
			&pos1, V_LJ, dV_LJ, ddV_LJ, V_offset, 0, 0};
		IntSetup start_interaction{};

		double V1, V0;
		double &V = interaction.V;
		double &Q = interaction.Q;
		double Q0, A, A_mean;
		uint32_t thermalization_n{0};
		bool thermalized{false};

		// STARTING

		// Initialize positions with uniform lattice conditions
		std::cout << "Init lattice" << std::endl;
		init_lattice(*pos0, L, n, q);
		apply_periodic_bounds(*pos0, L);

		// Compute new accelerations and quantities
		interaction.pos = &pos0;
		V = 0;
		Q = 0;
		do_interactions(*pos0,
						&interaction,
						&start_interaction,
						L);
		interaction.pos = &pos1;

		// Quantities analysis
		V0 = V;
		Q0 = Q / N;
		As[0] = 0.5; // neutral value, i don't have any other data

		// Output
		std::cout
			<< "First data:"
			<< "\tV: " << V0
			<< "\tQ: " << Q0
			<< "\tDelta: " << Delta
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
		std::cout << "Start monte carlo" << std::endl;
		for (uint32_t i = 1; i < M; i++)
		{
			// Compute all new positions
			for (uint32_t j = 0; j < N; j++)
			{
				(*pos1)[j].x = (*pos0)[j].x + Delta * (randu() - 0.5);
				(*pos1)[j].y = (*pos0)[j].y + Delta * (randu() - 0.5);
				(*pos1)[j].z = (*pos0)[j].z + Delta * (randu() - 0.5);
				// Periodic condition
				(*pos1)[j].x -= L * std::floor((*pos1)[j].x / L);
				(*pos1)[j].y -= L * std::floor((*pos1)[j].y / L);
				(*pos1)[j].z -= L * std::floor((*pos1)[j].z / L);
			}

			// Compute quantities
			Q = 0;
			V = 0;
			do_interactions(*pos1,
							&interaction,
							&start_interaction,
							L);

			// New potential
			V1 = V;
			// Transition probability (save it)
			A = std::exp(-beta * (V1 - V0));
			As[i % Delta_correction_n] = A;

			// corrections to Delta
			if (!thermalized)
			{
				if (i % Delta_correction_n == 0)
				{
					A_mean = mean(As);
					if (A_mean > 0.4 && A_mean < 0.6)
					{
						thermalized = true;
						thermalization_n = i;
					}
					else
					{
						// Correction
						Delta *= 1 + lambda * (A_mean - 0.5);
					}
				}
			}

			// Decide if accept the new status
			// and put it in _0
			if (randu() < A)
			{
				// accept the new configuration
				tmp = pos0;
				pos0 = pos1;
				pos1 = tmp;
				V0 = V1;
				Q0 = Q / N;
			}
			else
			{
				// else the new sample is the same as before
				// no need of doing anything
			}

			// Quantities analysis
			// Thermalization quantities
			Qs_all[i] = Q0;

			// Print progress
			if ((i % 1000) == 0)
			{
				A_mean = mean(As);
				std::cout
					<< "Step: " << i
					<< "\tP: " << (i * 100) / M
					<< "\tV: " << V0
					<< "\tQ: " << Q0
					<< "\tDelta: " << Delta
					<< "\tA_mean: " << A_mean
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
		}

		// Post analysis
		std::cout << "Post analysis" << std::endl;

		// Throw away not termalized Q
		Qs.resize(M - thermalization_n);
		for (uint32_t i{0}; i < Qs.size(); ++i)
		{
			Qs[i] = Qs_all[i + thermalization_n];
		}

		// B autocorrelation
		autocorrelation(Qac, Qs);
		// Correlation factor calculation with fit
		double par[3];
		fit_to_exp(par, tac, Qac);
		// Number of dependent points
		double tau = (1. / par[1]);

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
			<< "T: " << Temp
			<< " rho: " << rho
			<< " therm_n: " << thermalization_n
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
