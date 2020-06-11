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
#include <random>

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

double Vtot(
	const std::vector<Vec3D> &pos,
	const std::function<double(double)> &Vr,
	const double L)
{
	// Number of particles
	size_t N{pos.size()};

	double d;
	Vec3D alias;

	// Cumulative value to return
	double Vreturn{0};

	for (size_t i{0}; i < N - 1; i++)
	{
		// Acceleration on i
		// caused by all other particles
		for (size_t j{i + 1}; j < N; j++)
		{
			// Copy the particle j
			alias = pos[j];
			// Check if there is interaction
			// and also compute the correct alias
			if ((d = compute_alias(pos[i], alias, L)) > 0)
			{
				// add potential
				Vreturn += Vr(d);
			}
		}
	}
	return Vreturn;
}

int main(int argc, char const *argv[])
{
	// Number of cells on side
	constexpr uint32_t n = 5;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n;
	// Number of evo steps
	constexpr uint32_t M = 20000;
	// Indexes sections
	// How much to wait for thermalization from the initial regular lattice
	constexpr uint32_t i_thermalization = 2000;
	// Length of simulation
	constexpr uint32_t i_simulation = M - i_thermalization;
	// When to start compute the correlations
	constexpr uint32_t i_correlation_start = i_thermalization;
	// Total length of the correlation
	constexpr uint32_t i_correlation_length = 1500;
	// Length of the array that saves the data points
	constexpr uint32_t i_correlation_corrW0_length = M - i_correlation_start;
	// How many means to take prior to apply a correcton to Delta
	constexpr uint32_t i_Delta_correction_means = 100;
	// Start and stop of the period of corrections
	constexpr uint32_t i_Delta_correction_start = i_Delta_correction_means;
	constexpr uint32_t i_Delta_correction_stop = i_thermalization;

	// Physics dimensional quantities
	// Univ constants
	constexpr double kB_r = 1.38e-23; // J*K^-1
	// Problem constants
	constexpr double T_r = 119.76;								   // K
	constexpr double sigma_r = 3.405e-10;						   // m
	constexpr double eps_r = kB_r * T_r;						   // J
	constexpr double rhom_r = 0.001 / sigma_r / sigma_r / sigma_r; // #m^-3
	constexpr double rhoM_r = 1.2 / sigma_r / sigma_r / sigma_r;   // #m^-3

	/* Conversion
	*	Kl=sigma_r
	*	Krho=sigma_r^-3
	*	Ke=eps_r
	*	KT=T_r
	*	KF=Ke/Kl
	*	KP=KF/Kl^2
	*/

	// Problem adimensional quantities
	constexpr double T{1}; // Temperature
	constexpr double beta{1};

	// Print info data
	std::cout
		<< "M:\t" << M << '\n'
		<< "N:\t" << N << '\n'
		<< "T:\t" << T << '\n'
		<< "KP:" << eps_r / sigma_r / sigma_r / sigma_r << "Pa" << '\n'
		<< "KV:" << sigma_r * sigma_r * sigma_r << "m^3" << '\n'
		<< std::endl;

	// Used in the algorithm
	// Read the rho
	double rho;
	std::cout << "Rho: ";
	std::cin >> rho;
	double V;
	double Delta;
	double L;
	double A;
	double nu;
	double Vp0, Vp1;
	double W, W2, DW;
	double F;
	double d;
	double P, stdP;
	Vec3D alias;
	std::vector<double> As(i_Delta_correction_means);

	// Positions statuses
	std::vector<Vec3D> pos_a(N);
	std::vector<Vec3D> pos_b(N);

	std::vector<Vec3D> *pos0 = &pos_a;
	std::vector<Vec3D> *pos1 = &pos_b;
	std::vector<Vec3D> *tmp;

	// Initialize the rho for which i want the P
	// std::vector<double> rho(Nrho);
	// for (uint32_t i; i < Nrho; i++)
	// {
	// 	rho[i] = rhom * std::exp(std::log(rhoM / rhom) * i / (Nrho - 1));
	// 	// rho[i] = rhom + (rhoM / rhom) * i / (Nrho - 1);
	// }

	// Istantiate the random generator
	// std::random_device rd;
	// std::mt19937 gen(6789);
	// std::uniform_real_distribution<double> dis(0.0, 1.0);
	// const auto randu = [&gen, &dis]() -> double {
	const auto randu = []() -> double {
		constexpr double normalization =
			1. / (static_cast<double>(RAND_MAX) + 1);
		return rand() * normalization;
		// return dis(gen);
	};

	// Autocorrelation on virial
	std::vector<double> corrW(i_correlation_length);
	fill(corrW, 0);
	std::vector<double> corrW0(i_correlation_corrW0_length);
	fill(corrW0, 0);

	// Setup output file
	std::ofstream PV_file("output_data/PV.dat");
	PV_file << "#V P std[P] rho W DW Delta\n";

	// Save the correlations results
	std::ofstream cW_file("output_data/cVV.dat");
	cW_file << "# i F\n";

	std::cout << "Start" << std::endl;

	// Volume of the box
	V = N / rho;
	// Side length
	L = std::pow(V, 1.0 / 3);

	// Now start sampling P given V
	// really we only want the quantity W which is needed for P

	std::cout << n << "\tVol:" << V << std::endl;
	// Delta doesn't need to be constant
	// Ideal value is given by the fraction of accepted states
	// that should be 30-70%
	// I try a bounch of Delta, compute tau of correlation
	// and then use the Delta with tau minimum
	// D small -> p=1
	// D big -> p=0
	std::vector<double> del(5);
	// linspace(del, 2 * L / n / 23, 2 * L / n / 26);
	linspace(del, 2 * L / n / 100, 2 * L / n / 300);
	for (uint32_t de = 0; de < del.size(); de++)
	{
		Delta = del[de];
		std::cout << "de:" << de << " Delta: " << Delta << std::endl;

		// Initial uniform position
		init_lattice(*pos0, L, n, 1);
		apply_periodic_bounds(*pos0, L);

		// Get the corresponding potential
		Vp0 = Vtot(*pos0, V_LJ, L);

		// null the W to start adding the samples and later take the mean
		W = W2 = 0;

		// TODO: watch for equilibrium
		// Metropolis algorithm for generation of new samples of my system
		for (uint32_t i = 0; i < M; i++)
		{
			// Output some information on progress
			if (i % 1000 == 0)
			{
				std::cout << "Progress: " << i * 100 / M << std::endl;
			}
			// New state
			for (uint32_t j = 0; j < N; j++)
			{
				(*pos1)[j].x = (*pos0)[j].x + Delta * (randu() - 0.5);
				(*pos1)[j].y = (*pos0)[j].y + Delta * (randu() - 0.5);
				(*pos1)[j].z = (*pos0)[j].z + Delta * (randu() - 0.5);
			}
			apply_periodic_bounds(*pos1, L);

			// Compute this potential
			Vp1 = Vtot(*pos1, V_LJ, L);
			// A = std::exp(-beta * Vp1)/std::exp(-beta * Vp0);
			A = std::exp(-beta * (Vp1 - Vp0));
			// Not tecnically needed
			A = std::min(1.0, A);

			// Adjust Delta value
			// A=0.5 is the target
			// Save the value cyclic-ally
			As[i % i_Delta_correction_means] = A;
			// Check if in Correction of Delta window
			if (i >= i_Delta_correction_start &&
				i < i_Delta_correction_stop &&
				i % i_Delta_correction_means == 0)
			{
				double mean_As = mean(As);
				// Apply correction given the mean of A on the last elements
				Delta *= 1 + 0.01 * (mean_As - 0.3);
				std::cout << "Mean prob:\t" << mean_As << "\n";
			}
			else if (i == i_Delta_correction_stop)
			{
				// The resulting Delta after the corrections
				std::cout << "Simulation Delta:" << Delta
						  << " A mean: " << mean(As)
						  << std::endl;
			}

			// Decide if accept the new status
			nu = randu();
			if (nu < A)
			{
				// accept the new configuration
				tmp = pos0;
				pos0 = pos1;
				pos1 = tmp;

				Vp0 = Vp1;
			}
			// else the new sample is the same as before

			// *pos0 is now the new sample to calculate the observables on

			// Calculate the observables wanted only after thermalization
			if (i >= i_thermalization)
			{
				for (uint32_t ii{0}; ii < N; ii++)
				{
					for (uint32_t jj{0}; jj < N; jj++)
					{
						if (jj != ii)
						{
							// Copy the particle j
							alias = (*pos0)[jj];
							// Check if there is interaction
							// and also compute the correct alias
							if ((d = compute_alias((*pos0)[ii], alias, L)) > 0)
							{
								// They do interact
								F = 0.5 * dV_LJ(d) * d;
								W += F;
								W2 += F * F;

								// Save for correlations on F
								if (i >= i_correlation_start)
								{
									corrW0[i - i_correlation_start] = F;
								}
							}
						}
					}
				}
			}
		}

		// Compute the correlation
		autocorrelation(corrW, corrW0);

		// Save correlation for plot and finding tau
		cW_file << "\"Delta: " << Delta << "\"\n";
		for (uint32_t i = 0; i < i_correlation_length; i++)
		{
			cW_file << i << ' ' << corrW[i] << '\n';
		}
		cW_file << "\n\n";

		// Do the mean
		W /= i_simulation;
		W2 /= i_simulation;
		// Not ok but leave it here
		DW = 25 * 1.0 / (i_simulation - 1) * std::fabs(W2 - W * W);

		// Compute the sample of P at volume V
		P = N / V - 1.0 / 3 / V * W;
		stdP = std::sqrt(std::pow(1.0 / 3 / V, 2) * DW);
		// Write the new obtained P
		PV_file << V << ' '
				<< P << ' '
				<< stdP << ' '
				<< rho << ' '
				<< W << ' '
				<< DW << ' '
				<< Delta << ' '
				<< "\n";

	} //den
	cW_file.close();
	PV_file.close();

	std::cout << "End" << std::endl;
	return 0;
}
