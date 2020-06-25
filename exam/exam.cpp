#include "vec3d.h"
#include "simulation.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
extern "C"
{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
}

#define LJ_CORRECTION_OFFSET -(-0.00691026)

// Interaction potential
double V_LJ(double x)
{
#ifdef LJ_CORRECTION_OFFSET
	return 4 * (std::pow(1. / x, 12) - std::pow(1. / x, 6)) + LJ_CORRECTION_OFFSET;
#else
	return 4 * (std::pow(1. / x, 12) - std::pow(1. / x, 6));
#endif
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

int fit_to_exp(double lambda[3],
			   std::vector<double> &datax,
			   std::vector<double> &datay)
{
	const size_t n{datax.size()};

	struct data
	{
		size_t n;
		double *x;
		double *y;
	};

	auto exp_f{
		[](const gsl_vector *lambda, void *data,
		   gsl_vector *f) -> int {
			size_t n = static_cast<struct data *>(data)->n;
			double *x = static_cast<struct data *>(data)->x;
			double *y = static_cast<struct data *>(data)->y;
			double l0 = gsl_vector_get(lambda, 0);
			double l1 = gsl_vector_get(lambda, 1);
			double l2 = gsl_vector_get(lambda, 2);
			for (size_t i = 0; i < n; i++)
			{
				/* Model Yi = l0 * exp(-l1 * t_i) + l2 */
				double Yi = l0 * std::exp(-l1 * x[i]) + l2;
				gsl_vector_set(f, i, Yi - y[i]);
			}
			return GSL_SUCCESS;
		}};

	auto exp_df{
		[](const gsl_vector *lambda, void *data,
		   gsl_matrix *J) -> int {
			size_t n = static_cast<struct data *>(data)->n;
			double *x = static_cast<struct data *>(data)->x;
			double l0 = gsl_vector_get(lambda, 0);
			double l1 = gsl_vector_get(lambda, 1);
			for (size_t i = 0; i < n; i++)
			{
				/* Jacobian matrix J(i,j) = dfi / dxj, */
				/* where fi = (Yi - yi)/sigma[i],*/
				/* Yi = A * exp(-lambda * t_i) + b */
				/* and the xj are the parameters (A,lambda,b) */
				double e = std::exp(-l1 * x[i]);
				gsl_matrix_set(J, i, 0, e);
				gsl_matrix_set(J, i, 1, -x[i] * l0 * e);
				gsl_matrix_set(J, i, 2, 1.0);
			}

			return GSL_SUCCESS;
		}};

	// number of parameters
	const size_t p = 3;
	const gsl_multifit_nlinear_type *multifittype = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters fdf_params =
		gsl_multifit_nlinear_default_parameters();
	fdf_params.trs = gsl_multifit_nlinear_trs_lm;
	/* allocate workspace with default parameters */
	gsl_multifit_nlinear_workspace *w =
		gsl_multifit_nlinear_alloc(multifittype, &fdf_params, n, p);
	if (w == NULL)
	{
		std::cerr << "gsl_multifit_nlinear_alloc: " << w << std::endl;
		return -1;
	}
	double lambda_init[p] = {1.0, 10., 0.0};
	gsl_vector_view lambda_init_wiew = gsl_vector_view_array(lambda_init, p);
	gsl_multifit_nlinear_fdf fdf;
	// define the function to be minimized
	fdf.f = exp_f;
	// set to NULL for finite-difference Jacobian
	fdf.df = exp_df;
	// not using geodesic acceleration
	fdf.fvv = NULL;
	fdf.n = n;
	fdf.p = p;
	struct data d = {n, datax.data(), datay.data()};
	fdf.params = &d;
	gsl_multifit_nlinear_init(
		&lambda_init_wiew.vector,
		&fdf,
		w);

	int status, info;

	/* solve the system with a maximum of 100 iterations */
	const double xtol = 1e-8;
	const double gtol = 1e-8;
	const double ftol = 0.0;
	status = gsl_multifit_nlinear_driver(
		100, xtol, gtol, ftol,
		NULL, NULL, &info, w);

	// std::cout
	// 	<< "status: " << status << '\n'
	// 	<< "niter:" << gsl_multifit_nlinear_niter(w) << '\n'
	// 	<< "reason stop: "
	// 	<< ((info == 1) ? "small step size" : "small gradient") << '\n'
	// 	<< "l0: " << gsl_vector_get(w->x, 0) << '\n'
	// 	<< "l1: " << gsl_vector_get(w->x, 1) << '\n'
	// 	<< "l2: " << gsl_vector_get(w->x, 2) << '\n'
	// 	<< std::endl;
	lambda[0] = gsl_vector_get(w->x, 0);
	lambda[1] = gsl_vector_get(w->x, 1);
	lambda[2] = gsl_vector_get(w->x, 2);

	gsl_multifit_nlinear_free(w);

	return 0;
}

int main()
{
	// Number of cells on side
	constexpr uint32_t n = 5;
	// Number of particles in the box
	constexpr uint32_t N = n * n * n;
	// Number of time steps
	constexpr uint32_t M = 10000;

	// Physics constants
	constexpr double kB_r = 1.3806503e-23;	  // J
	constexpr double hbar_r = 1.05457148e-34; // J s

	// Argon data
	constexpr double sigma_r = 3.398e-10;			// m
	constexpr double T_r = 122;						// K
	constexpr double eps_r = T_r * kB_r;			// J
	constexpr double m_r = 39.95 * 1e-3 / 6.022e23; // Kg
	constexpr double dt_r = 1e-14;					// s
	// vary rho from 0.65 to 0.9
	std::vector<double> rhos(10);
	linspace(rhos, 0.65, 0.9);
	const double rho = rhos[0];		   //
	constexpr double lambda = 0.86e-3; //

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
#ifdef LJ_CORRECTION_OFFSET
		<< "\nV_LJ off:\t" << LJ_CORRECTION_OFFSET
#endif
		<< "\nV_LJ(L/2):\t" << V_LJ(L / 2)
		<< "\nvstd:\t" << vstd
		<< '\n'
		<< std::endl;

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

	// Wanted quantities
	std::vector<double> Qs(M);

	// Boilerplate variables
	struct IntCall
	{
		std::vector<Vec3D> **pos = nullptr;
		std::vector<Vec3D> **acc = nullptr;
		double (*Vr)(double) = nullptr;
		double (*dVr)(double) = nullptr;
		double (*ddVr)(double) = nullptr;
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
				V += Vr(d);
			// Compute the quantity
			Q += ddVr(d) + dVr(d) * 2 / d;
		}
	} interaction;
	interaction.pos = &pos1;
	interaction.acc = &acc1;
	interaction.Vr = V_LJ;
	interaction.dVr = dV_LJ;
	interaction.ddVr = ddV_LJ;

	struct IntSetup
	{
		std::vector<Vec3D> **acc = nullptr;
		void operator()(const size_t i)
		{
			// Clear the acceleration before start adding contributions
			(**acc)[i].clear();
		}
	} start_interaction;
	start_interaction.acc = &acc1;

	double E, T;
	double Qsum = 0, Q2sum = 0;
	double &V = interaction.V;
	double &Q = interaction.Q;
	std::ofstream data_file("output_data/data_evolution.dat");
	data_file << "#i,t,T,V,E,Q\n";

	//******************************* STARTING *********************************

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
	data_file
		<< 0 << ' '
		<< 0 << ' '
		<< T << ' '
		<< V << ' '
		<< E << ' '
		<< Q << '\n';

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
		data_file
			<< i << ' '
			<< i * dt << ' '
			<< T << ' '
			<< V << ' '
			<< E << ' '
			<< Q << '\n';

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
	data_file.close();

	// Post analysis
	std::cout << "Post analysis" << std::endl;

	// B autocorrelation
	const uint32_t i_ac_length = 500;
	std::vector<double> Qac(i_ac_length);
	std::vector<double> tac(i_ac_length);
	autocorrelation(Qac, Qs);
	std::ofstream Qac_file("output_data/Qautocorr.dat");
	Qac_file << "#i,t,Qac\n";
	for (uint32_t i = 0; i < Qac.size(); ++i)
	{
		tac[i] = i * dt;
		Qac_file << i << ' ' << tac[i] << ' ' << Qac[i] << '\n';
	}
	Qac_file.close();

	// Correlation factor calculation with fit
	double par[3];
	fit_to_exp(par, tac, Qac);
	double tau = 1. / (par[1]);

	double meanQ = mean(Qs);
	double DmeanQ = tau / (M - 1) * variance(Qs);

	std::cout
		<< "B: " << meanQ / 48
		<< " +/- " << std::sqrt(DmeanQ) / 48
		<< std::endl;

	return 0;
}
