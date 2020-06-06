#include "simulation.h"

int init_lattice(std::vector<Vec3D> &pos, double L, uint32_t n, uint32_t q)
{
	// Get the number of particles
	size_t N{pos.size()};
	// Length of a single cell
	double dl{L / n};
	if (q * n * n * n > N)
	{
		// too few particles spaces
		return INIT_LATTICE_N_TOO_SMALL;
	}
	if (!(q == 1))
	{
		// Only q=1 is supportes as of now
		return INIT_LATTICE_Q_NOT_SUPPORTED;
	}

	// Place the particles
	for (uint32_t i{0}; i < n; i++)
	{
		for (uint32_t j{0}; j < n; j++)
		{
			for (uint32_t k{0}; k < n; k++)
			{
				if (q == 1)
				{
					pos[q * (i * n * n + j * n + k)] =
						Vec3D{i * dl, j * dl, k * dl};
				}
			}
		}
	}
	return 0;
}

void init_distribute_maxwell_boltzmann(Vec3D &vec, double std)
{
	// Give a normal distribution to the coponents
	vec.x = std * randn();
	vec.y = std * randn();
	vec.z = std * randn();
}

double new_acceleration(
	std::vector<Vec3D> &acc,
	const std::vector<Vec3D> &pos,
	const std::function<double(double)> &Vr,
	const double m,
	const double L,
	const std::function<double(double)> &dVr_arg)
{
	// Needed for evenual 5 points differentiation
	constexpr double diff_eps{1e-10};

	// Number of particles
	size_t N{pos.size()};

	double d;
	Vec3D tmp_acc;
	Vec3D alias;
	std::function<double(double)> dVr;
	// Cumulative value to return
	double Vreturn{0};

	// Check if derivative has been supplied
	if (dVr_arg)
	{
		dVr = dVr_arg;
	}
	else
	{
		// Computational derivative
		dVr = [&](double d) -> double {
			return derive_5points(Vr, d, diff_eps);
		};
	}

	for (size_t i{0}; i < N; i++)
	{
		// Clear the acceleration before start adding contributions
		acc[i].clear();

		// Acceleration on i
		// caused by all other particles
		for (size_t j{0}; j < N; j++)
		{
			// Skip the self case
			if (i != j)
			{
				// Copy the particle j
				alias = pos[j];

				// Check if there is interaction
				// and also compute the correct alias
				if ((d = compute_alias(pos[i], alias, L)) > 0)
				{
					// Compute the derivative on function v around d
					// Compute the acc contribute for the pair ij
					// Add the contribute to the acceleration
					// using the difference versor
					acc[i] += -dVr(d) / m / d * (pos[i] - alias);
					// add potential on half of the cases
					if (i < j)
						Vreturn += Vr(d);
				}
				// Else: no interaction, un-perfect packing of spheres
			}
		}
	}
	return Vreturn;
}

int sign(double x)
{
	return x > 0 ? 1 : -1;
}
double compute_alias(const Vec3D &me, Vec3D &other, double L)
{
	if (std::fabs(me.x - other.x) > L / 2)
	{
		other.x += sign(me.x - other.x) * L;
	}
	if (std::fabs(me.y - other.y) > L / 2)
	{
		other.y += sign(me.y - other.y) * L;
	}
	if (std::fabs(me.z - other.z) > L / 2)
	{
		other.z += sign(me.z - other.z) * L;
	}

	double d{(me - other).norm()};
	// Return the distance
	// but a negative value if is >= L/2
	return d < L / 2 ? d : -d;

	// Old brute-force version
	// Vec3D alias;
	// double d;
	// // Iterate on all possible alias of j
	// // it will interact with at most one of the aliases
	// for (uint8_t a{0}; a < 27; a++)
	// {
	// 	alias = other + L * aliaser[a];
	// 	d = (me - alias).norm();
	// 	// Check interaction
	// 	if (d < L / 2)
	// 	{
	// 		// Now that i have found the alias
	// 		// i can break iterating on them
	// 		other = alias;
	// 		return d;
	// 	}
	// 	// else ignore interaction
	// }
	// return -1;
}

void apply_periodic_bounds(std::vector<Vec3D> &pos, double L)
{
	size_t N{pos.size()};
	for (size_t i{0}; i < N; i++)
	{
		pos[i].x -= L * std::floor(pos[i].x / L);
		pos[i].y -= L * std::floor(pos[i].y / L);
		pos[i].z -= L * std::floor(pos[i].z / L);
	}
}

const std::vector<Vec3D> aliaser{
	// in order of feasability
	Vec3D{0, 0, 0},
	Vec3D{0, 0, -1},
	Vec3D{0, 0, 1},
	Vec3D{0, -1, 0},
	Vec3D{0, 1, 0},
	Vec3D{-1, 0, 0},
	Vec3D{1, 0, 0},

	Vec3D{-1, -1, 0},
	Vec3D{1, 1, 0},
	Vec3D{-1, 0, -1},
	Vec3D{1, 0, 1},
	Vec3D{-1, 0, 1},
	Vec3D{-1, 1, 0},
	Vec3D{0, -1, -1},
	Vec3D{0, -1, 1},
	Vec3D{0, 1, -1},
	Vec3D{0, 1, 1},
	Vec3D{1, -1, 0},
	Vec3D{1, 0, -1},

	Vec3D{-1, -1, -1},
	Vec3D{-1, -1, 1},
	Vec3D{-1, 1, -1},
	Vec3D{-1, 1, 1},
	Vec3D{1, -1, -1},
	Vec3D{1, -1, 1},
	Vec3D{1, 1, -1},
	Vec3D{1, 1, 1},
};

// AutoCorr constructor
AutoCorr::AutoCorr(
	uint32_t c_istart,
	uint32_t c_istop,
	uint32_t c_length,
	uint32_t Nq)
	: c_istart{c_istart},
	  c_istop{c_istop},
	  c_length{c_length},
	  Nq{Nq},
	  c_means_number{c_istop - c_istart - c_length + 1},

	  c(c_length),
	  c_0(Nq * c_length),
	  c_02m(c_length)
{
	// Clear all the initial conditions
	for (uint32_t i{0}; i < Nq * c_length; i++)
	{
		c_0[i].clear();
	}
	// And also the results vector
	fill(c, 0);
}

// Handles the single addition of new data
int AutoCorr::add_data(
	const std::vector<Vec3D> &q,
	uint32_t i)
{
	// Check temporal position in simulation
	if (i >= c_istart && i < c_istop)
	{
		// wrapping index needed for cvv_v0 and cvv_v02m
		// it shouldn't really be necessary to subtract c_istart
		uint32_t c_i = (i - c_istart) % c_length;

		// Setup initial conditions of correlations
		// until the incomplete sequences at the end
		if (i < c_istop - c_length + 1)
		{
			// Normal operation
			// Calculate already the RMS
			c_02m[c_i] = 0;
			for (uint32_t j = 0; j < Nq; j++)
			{
				// Save the new c_0 offsetting location
				// using c_i and calculate the norm2
				c_0[c_i * Nq + j] = q[j];
				c_02m[c_i] += q[j].norm2();
			}

			// Mean
			// Can be avoided because cancels out
			// with the other mean at the numerator
			// c_02m[i] /= Nq;
		}
		else
		{
			// Last cvv_length steps
			// setting to zero will prevent
			// adding the corresponding incomplete contributes
			c_02m[c_i] = 0;
		}

		// Contribution on c calculations
		double c_tmp;
		uint32_t jjj;
		for (uint32_t jj = 0; jj < c_length; jj++)
		{
			// Contribution on c[jj]

			// The index of c_02m and c_0 corresponding to jj
			jjj = (c_i - jj + c_length) % c_length;

			// Skip the steps for which c_02m==0
			// they happens if I'm in the last steps for which
			// it has been hard set above
			// or can happen if the mean velocity of *ALL* the particles
			// is null either way it would skyrocket the correlation
			if (c_02m[jjj] != 0)
			{
				// Incremental mean
				c_tmp = 0;
				for (uint32_t j = 0; j < Nq; j++)
				{
					// Scalar product
					c_tmp += q[j] * c_0[jjj * Nq + j];
				}
				// Mean on particles
				// Avoided as stated before on computing c_02m,
				// decommenting *BOTH* lines would not change the result
				// c_tmp /= Nq;
				// Normalize on the <c0^2>
				c_tmp /= c_02m[jjj];
				// Add contribution
				c[jj] += c_tmp;
			}
		}
	}
	else if (i < c_istart)
	{
		return AUTCO_EARLY_DATA;
	}
	else if (i >= c_istop)
	{
		return AUTCO_LATE_DATA;
	}
	return 0;
}

// Getter for the data
double AutoCorr::operator[](size_t i) const
{
	// Means of the final correlation result
	return c[i] / c_means_number;
}

uint32_t AutoCorr::get_means_number() const
{
	return c_means_number;
}

void autocorrelation(std::vector<double> &corr,
					 const std::vector<double> &x)
{
	double m = mean(x);

	for (int t = 0; t < corr.size(); t++)
	{
		double n = 0; // Numerator
		double d = 0; // Denominator

		for (int i = 0; i < x.size() - t; i++)
		{
			double xim = x[i] - m;
			n += xim * (x[i + t] - m);
			d += xim * xim;
		}

		corr[t] = n / d;
	}
}
