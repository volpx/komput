#include "simulation.h"

Vec3D::Vec3D(double x, double y, double z) : x{x}, y{y}, z{z}
{
}

double Vec3D::norm() const
{
	return std::sqrt(this->x * this->x + this->y * this->y + this->y * this->y);
}

double Vec3D::norm2() const
{
	return this->x * this->x + this->y * this->y + this->y * this->y;
}

Vec3D Vec3D::normalize()
{
	double n{this->norm()};
	this->x /= n;
	this->y /= n;
	this->z /= n;

	return *(this);
}

Vec3D Vec3D::getNormalized() const
{
	// Use copy constructor to create a copy instance
	Vec3D res{*this};
	res.normalize();
	return res;
}

Vec3D Vec3D::operator+=(const Vec3D &a)
{
	this->x += a.x;
	this->y += a.y;
	this->z += a.z;
	return *(this);
}

Vec3D Vec3D::operator-=(const Vec3D &a)
{
	this->x -= a.x;
	this->y -= a.y;
	this->z -= a.z;
	return *(this);
}

Vec3D operator-(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x - b.x, a.y - b.y, a.y - b.y};
}

Vec3D operator+(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x + b.x, a.y + b.y, a.y + b.y};
}

Vec3D operator*(const double k, const Vec3D &a)
{
	return Vec3D{a.x * k, a.y * k, a.y * k};
}
Vec3D operator*(const Vec3D &a, const double k)
{
	return k * a;
}

int init_lattice(std::vector<Vec3D> &pos, double dl, uint32_t n, uint32_t q)
{
	size_t N{pos.size()};
	if (q * n * n * n > N)
	{
		// too few particles spaces
		return INIT_LATTICE_N_TOO_SMALL;
	}
	if (!(q == 1))
	{
		// q not supported
		return INIT_LATTICE_Q_NOT_SUPPORTED;
	}

	for (uint32_t i{0}; i < n; i++)
	{
		printf("I:%d", i);
		for (uint32_t j{0}; j < n; j++)
		{
			for (uint32_t k{0}; k < n; k++)
			{
				if (q == 1)
				{
					// printf("I:%d,J:%d,K:%d\n", i, j, k);
					pos[q * (k + j * n + i * n * n)] = Vec3D{0, 0, 0} + Vec3D{i * dl, j * dl, k * dl};
				}
			}
		}
	}
	return 0;
}

void init_distribute_maxwell_boltzmann(Vec3D &vec, double std)
{
	double q, w, e, r;
	do
	{
		q = rand() * (1.0 / RAND_MAX);
	} while (q == 0);
	do
	{
		e = rand() * (1.0 / RAND_MAX);
	} while (e == 0);

	w = rand() * (1.0 / RAND_MAX);
	r = rand() * (1.0 / RAND_MAX);

	vec.x = std * std::sqrt(-2.0 * std::log(q)) * std::cos(2 * M_PI * w);
	vec.y = std * std::sqrt(-2.0 * std::log(q)) * std::sin(2 * M_PI * w);
	vec.z = std * std::sqrt(-2.0 * std::log(e)) * std::cos(2 * M_PI * r);
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

void apply_periodic_bounds_unitary(std::vector<Vec3D> &pos)
{
	size_t N{pos.size()};
	for (size_t i{0}; i < N; i++)
	{
		pos[i].x -= std::floor(pos[i].x);
		pos[i].y -= std::floor(pos[i].y);
		pos[i].z -= std::floor(pos[i].z);
	}
}
