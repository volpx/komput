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

void Vec3D::normalize()
{
	double n{this->norm()};
	this->x /= n;
	this->y /= n;
	this->z /= n;
}

Vec3D Vec3D::getNormalized() const
{
	// Use copy constructor to create a copy instance
	Vec3D res{*this};
	res.normalize();
	return res;
}

Vec3D operator-(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x - b.x, a.y - b.y, a.y - b.y};
}

Vec3D operator+(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x + b.x, a.y + b.y, a.y + b.y};
}

int init_born_von_karman_3d(double L, std::vector<Vec3D> q)
{
	size_t N = q.size();
	double dl = L / std::pow(N, 1.0 / 3);

	for (size_t i = 0;; i++)
	{
		for (size_t j = 0;; j++)
		{
			for (size_t k = 0;; k++)
			{
			}
		}
	}
}
