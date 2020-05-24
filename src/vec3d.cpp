#include "vec3d.h"

Vec3D::Vec3D(double x, double y, double z)
	: x{x},
	  y{y},
	  z{z}
{
}

double Vec3D::norm() const
{
	return std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

double Vec3D::norm2() const
{
	return this->x * this->x + this->y * this->y + this->z * this->z;
}

Vec3D &Vec3D::normalize()
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

Vec3D &Vec3D::clear()
{
	this->x = this->y = this->z = 0;
	return *this;
}

Vec3D &Vec3D::operator=(const Vec3D &a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
	return *(this);
}

Vec3D &Vec3D::operator+=(const Vec3D &a)
{
	this->x += a.x;
	this->y += a.y;
	this->z += a.z;
	return *(this);
}

Vec3D &Vec3D::operator-=(const Vec3D &a)
{
	this->x -= a.x;
	this->y -= a.y;
	this->z -= a.z;
	return *(this);
}

Vec3D &Vec3D::operator*=(const double k)
{
	this->x *= k;
	this->y *= k;
	this->z *= k;
	return *(this);
}

Vec3D &Vec3D::operator/=(const double k)
{
	this->x /= k;
	this->y /= k;
	this->z /= k;
	return *(this);
}

Vec3D operator-(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3D operator+(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3D operator*(const double k, const Vec3D &a)
{
	return Vec3D{a.x * k, a.y * k, a.z * k};
}
Vec3D operator*(const Vec3D &a, const double k)
{
	return k * a;
}

double operator*(const Vec3D &a, const Vec3D &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

#if __cplusplus >= 202002L
double operator<=>(const Vec3D &a, const Vec3D &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
#endif

std::ostream &operator<<(std::ostream &out, const Vec3D &vec)
{
	out << "[" << vec.x << "," << vec.y << "," << vec.z << "]";
	return out;
}
