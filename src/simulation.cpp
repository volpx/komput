#include "simulation.h"

Vec3D::Vec3D(double x,double y,double z):
	x{x}, y{y}, z{z}
{

}

double Vec3D::norm() const
{
	return std::sqrt(this->x*this->x+this->y*this->y+this->y*this->y);
}

double Vec3D::norm2() const
{
	return this->x*this->x+this->y*this->y+this->y*this->y;
}

void Vec3D::normalize()
{
	double n{this->norm()};
	this->x/=n;
	this->y/=n;
	this->z/=n;
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
	return Vec3D{a.x-b.x, a.y-b.y, a.y-b.y};
}
