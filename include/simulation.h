#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <cmath>

class Vec3D{
public:
	Vec3D(double x=0,double y=0,double z=0);

	double x;
	double y;
	double z;

	double norm() const;
	double norm2() const;

	void normalize();
	Vec3D Vec3D::getNormalized() const;

	friend Vec3D operator-(const Vec3D &a, const Vec3D &b);
};

#endif // __SIMULATION_H__
