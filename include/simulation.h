#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <cstdint>
#include <cmath>
#include <vector>
#include <cstdlib>

//debug only
#include <iostream>

class Vec3D
{
public:
	Vec3D(double x = 0, double y = 0, double z = 0);

	double x;
	double y;
	double z;

	double norm() const;
	double norm2() const;

	Vec3D normalize();
	Vec3D getNormalized() const;

	Vec3D operator+=(const Vec3D &a);
	Vec3D operator-=(const Vec3D &a);

	friend Vec3D operator-(const Vec3D &a, const Vec3D &b);
	friend Vec3D operator+(const Vec3D &a, const Vec3D &b);
	friend Vec3D operator*(const Vec3D &a, const double k);
	friend Vec3D operator*(const double k, const Vec3D &a);
};

/* Init the cubic lattice positions with:
 *  dl
 *  n cells along the side
 *  q particles for every cell
 *
 * Return:
 *  0 on succesful init
 *
 * NOTE:
 *  q determines the configuration used
 */
#define INIT_LATTICE_N_TOO_SMALL -3
#define INIT_LATTICE_Q_NOT_SUPPORTED -4
int init_lattice(std::vector<Vec3D> &pos, double dl, uint32_t n = 1, uint32_t q = 1);

void init_distribute_maxwell_boltzmann(Vec3D &vec, double std = 1);

void apply_periodic_bounds(std::vector<Vec3D> &pos, double L = 1);
void apply_periodic_bounds_unitary(std::vector<Vec3D> &pos);

#endif // __SIMULATION_H__
