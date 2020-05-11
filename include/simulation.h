#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <cstdint>
#include <cmath>
#include <vector>

class Vec3D
{
public:
	Vec3D(double x = 0, double y = 0, double z = 0);

	double x;
	double y;
	double z;

	double norm() const;
	double norm2() const;

	void normalize();
	Vec3D getNormalized() const;

	friend Vec3D operator-(const Vec3D &a, const Vec3D &b);
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
#define INIT_LATTICE_N_TOO_SMALL 3
int init_lattice(std::vector<Vec3D> &pos, double dl, uint32_t n = 1, uint32_t q = 1);

#endif // __SIMULATION_H__
