#include "simulation.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

double V_LJ(const Vec3D &a, const Vec3D &b)
{
	double d2 = (a - b).norm2();
	return 4 * (std::pow(1 / d2, 6) - std::pow(1 / d2, 3));
}

double V(std::vector<Vec3D> &pos)
{
	double res{0};
	uint32_t l{pos.size()};
	// TODO: check ranges
	for (uint32_t i{0}; i < l - 1; i++)
	{
		for (uint32_t j{i + 1}; j < l; j++)
		{
			res += V_LJ(pos[i], pos[j]);
		}
	}
	return res;
}

int main(int argc, char const *argv[])
{
	int N = 1000;
	std::vector<Vec3D> pos(N);
	std::vector<Vec3D> vel(N);

	// TODO: Initialize with Born-Von Karman conditions

	// Pause
	getchar();

	return 0;
}
