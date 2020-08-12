#pragma once

#include <cmath>
#include <cstdint>
#include <ostream>

template <typename T = float,
		  uint8_t d = 3>
struct VecD
{
	T q[d];

	// Empty constructor
	VecD(){};
	// Normal constructor (could overload the above)
	VecD(T x[d])
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] = x[i];
		}
	};
	// Copying constructor
	VecD(const VecD<T, d> &other);

	// Getters for the norm and norm^2
	template <typename F>
	F norm2() const
	{
		F sum{0};
		for (uint8_t i{0}; i < d; ++i)
		{
			sum += q[i] * q[i];
		}
		return sum;
	}
	template <typename F>
	F norm() const
	{
		return std::sqrt(norm2<F>());
	}

	T norm2() const
	{
		return norm2<T>();
	};
	T norm() const
	{
		return norm2();
	};

	// Make a versor by either modify myself or create a copy
	VecD<T, d> &normalize()
	{
		const T n{norm()};
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] /= n;
		}
		return *this;
	};
	VecD<T, d> getNormalized() const
	{
		VecD<T, d> ret{*this};
		ret.normalise();
		return ret;
	};

	// Set to null vector
	VecD<T, d> &clear()
	{
		(*this) *= 0;
		return *this;
	};

	// Overload member access operator
	T &operator[](uint8_t i)
	{
		return q[i];
	};
	const T &operator[](uint8_t i) const
	{
		return q[i];
	};

	// Overload of arithmetic operators
	VecD<T, d> &operator=(const VecD<T, d> &a)
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] = a[i];
		}
		return *this;
	};

	VecD<T, d> &operator+=(const VecD<T, d> &a)
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] += a[i];
		}
		return *this;
	};
	VecD<T, d> &operator-=(const VecD<T, d> &a)
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] -= a[i];
		}
		return *this;
	};
	VecD<T, d> &operator*=(const T k)
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] *= k;
		}
		return *this;
	};
	VecD<T, d> &operator/=(const T k)
	{
		for (uint8_t i{0}; i < d; ++i)
		{
			q[i] /= k;
		}
		return *this;
	};
};

template <typename T,
		  uint8_t d>
VecD<T, d> operator+(const VecD<T, d> &a, const VecD<T, d> &b)
{
	VecD<T, d> ret{};
	for (uint8_t i{0}; i < d; ++i)
	{
		ret[i] = a[i] + b[i];
	}
	return ret;
}
template <typename T,
		  uint8_t d>
VecD<T, d> operator-(const VecD<T, d> &a, const VecD<T, d> &b)
{
	VecD<T, d> ret{};
	for (uint8_t i{0}; i < d; ++i)
	{
		ret[i] = a[i] - b[i];
	}
	return ret;
}
template <typename T,
		  uint8_t d,
		  typename K>
VecD<T, d> operator*(const K k, const VecD<T, d> &a)
{
	VecD<T, d> ret{};
	for (uint8_t i{0}; i < d; ++i)
	{
		ret[i] = k * a[i];
	}
	return ret;
}

template <typename T,
		  uint8_t d>
std::ostream &operator<<(std::ostream &out, const VecD<T, d> &vec)
{
	out << "Vec" << (int)d << "D[";
	for (uint8_t i{0}; i < d; ++i)
	{
		out << vec[i];
	}
	out << ']';
	return out;
}
