#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <assert.h>
#include "random.hpp"
#include <random>

template<typename REAL>
class Vector : public std::vector<REAL> // inherit from the STL vector
{
public:
	typedef std::size_t size_type; // type used for array indices
	Vector() : std::vector<REAL>() // default constructor, inherited
	{
	}

	/* constructor to create an array of given size and filled by a value */
	Vector( const size_type size, const REAL defaultvalue = 0) 
	: std::vector<REAL>(size, defaultvalue)
	{

	}

	/* fill vector with one value */
	Vector & operator= (const REAL value)
	{
		const size_type s = this->size();
		Vector & self = *this;
		for (size_type i = 0; i<s; ++i)
		{
			self[i] = value;
		}
		return *this;
	}

	/* Multiplication by a scalar */
	Vector & operator*= (const REAL value)
	{
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] *= value;
		}
		return *this;
	}

	/* Division by a scalar */
	Vector & operator/= (const REAL value)
	{
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] /= value;
		}
		return *this;
	}

	/* Add another vector */
	Vector & operator+= (const Vector & y)
	{
		assert(this->size() == y.size());
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] += y[i];
		}
		return *this;
	}
	/* Subtract another vector */
	Vector & operator-= (const Vector & y)
	{
		assert(this->size() == y.size());
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] -= y[i];
		}
		return *this;
	}

	/* Adding two vectors y+x */
	Vector operator+ (const Vector & x) const
	{
		assert(x.size() == this->size()); // checks if dimensions of the two vectors match
		Vector sum(*this);
		sum += x;
		return sum;
	}

	/* Adding a constant y+alpha (component-wise) */
	Vector operator+ (const REAL & x) const
	{
		Vector sum(*this);
		for (size_type i = 0; i < sum.size(); ++i)
		{
			sum[i] += x;
		}
		return sum;
	}

	/* Subtraction of two vectors y-x */
	Vector operator- (const Vector & x) const
	{
		assert(x.size() == this->size()); // checks if dimensions of the two vectors match
		Vector sum(*this);
		sum -= x;
		return sum;
	}

	/* Negation of a vector */
	Vector & operator-()
	{
		Vector & neg(*this);
		neg *= -1.0;
		return *this;
	}


	/* Square of the Euclidean norm */
	REAL two_norm_2 () const
	{
		REAL sum(0);
		const Vector & self = *this;
		for (size_type i = 0; i < (size_type) this->size(); ++i)
		{
			sum += self[i] * self[i];
		}
		return sum;
	}

};

// some free functions

/* Multiplication from the left hand side by scalar */
	template<typename REAL>
	Vector<REAL> operator* (const REAL & alpha, const Vector<REAL> x)
	{
		Vector<REAL> result(x);
		result *= alpha;
		return result;
	}
/* Multiplication from the right hand side by scalar */
	template<typename REAL>
	Vector<REAL> operator* (const Vector<REAL> x, const REAL & alpha)
	{
		Vector<REAL> result(x);
		result *= alpha;
		return result;
	}

/* Division from the right hand side by scalar */
	template<typename REAL1, typename REAL2>
	Vector<REAL1> operator/ (const Vector<REAL1> x, const REAL2 & alpha)
	{
		Vector<REAL1> result(x);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] = x[i]/alpha;
		}
		return result;
	}

/* Fill vector with random entries following a probability distribution */
	template<typename REAL, class Distribution>
	void fill_random(Vector<REAL> & vec, Distribution distribution)
	{
		for (std::size_t i = 0; i < vec.size(); ++i)
		{
			vec[i] = distribution.draw();
		}

	}

/* Fill vector with random entries from a specific region */
	template<typename REAL>
	void fill_from_region(Vector<REAL> &vec, Vector<REAL> region_min, Vector<REAL> region_max)
	{
		assert(vec.size() == region_min.size() && vec.size() == region_max.size());
		for (std::size_t i = 0; i < vec.size(); ++i)
		{
			Uniform_real_distribution<REAL> dis_unif(region_min[i], region_max[i]);
			vec[i] = dis_unif.draw();
		}

	}


#endif