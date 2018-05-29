#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include <assert.h>

template<typename REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Model(Vector<number_type> x_data, Vector<number_type> y_data, Vector<number_type> dy_data) // default constructor
	: x_data_(x_data)
	, y_data_(y_data)
	, dy_data_(dy_data)
	{}

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 10;
	}

	/* degrees of freedom */
	size_type d_of_freedom() const
	{
		return x_data_.size() - n_parameters();
	}

	/* model function definition used for fitting */
	number_type f(number_type x, const Vector<number_type> & popt)
	{
		assert(n_parameters() == popt.size());
		return popt[0] * atan(popt[1]*x +popt[2]) + popt[3]+ popt[4]*exp(-(x-popt[5])*(x-popt[5])/2/popt[6]/popt[6]) + popt[7]*exp(-(x-popt[8])*(x-popt[8])/2/popt[9]/popt[9]);
	}

	/* definition of chi2 sum */
	number_type chi2(const Vector<number_type> & q)
	{
		/* checks if data arrays have equal lengths */
		assert(x_data_.size() == y_data_.size());
		assert(x_data_.size() == dy_data_.size());
		number_type chi2 = 0;
		for (size_type i = 0; i<x_data_.size(); ++i)
		{
			chi2 += (y_data_[i] - f(x_data_[i], q))*(y_data_[i] - f(x_data_[i], q))/dy_data_[i]/dy_data_[i];
		}

		return chi2;

	}

private:
	Vector<number_type> x_data_;
	Vector<number_type> y_data_;
	Vector<number_type> dy_data_;




};


#endif