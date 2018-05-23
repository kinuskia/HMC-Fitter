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
	Model() // default constructor
	{}

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 5;
	}

	/* model function definition used for fitting */
	number_type f(number_type x, const Vector<number_type> & popt)
	{
		assert(n_parameters() == popt.size());
		return popt[0] + popt[1]*exp(-popt[2]*x)*sin(popt[3]*x + popt[4]);
	}

private:


};


#endif