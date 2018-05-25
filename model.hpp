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
		return 10;
	}

	/* model function definition used for fitting */
	number_type f(number_type x, const Vector<number_type> & popt)
	{
		assert(n_parameters() == popt.size());
		return popt[0] * atan(popt[1]*x +popt[2]) + popt[3]+ popt[4]*exp(-(x-popt[5])*(x-popt[5])/2/popt[6]/popt[6]) + popt[7]*exp(-(x-popt[8])*(x-popt[8])/2/popt[9]/popt[9]);
	}

private:


};


#endif