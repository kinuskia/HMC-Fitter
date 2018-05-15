#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>

template<typename REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Model() // default constructor
	{}

	/* model function definition used for fitting */
	number_type f(number_type x, const Vector<number_type> & popt)
	{
		return popt[0] + popt[1]*sin(popt[2]*x);
	}


};


#endif