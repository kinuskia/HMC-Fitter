#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/read_data.hpp"
#include "../auxiliary_files/matrix.hpp"

#include <random>
#include <assert.h>
#include <fstream>

template <typename REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Model(size_type n_parameters)
	: n_parameters_(n_parameters)
	{}

	size_type n_parameters() const
	{
		return n_parameters_;
	}

	size_type d_of_freedom() const
	{
		return 1;
	}

	number_type potential(const Vector<number_type> & q)
	{
		number_type sum = 0;
		for (size_type i = 0; i < q.size(); ++i)
		{
			sum += (q[i] - number_type(i))*(q[i] - number_type(i))/(number_type(i)+1)/(number_type(i)+1)/2.;
		}
		return sum;
	}
	bool constraints_respected(const Vector<number_type> &q)
	{
		return true;
	}
private:
	size_type n_parameters_;



};

typedef double number_type;
typedef std::size_t size_type;

#include "../auxiliary_files/hmc.hpp"
#include <vector>
#include <ctime>


int main ()
{
	// size_type n_param = 12;
	// // Estimated search region
	// Vector<number_type> range_min(n_param);
	// Vector<number_type> range_max(n_param);
	// range_min = -1.;
	// range_max = 1.+ n_param;

	// // Characteristic length scales
	// Vector<number_type> c_lengths(n_param, 1);
	// c_lengths = range_max - range_min;

	// // Setting up model and sampler
	// Model<number_type> gaussian(n_param);
	// HMC<number_type> sampler(gaussian, range_min, range_max, c_lengths, 5e-3, 40, 50, 1e-1);

	// // find minimum
	// sampler.walk_automatic();


	

	return 0;
}