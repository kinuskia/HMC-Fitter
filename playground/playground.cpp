#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/storage.hpp"
//#include "../auxiliary_files/read_data.hpp"
//#include "../auxiliary_files/matrix.hpp"

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

//#include "../auxiliary_files/hmc.hpp"

#include <vector>
#include <ctime>


int main ()
{
	size_type N = 1e5;
	Vector<number_type> popt(3);
	Storage<number_type> data(popt);

	size_type counter = 0;
	for (size_type i = 0; i < N; ++i)
	{
		for (size_type j = 0; j < popt.size(); ++j)
		{
			data.read_in(counter);
			counter++;
		}
	}

	data.write("output.txt");


	

	return 0;
}