#ifndef STORAGE_HPP
#define STORAGE_HPP

#include "vector.hpp"
#include <string>
#include <fstream>

template<typename REAL>
class Storage
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Storage (Vector<number_type> popt, size_type n_of_extra)
	: n_popt_(popt.size())
	, n_extra_(n_of_extra) 
	, n_variables_(popt.size() + n_of_extra)
	, data_()
	{
	}

	/* receive and save input value */
	void read_in (number_type value)
	{
		data_.push_back(value);
	}
	/* receive and save input vector */
	void read_in (Vector<number_type> popt)
	{
		for (size_type i = 0; i < popt.size(); ++i)
		{
			data_.push_back(popt[i]);
		}
	}

	/* Write data in an output file */
	void write(std::string filename)
	{
		std::ofstream outfile(filename);
		size_type N = data_.size();
		for (size_type i = 0; i < N; ++i)
		{
			if ((i%n_variables_) == 0) // add counter at beginning of line
			{
				outfile << i/n_variables_ + 1;
				outfile << " ";
			}

			outfile << data_[i];
			
			if((i%n_variables_) == (n_variables_ - 1)) // new line
			{
				outfile << "\n";
				continue;
			}
			outfile << " ";
		}

	}


private:
	Vector<number_type> data_; // data vector
	size_type n_variables_; // number of variables stored in the vector
	size_type n_popt_; // number of fitting parameters
	size_type n_extra_; // number of additional variables
};


#endif 