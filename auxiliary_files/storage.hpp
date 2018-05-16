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
	Storage (Vector<number_type> popt, size_type n_of_extra = 0)
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

	/* calculate Average vector */
	void average(Vector<number_type> & average_vector)
	{
		average_vector = 0;
		assert(average_vector.size() == n_variables_);
		for (size_type i = 0; i < data_.size(); ++i)
		{
			average_vector[i%n_variables_] += data_[i];
		}
		number_type entries_per_variable = data_.size()/n_variables_;
		average_vector = average_vector / entries_per_variable;
	}

	/* Renormalize by subtraction of the average vector  */
	void renormalize()
	{
		// Calculate average vector and subtract it from the data
		Vector<number_type> average_vector(n_variables_);
		average(average_vector);
		for (size_type i = 0; i < data_.size(); ++i)
		{
			data_[i] -= average_vector[i%n_variables_];
		}
	}

	/* Determine a vector containing the period of each signal */
	void period(Vector<size_type> & period_vector)
	{
		assert(period_vector.size() == n_variables_);
		renormalize();
		for (size_type i = 0; i < n_variables_; ++i)
		{
			// calculate autocorrelation function for variable i
			size_type entries_per_variable = data_.size() / n_variables_;
			Vector<number_type> gamma(entries_per_variable);

			// period corresponds to two sign changes of the autocorrelation function
			size_type counter_period = 0;
			size_type sign_changes = 0;
			bool sign = true; // true = +, false = -
			bool sign_old = sign;
			for (size_type j = 0; j < gamma.size(); ++j)
			{
				// calculate j-th term of the autocorrelation function
				number_type sum = 0;
				for (size_type k = 0; k < gamma.size(); ++k)
				{
					if (k > j)
					{
						break;
					}
					sum += data_[i + n_variables_*k]*data_[i + n_variables_*(j-k)];
				}
				gamma[j] = sum/gamma.size();

				sign = (gamma[j] >= 0);
				if (j > 0)
				{
					if (sign_changes == 2)
					{
						break; // stop calculation autocorrelation entries as soon as period found
					}
					else if (sign != sign_old)
					{
						sign_changes++;
					}
					counter_period++;
				}
				sign_old = sign;
			}
			period_vector[i] = counter_period;

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