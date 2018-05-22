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

	// TODO: Include bool for whether data has been read in


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

	/* calculate average of variable i */
	number_type mean(size_type index)
	{
		assert((index >= 0) && (index < n_variables_)); //check if index within valid bounds
		number_type result = 0;
		for (size_type i = index; i < data_.size(); i+=n_variables_)
		{
			result += data_[i];
		}
		number_type entries_per_variable = data_.size()/n_variables_;
		result /= entries_per_variable;
		return result;

	}

	/* calculate mean of each variable at once and save it into a vector */
	/* as well as uncertitude vector, taking into account autocorrelation */
	void mean(Vector<number_type> & average_vector, Vector<number_type> & err_vector)
	{
		assert(average_vector.size() == err_vector.size());
		assert(average_vector.size() == n_popt_);
		// Calculation of the mean vector
		average_vector = 0;
		for (size_type i = 0; i < data_.size(); ++i)
		{
			average_vector[i%n_variables_] += data_[i];
		}
		number_type entries_per_variable = data_.size()/n_variables_;
		average_vector = average_vector / entries_per_variable;

		// Calculation of its uncertainty
		Vector<number_type> autocorr_times(n_popt_);
		Vector<number_type> autocorr_times_err(n_popt_);
		autocorr_time(autocorr_times, autocorr_times_err); // calculate integrated autocorrelation times
		for (size_type i = 0; i<n_popt_; ++i)
		{
			number_type ESS = entries_per_variable/(2*autocorr_times[i]); //effective sample size
			number_type variance = this->variance(i);
			err_vector[i] = sqrt(variance/ESS);
		}
	}

	/* 
	Autocorrelation function \Gamma_{\alpha, \beta}(n)
	for correlation between parameters with indices alpha and beta with time lag n 
	*/
	number_type gamma(size_type alpha, size_type beta, size_type lag)
	{
		size_type entries_per_variable = data_.size()/n_variables_;
		number_type mean_alpha = this->mean(alpha);
		number_type mean_beta = this->mean(beta);
		assert(entries_per_variable > lag); // condition for validity of this estimator
		number_type result = 0;
		for (size_type i = 0; i < (entries_per_variable-lag); ++i)
		{
			result += (data_[alpha + i*n_variables_] - mean_alpha) * (data_[beta + (i+lag)*n_variables_] - mean_beta);
		}
		result /= (entries_per_variable-lag);
		return result;
	}


	/* calculate a integrated autocorrelation time vector (one entry for each fitting parameter) */
	void autocorr_time(Vector<number_type> & times, Vector<number_type> & times_err)
	{
		assert(times.size() == n_popt_); // check for correct dimensions
		assert(times_err.size() == n_popt_);
		size_type entries_per_variable = data_.size()/n_variables_;
		for (size_type i = 0; i < times.size(); ++i)
		{
			// Calculate integrated autocorrelation time for parameter i
			number_type C = abs(this->gamma(i, i, 0));
			size_type W = 0;
			number_type gamma_old;
			number_type gamma_current;
			for (size_type t = 1; t<=entries_per_variable; ++t)
			{
				gamma_current = abs(this->gamma(i, i, t));
				C+= 2*gamma_current;

				/* 
				procedure to find optimal summation window W: 
				break as soon as module of autocorrelation function increases 
				*/
				if (t>1)
				{
					if (gamma_old < gamma_current)
					{
						W = t;
						break;
					}
				}

				gamma_old = abs(gamma_current);
			}
			number_type v = abs(this->gamma(i, i, 0));
			times[i] = C/2.0/v;
			times_err[i] = times[i]*sqrt(4./entries_per_variable*(W+1./2-times[i])) + C*exp(-1.*W/times[i]);
		}
	}

	/* calculate (unbiased) variance of variable i */
	number_type variance(size_type index)
	{
		number_type mean = this->mean(index);
		number_type result = 0;
		for (size_type i = index; i < data_.size(); i+=n_variables_)
		{
			result += (data_[i] - mean)*(data_[i] - mean);
		}
		number_type entries_per_variable = data_.size()/n_variables_;
		result /= (entries_per_variable -1);
		return result;
	}

	/* Calculate unbiased standard deviation of variable i */
	number_type std_deviation(size_type index)
	{
		number_type result = sqrt(variance(index));
		return result;
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
		assert(period_vector.size() == n_popt_);

		for (size_type i = 0; i < n_popt_; ++i)
		{
			size_type entries_per_variable = data_.size()/n_variables_;
			/* 
			To calculate the period of the signal, we determine the
			number of steps for two sign changes
			*/
			size_type counter_period = 0;
			size_type sign_changes = 0;
			bool sign = true; // true = '+', false = '-'
			bool sign_old = sign;
			for (size_type j = 0; j < entries_per_variable; ++j)
			{
				number_type gamma = this->gamma(i, i, j);
				sign = (gamma >= 0);
				if (j > 0)
				{
					if (sign_changes == 2)
					{
						break; // stop calculating once period found
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