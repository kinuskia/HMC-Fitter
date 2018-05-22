#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include <random>
#include <assert.h>

typedef double number_type;
typedef std::size_t size_type;

number_type gamma (const Vector<number_type> & q, size_type t)
{
	size_type N = q.size();
	number_type q_mean = q.mean();
	assert(N > t);
	number_type result = 0;
	for (size_type i = 0; i < (N-t); ++i)
	{
		result += (q[i] - q_mean)*(q[i+t] - q_mean);
	}
	result /= (N-t);
	return result;
}

number_type potential (number_type x)
{
	return 0.25*(x-1)*(x-1);
}

number_type correlation_time(const Vector<number_type> & q, size_type W)
{
	number_type C = gamma(q, 0);
	for (size_type t = 1; t <= W; ++t)
	{
		C += 2*gamma(q, t);
	}
	number_type v = gamma(q, 0);

	number_type result = C/2/v;
	return result;
}

void correlation_time(const Vector<number_type> & q, Vector<number_type> & result)
{
	assert(result.size() == 2);
	number_type C = gamma(q, 0);
	size_type W = 0;
	number_type gamma_old;
	for (size_type t = 1; t <= q.size(); ++t)
	{
		C += 2*gamma(q, t);
		if (t>1)
		{
			if (gamma_old < gamma(q, t) || gamma(q, t) < 0)
			{
				W = t;
				std::cout << W << "\n";
				break;
			}
			
		}
		gamma_old = gamma(q, t);
	}
	number_type v = gamma(q, 0);
	result[0] = C/2/v;
	result[1] = result[0]*sqrt(4./q.size()*(W+1./2-result[0])) + C*exp(-1.*W/result[0]);
	std::cout << "(stat:) " << result[0]*sqrt(4./q.size()*(W+1./2-result[0])) << " + (syst.) " << C*exp(-1.*W/result[0]) << "\n";
}




int main ()
{

	
	Vector<number_type> q(10000);
	std::random_device rd;
	std::mt19937 gen(rd());

	q[0] = 1;
	// number_type p = -2;
	// for (size_type i = 1; i < q.size(); ++i)
	// {
	// 	p += -q[i-1]/2;
	// 	q[i] = q[i-1] + p;
	// }
	for (size_type i = 1; i < q.size(); ++i)
	{
		std::normal_distribution<> dis(q[i-1], 1);
		number_type candidate = dis(gen);
		number_type energy_diff = potential(candidate) - potential(q[i-1]);
		if (energy_diff < 0)
		{
			q[i] = candidate;
		}
		else
		{
			std::uniform_real_distribution<> dis_unif(0, 1);
			if (dis_unif(gen) < exp(-energy_diff))
			{
				q[i] = candidate;
			}
			else
			{
				q[i] = q[i-1];
			}
		}
	}
	q.save_as("sequence.txt");
	Vector<number_type> gammas(0);

	for (size_type t = 0; t < 0.05*q.size(); ++t)
	{
		gammas.push_back(gamma(q, t));
	}
	gammas.save_as("gammas.txt");

	Vector<number_type> times(gammas.size());
	for (size_type i = 0; i<times.size(); ++i)
	{
		times[i] = correlation_time(q, i+1);
	}
	times.save_as("times.txt");

	Vector<number_type> result(2);
	correlation_time(q, result);
	std::cout << "tau_int: " << result[0] << " + - " << result[1] << "\n";



	return 0;
}