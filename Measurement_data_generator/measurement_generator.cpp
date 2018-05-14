#include <iostream>
#include <cmath>
#include <random>
#include "../vector.hpp"

double f(double x, Vector<double> popt)
{
	return popt[0] + popt[1]*sin(popt[2]*x);
}


int main ()
{
	std::size_t n_data = 30;
	Vector<double> popt(3);
	popt[0] = 2.3;
	popt[1] = 1.2;
	popt[2] = 1.1;
	//popt[3] = 0.5;
	Vector<double> pstd(popt.size());
	double rel = 0.03;
	double abs = 0.03;
	pstd = rel*popt + abs;

	double x_min = 0;
	double x_max = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis_unif(x_min, x_max);

	Vector<double> x(n_data);
	Vector<double> y(n_data);
	Vector<double> dy(n_data);

	// Create x values
	for (std::size_t i = 0; i < n_data; ++i)
	{
		x[i] = dis_unif(gen);
	}
	std::sort(x.begin(), x.end());
	

	// Create y and dy values
	for (std::size_t i = 0; i<x.size(); ++i)
	{
		Vector<double> fit_parameters(popt.size());
		for (int j = 0; j < popt.size(); ++j)
		{
			std::normal_distribution<> dis_norm(popt[j], pstd[j]);
			//std::uniform_real_distribution<> dis_unif(popt[j]-sqrt(3)*pstd[j], popt[j]+sqrt(3)*pstd[j]);
			//std::gamma_distribution<> dis_norm(popt[j]*popt[j]/pstd[j]/pstd[j], pstd[j]*pstd[j]/popt[j]);
			fit_parameters[j] = dis_norm(gen);
		}
		y[i] = f(x[i], fit_parameters);
		dy[i] = rel*y[i] + abs;
	}

	// Print output
	for (std::size_t i = 0; i < n_data; ++i)
	{
		std::cout << x[i] << "  " << y[i] << " " << dy[i] << "\n";
	}





	return 0;
}