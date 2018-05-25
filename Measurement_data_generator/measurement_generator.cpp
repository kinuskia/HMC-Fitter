#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "../auxiliary_files/vector.hpp"

double f(double x, Vector<double> popt)
{
	return popt[0] * atan(popt[1]*x +popt[2]) + popt[3]+ popt[4]*exp(-(x-popt[5])*(x-popt[5])/2/popt[6]/popt[6]) + popt[7]*exp(-(x-popt[8])*(x-popt[8])/2/popt[9]/popt[9]);
}


int main ()
{
	std::size_t n_data = 60;
	Vector<double> popt(10);
	popt[0] = -1.0;
	popt[1] = 1.2;
	popt[2] = -3.0;
	popt[3] = 3.0;
	popt[4] = 4.;
	popt[5] = 4.;
	popt[6] = 1.2;
	popt[7] = 5.;
	popt[8] = 7.;
	popt[9] = 0.7;
	
	Vector<double> pstd(popt.size());
	double rel = 0.02;
	double abs = 0.02;
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

	// save output to text file
	std::ofstream outfile("measurements.txt");
	for (std::size_t i = 0; i < n_data; ++i)
	{
		outfile << x[i] << "  " << y[i] << " " << dy[i] << "\n";
	}





	return 0;
}