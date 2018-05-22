#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "auxiliary_files/vector.hpp"
#include "auxiliary_files/read_data.hpp"
#include "auxiliary_files/hmc.hpp"
#include "model.hpp"
#include <fstream>


int main ()
{
	Vector<double> x_data(0);
	Vector<double> y_data(0);
	Vector<double> dy_data(0);

	// Read in measured data
	read_data("Measurement_data_generator/measurements.txt", x_data, y_data, dy_data);

	// Vector for fitting parameters
	Vector<double> popt(4);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	range_min[0] = 0;
	range_max[0] = 6;
	range_min[1] = -1;
	range_max[1] = 5;
	range_min[2] = 0;
	range_max[2] = 2;
	range_min[3] = -0.5;
	range_max[3] = 1.5;

	//initialize HMC opbject
	HMC<double> sampler(x_data, y_data, dy_data, 6e-3, 25, 60, 2);
	
	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 200, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 400, 200, "preliminary_tools/correlation_times.txt");


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	//fill_from_region(popt, range_min, range_max);
	popt[0] = 2.31;
	popt[1] = 1.14;
	popt[2] = 1.10;
	popt[3] = 0.57;
	sampler.walk(3e4, 60*8, popt, 10);
	

	return 0;
}