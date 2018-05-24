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
	Vector<double> popt(5);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	range_min[0] = 2.24;
	range_max[0] = 2.35;
	range_min[1] = -1.4;
	range_max[1] = -0.9;
	range_min[2] = 0.2;
	range_max[2] = 0.4;
	range_min[3] = 1.88;
	range_max[3] = 2.08;
	range_min[4] = -2.58;
	range_max[4] = -2.08;

	//initialize HMC opbject
	HMC<double> sampler(x_data, y_data, dy_data, 5e-3, 10, 20, 0.01);
	
	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 400, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 500, 120, "preliminary_tools/correlation_times.txt");


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	fill_from_region(popt, range_min, range_max);
	 // popt[0] = 2.20;
	 // popt[1] = 1.03;
	 // popt[2] = 1.10;
	 // popt[3] = 0.58;
	sampler.walk(1e5, 60*8, popt, 10);
	

	return 0;
}