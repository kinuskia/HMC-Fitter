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
	Vector<double> popt(10);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	// range_min[0] = -4.;
	// range_max[0] = 4.;
	// range_min[1] = 0.;
	// range_max[1] = 5.;
	// range_min[2] = 2.;
	// range_max[2] = 8.;
	// range_min[3] = 0.;
	// range_max[3] = 4.;
	// range_min[4] = 0.5;
	// range_max[4] = 4.;
	// range_min[5] = 2.;
	// range_max[5] = 6.;
	// range_min[6] = 0.01;
	// range_max[6] = 2.;
	// range_min[7] = 0.5;
	// range_max[7] = 5.;
	// range_min[8] = 6.;
	// range_max[8] = 8.;
	// range_min[9] = 0.01;
	// range_max[9] = 2.;
	range_min[0] = -1.1;
	range_max[0] = -1.0;
	range_min[1] = 1.0;
	range_max[1] = 1.2;
	range_min[2] = -1.6;
	range_max[2] = -1.4;
	range_min[3] = 2.9;
	range_max[3] = 3.1;
	range_min[4] = 3.9;
	range_max[4] = 4.1;
	range_min[5] = 3.8;
	range_max[5] = 4.0;
	range_min[6] = 1.6;
	range_max[6] = 1.8;
	range_min[7] = 4.2;
	range_max[7] = 4.4;
	range_min[8] = 7.0;
	range_max[8] = 7.2;
	range_min[9] = 0.5;
	range_max[9] = 0.7;

	// Characteristic length scales for the parameters // default 1
	Vector<double> c_lengths(popt.size(), 1);
	


	//initialize HMC opbject
	HMC<double> sampler(x_data, y_data, dy_data, range_min, range_max, c_lengths, 5e-4, 20, 30, 0.01);
	//sampler.bounds_fixed(true);
	//sampler.do_analysis(false);
	//sampler.discard_from(50);

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 200, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 300, 500, "preliminary_tools/correlation_times.txt");


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	fill_from_region(popt, range_min, range_max);
	
	sampler.walk(1e4, 60*30, popt, 10);
	

	return 0;
}