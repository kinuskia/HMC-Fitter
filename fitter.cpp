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
	range_min[0] = -4.;
	range_max[0] = 4.;
	range_min[1] = 0.;
	range_max[1] = 5.;
	range_min[2] = 2.;
	range_max[2] = 8.;
	range_min[3] = 0.;
	range_max[3] = 4.;
	range_min[4] = 0.5;
	range_max[4] = 4.;
	range_min[5] = 2.;
	range_max[5] = 6.;
	range_min[6] = 0.01;
	range_max[6] = 2.;
	range_min[7] = 0.5;
	range_max[7] = 5.;
	range_min[8] = 6.;
	range_max[8] = 8.;
	range_min[9] = 0.01;
	range_max[9] = 2.;
	// range_min[0] = -1.05;
	// range_max[0] = -1.0;
	// range_min[1] = 1.04;
	// range_max[1] = 1.17;
	// range_min[2] = -1.6;
	// range_max[2] = -1.43;
	// range_min[3] = 2.96;
	// range_max[3] = 3.04;
	// range_min[4] = 3.95;
	// range_max[4] = 4.0;
	// range_min[5] = 3.88;
	// range_max[5] = 3.91;
	// range_min[6] = 1.65;
	// range_max[6] = 1.70;
	// range_min[7] = 4.25;
	// range_max[7] = 4.31;
	// range_min[8] = 7.1275;
	// range_max[8] = 7.1375;
	// range_min[9] = 0.61;
	// range_max[9] = 0.62;

	// Characteristic length scales for the parameters // default 1
	Vector<double> c_lengths(popt.size(), 1);
	// c_lengths[0] = 0.05;
	// c_lengths[1] = 0.13;
	// c_lengths[2] = 0.17;
	// c_lengths[3] = 0.08;
	// c_lengths[4] = 0.05;
	// c_lengths[5] = 0.03;
	// c_lengths[6] = 0.05;
	// c_lengths[7] = 0.06;
	// c_lengths[8] = 0.01;
	// c_lengths[9] = 0.01;
	


	//initialize HMC opbject
	HMC<double> sampler(x_data, y_data, dy_data, range_min, range_max, c_lengths, 2e-4, 50, 70, 1e-1);
	sampler.bounds_fixed(true);
	sampler.do_analysis(false);
	sampler.discard_from(50);

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 50, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 300, 500, "preliminary_tools/correlation_times.txt");


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	
	fill_from_region(popt, range_min, range_max);
	// popt[0] = -1.02813235022;
	// popt[1] = 1.09070131382;
	// popt[2] = -1.49770662438;
	// popt[3] = 3.00762824513;
	// popt[4] = 3.97095733249;
	// popt[5] = 3.89567875037;
	// popt[6] = 1.67076043721;
	// popt[7] = 4.27963797423;
	// popt[8] = 7.13321500346;
	// popt[9] = 0.615462284398;
	Vector<double> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	sampler.walk(1e5, 60*30, popt, 10);
	

	return 0;
}