#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "model.hpp"
#include "auxiliary_files/vector.hpp"
#include "auxiliary_files/read_data.hpp"
#include "auxiliary_files/write_scripts.hpp"
#include "auxiliary_files/hmc.hpp"
#include <fstream>


int main (int argc, char* argv[])
{
	typedef double number_type;
	typedef std::size_t size_type;
	Vector<number_type> t(0);
	Vector<number_type> Vl_Vl(0);
	Vector<number_type> Vl_Vs1(0);
	Vector<number_type> Vl_Vs2(0);
	Vector<number_type> Vl_Vs3(0);
	Vector<number_type> Vs1_Vl(0);
	Vector<number_type> Vs1_Vs1(0);
	Vector<number_type> Vs1_Vs2(0);
	Vector<number_type> Vs1_Vs3(0);
	Vector<number_type> Vs2_Vl(0);
	Vector<number_type> Vs2_Vs1(0);
	Vector<number_type> Vs2_Vs2(0);
	Vector<number_type> Vs2_Vs3(0);
	Vector<number_type> Vs3_Vl(0);
	Vector<number_type> Vs3_Vs1(0);
	Vector<number_type> Vs3_Vs2(0);
	Vector<number_type> Vs3_Vs3(0);
	Vector<number_type> d_Vl_Vl(0);
	Vector<number_type> d_Vl_Vs1(0);
	Vector<number_type> d_Vl_Vs2(0);
	Vector<number_type> d_Vl_Vs3(0);
	Vector<number_type> d_Vs1_Vl(0);
	Vector<number_type> d_Vs1_Vs1(0);
	Vector<number_type> d_Vs1_Vs2(0);
	Vector<number_type> d_Vs1_Vs3(0);
	Vector<number_type> d_Vs2_Vl(0);
	Vector<number_type> d_Vs2_Vs1(0);
	Vector<number_type> d_Vs2_Vs2(0);
	Vector<number_type> d_Vs2_Vs3(0);
	Vector<number_type> d_Vs3_Vl(0);
	Vector<number_type> d_Vs3_Vs1(0);
	Vector<number_type> d_Vs3_Vs2(0);
	Vector<number_type> d_Vs3_Vs3(0);

	// Read in measured data, skip first row
	read_data("correlators/correlator_Vl-Vl", t, Vl_Vl, d_Vl_Vl, 0);
	read_data("correlators/correlator_Vl-Vs1", t, Vl_Vs1, d_Vl_Vs1, 0);
	read_data("correlators/correlator_Vl-Vs2", t, Vl_Vs2, d_Vl_Vs2, 0);
	read_data("correlators/correlator_Vl-Vs3", t, Vl_Vs3, d_Vl_Vs3, 0);
	read_data("correlators/correlator_Vs1-Vl", t, Vs1_Vl, d_Vs1_Vl, 0);
	read_data("correlators/correlator_Vs1-Vs1", t, Vs1_Vs1, d_Vs1_Vs1, 0);
	read_data("correlators/correlator_Vs1-Vs2", t, Vs1_Vs2, d_Vs1_Vs2, 0);
	read_data("correlators/correlator_Vs1-Vs3", t, Vs1_Vs3, d_Vs1_Vs3, 0);
	read_data("correlators/correlator_Vs2-Vl", t, Vs2_Vl, d_Vs2_Vl, 0);
	read_data("correlators/correlator_Vs2-Vs1", t, Vs2_Vs1, d_Vs2_Vs1, 0);
	read_data("correlators/correlator_Vs2-Vs2", t, Vs2_Vs2, d_Vs2_Vs2, 0);
	read_data("correlators/correlator_Vs2-Vs3", t, Vs2_Vs3, d_Vs2_Vs3, 0);
	read_data("correlators/correlator_Vs3-Vl", t, Vs3_Vl, d_Vs3_Vl, 0);
	read_data("correlators/correlator_Vs3-Vs1", t, Vs3_Vs1, d_Vs3_Vs1, 0);
	read_data("correlators/correlator_Vs3-Vs2", t, Vs3_Vs2, d_Vs3_Vs2, 0);
	read_data("correlators/correlator_Vs3-Vs3", t, Vs3_Vs3, d_Vs3_Vs3, 0);

	// Set up fitting model
	Model<number_type> correlators(
	t,
	Vl_Vl,
	Vl_Vs1,
	Vl_Vs2,
	Vl_Vs3,
	Vs1_Vl,
	Vs1_Vs1,
	Vs1_Vs2,
	Vs1_Vs3,
	Vs2_Vl,
	Vs2_Vs1,
	Vs2_Vs2,
	Vs2_Vs3,
	Vs3_Vl,
	Vs3_Vs1,
	Vs3_Vs2,
	Vs3_Vs3,
	d_Vl_Vl,
	d_Vl_Vs1,
	d_Vl_Vs2,
	d_Vl_Vs3,
	d_Vs1_Vl,
	d_Vs1_Vs1,
	d_Vs1_Vs2,
	d_Vs1_Vs3,
	d_Vs2_Vl,
	d_Vs2_Vs1,
	d_Vs2_Vs2,
	d_Vs2_Vs3,
	d_Vs3_Vl,
	d_Vs3_Vs1,
	d_Vs3_Vs2,
	d_Vs3_Vs3
	);

	//correlators.print_content();

	// Vector for fitting parameters
	Vector<number_type> popt(12);


	// Estimated search region
	Vector<number_type> range_min(popt.size());
	Vector<number_type> range_max(popt.size());
	// Characteristic length scales for the parameters // default 1
	Vector<number_type> c_lengths(popt.size(), 1);
	// characteristic length scales are here relative to range_max-range-min ...
	range_min[0] = 0.4;
	range_max[0] = 2.5;
	c_lengths[0] = 1;
	range_min[1] = 1.0;
	range_max[1] = 5.0;
	c_lengths[1] = 1;
	range_min[2] = 0;
	range_max[2] = 1;
	c_lengths[2] = 1.0;
	range_min[3] = 0;
	range_max[3] = 1;
	c_lengths[3] = 1.0;
	range_min[4] = -1;
	range_max[4] = 1;
	c_lengths[4] = 1.0;
	range_min[5] = -1;
	range_max[5] = 1;
	c_lengths[5] = 1.0;
	range_min[6] = -0.1;
	range_max[6] = 0.1;
	c_lengths[6] = 1.0;
	range_min[7] = -0.1;
	range_max[7] = 0.1;
	c_lengths[7] = 1.0;
	range_min[8] = -0.1;
	range_max[8] = 0.1;
	c_lengths[8] = 1.0;
	range_min[9] = 1.6;
	range_max[9] = 3;
	c_lengths[9] = 1.0;
	range_min[10] = 1.6;//
	range_max[10] = 3;
	c_lengths[10] = 1.0;
	range_min[11] = 1.6;//
	range_max[11] = 3;
	c_lengths[11] = 1.0;

	

	
	
	// ... and are now made absolute
	c_lengths *= (range_max-range_min);


	//initialize HMC opbject
	HMC<number_type> sampler(correlators, range_min, range_max, c_lengths, 1e-3, 90, 130, 1e6);
	//sampler.bounds_fixed(false);
	//sampler.do_analysis(true);
	

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns the discretization error of H in a file (to adjust leapfrog step size)
	// sampler.get_H_errors(1e3, "preliminary_tools/H_errors.txt");


	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 50, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 300, 700, "preliminary_tools/correlation_times.txt");

	/*
	For the estimation of the lower and upper bounds of the parameters, values with a potential above
	the value below are discarded
	*/
	//sampler.discard_from(150.);
	//fill_from_region(popt, range_min, range_max);
	//sampler.walk(2e3, 1, 60*55, popt, 10);
	
	

	sampler.walk_automatic();
	//write_scripts(1, "script");


	/* 	
		Routine for the computer cluster: Produce Markov chain of specific length and save to output file
		of the following format dataX.txt, where X is given in the console when executing, e.g. ./fitter 5
	*/
	std::string filenumber;
	if (argc > 1)
	{
		filenumber = argv[1];
	}
	//sampler.walk_silently(2e3, "data", filenumber);
	//sampler.walk_silently_disregarding(1e4, 40.225, "data", filenumber);

	

	


	/* ACTUAL RUN */


	Vector<number_type> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	//fill_from_region(popt, range_min, range_max);
	//sampler.walk(2e4, 60*30, popt, 10);

	

	

	return 0;
}