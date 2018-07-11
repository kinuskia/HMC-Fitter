#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "model_24.hpp"
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
	Vector<number_type> A0l_A0l(0);
	Vector<number_type> A0l_A0s(0);
	Vector<number_type> A0l_Pl(0);
	Vector<number_type> A0l_Ps(0);
	Vector<number_type> A0s_A0l(0);
	Vector<number_type> A0s_A0s(0);
	Vector<number_type> A0s_Pl(0);
	Vector<number_type> A0s_Ps(0);
	Vector<number_type> Pl_A0l(0);
	Vector<number_type> Pl_A0s(0);
	Vector<number_type> Pl_Pl(0);
	Vector<number_type> Pl_Ps(0);
	Vector<number_type> Ps_A0l(0);
	Vector<number_type> Ps_A0s(0);
	Vector<number_type> Ps_Pl(0);
	Vector<number_type> Ps_Ps(0);
	Vector<number_type> d_A0l_A0l(0);
	Vector<number_type> d_A0l_A0s(0);
	Vector<number_type> d_A0l_Pl(0);
	Vector<number_type> d_A0l_Ps(0);
	Vector<number_type> d_A0s_A0l(0);
	Vector<number_type> d_A0s_A0s(0);
	Vector<number_type> d_A0s_Pl(0);
	Vector<number_type> d_A0s_Ps(0);
	Vector<number_type> d_Pl_A0l(0);
	Vector<number_type> d_Pl_A0s(0);
	Vector<number_type> d_Pl_Pl(0);
	Vector<number_type> d_Pl_Ps(0);
	Vector<number_type> d_Ps_A0l(0);
	Vector<number_type> d_Ps_A0s(0);
	Vector<number_type> d_Ps_Pl(0);
	Vector<number_type> d_Ps_Ps(0);

	// Read in measured data, skip first row
	read_data("correlators/correlator_A0l-A0l", t, A0l_A0l, d_A0l_A0l, 1);
	read_data("correlators/correlator_A0l-A0s", t, A0l_A0s, d_A0l_A0s, 1);
	read_data("correlators/correlator_A0l-Pl", t, A0l_Pl, d_A0l_Pl, 1);
	read_data("correlators/correlator_A0l-Ps", t, A0l_Ps, d_A0l_Ps, 1);
	read_data("correlators/correlator_A0s-A0l", t, A0s_A0l, d_A0s_A0l, 1);
	read_data("correlators/correlator_A0s-A0s", t, A0s_A0s, d_A0s_A0s, 1);
	read_data("correlators/correlator_A0s-Pl", t, A0s_Pl, d_A0s_Pl, 1);
	read_data("correlators/correlator_A0s-Ps", t, A0s_Ps, d_A0s_Ps, 1);
	read_data("correlators/correlator_Pl-A0l", t, Pl_A0l, d_Pl_A0l, 1);
	read_data("correlators/correlator_Pl-A0s", t, Pl_A0s, d_Pl_A0s, 1);
	read_data("correlators/correlator_Pl-Pl", t, Pl_Pl, d_Pl_Pl, 1);
	read_data("correlators/correlator_Pl-Ps", t, Pl_Ps, d_Pl_Ps, 1);
	read_data("correlators/correlator_Ps-A0l", t, Ps_A0l, d_Ps_A0l, 1);
	read_data("correlators/correlator_Ps-A0s", t, Ps_A0s, d_Ps_A0s, 1);
	read_data("correlators/correlator_Ps-Pl", t, Ps_Pl, d_Ps_Pl, 1);
	read_data("correlators/correlator_Ps-Ps", t, Ps_Ps, d_Ps_Ps, 1);

	// Set up fitting model
	Model<number_type> correlators(
	t,
	A0l_A0l,
	A0l_A0s,
	A0l_Pl,
	A0l_Ps,
	A0s_A0l,
	A0s_A0s,
	A0s_Pl,
	A0s_Ps,
	Pl_A0l,
	Pl_A0s,
	Pl_Pl,
	Pl_Ps,
	Ps_A0l,
	Ps_A0s,
	Ps_Pl,
	Ps_Ps,
	d_A0l_A0l,
	d_A0l_A0s,
	d_A0l_Pl,
	d_A0l_Ps,
	d_A0s_A0l,
	d_A0s_A0s,
	d_A0s_Pl,
	d_A0s_Ps,
	d_Pl_A0l,
	d_Pl_A0s,
	d_Pl_Pl,
	d_Pl_Ps,
	d_Ps_A0l,
	d_Ps_A0s,
	d_Ps_Pl,
	d_Ps_Ps
	);

	//correlators.print_content();

	// Vector for fitting parameters
	Vector<number_type> popt(24);


	// Estimated search region
	Vector<number_type> range_min(popt.size());
	Vector<number_type> range_max(popt.size());
	// Characteristic length scales for the parameters // default 1
	Vector<number_type> c_lengths(popt.size(), 1);
	// characteristic length scales are here relative to range_max-range-min ...
	range_min[0] = 0.98528;
	range_max[0] = 0.98536;
	c_lengths[0] = 1.0;
	range_min[1] = 1.3508;
	range_max[1] = 1.3525;
	c_lengths[1] = 1.0;
	range_min[2] = 0.26475;
	range_max[2] = 0.26500;
	c_lengths[2] = 1.0;
	range_min[3] = 0.5193;
	range_max[3] = 0.5218;
	c_lengths[3] = 1.0;
	range_min[4] = 0.06520;
	range_max[4] = 0.06525;
	c_lengths[4] = 1.0;
	range_min[5] = 0.04064;
	range_max[5] = 0.04078;
	c_lengths[5] = 1.0;
	range_min[6] = 0.17866;
	range_max[6] = 0.17886;
	c_lengths[6] = 1.0;
	range_min[7] = 0.5816;
	range_max[7] = 0.5849;
	c_lengths[7] = 1.0;
	range_min[8] = 0.022470;
	range_max[8] = 0.022492;
	c_lengths[8] = 1.0;
	range_min[9] = 0.03639;
	range_max[9] = 0.03649;
	c_lengths[9] = 1.0;
	range_min[10] = 0.5444;//
	range_max[10] = 0.5458;
	c_lengths[10] = 1.0;
	range_min[11] = -0.00458;//
	range_max[11] = -0.00450;
	c_lengths[11] = 1.0;
	range_min[12] = 0.8053;//
	range_max[12] = 0.8063;
	c_lengths[12] = 1.0;
	range_min[13] = -0.00309;
	range_max[13] = -0.003034;//
	c_lengths[13] = 1.0;
	range_min[14] = -0.00637;//
	range_max[14] = -0.00629;
	c_lengths[14] = 1.0;
	range_min[15] = 1.2912;//
	range_max[15] = 1.2929;
	c_lengths[15] = 1.0;
	range_min[16] = -0.004787;//
	range_max[16] = -0.004716;
	c_lengths[16] = 1.0;
	range_min[17] = 2.181;
	range_max[17] = 2.185;//
	c_lengths[17] = 1.;
	range_min[18] = 1.543;
	range_max[18] = 1.553;//
	c_lengths[18] = 1.;
	range_min[19] = 2.1274;
	range_max[19] = 2.1296;//
	c_lengths[19] = 1.0;
	range_min[20] = 1.643;
	range_max[20] = 1.6505;//
	c_lengths[20] = 1.;
	range_min[21] = 1.579;
	range_max[21] = 1.585;//
	c_lengths[21] = 1.0;
	range_min[22] = 2.1043;
	range_max[22] = 2.1063;//
	c_lengths[22] = 1.0;
	range_min[23] = 1.711;
	range_max[23] = 1.717;//
	c_lengths[23] = 1.0;

	
	
	// ... and are now made absolute
	c_lengths *= (range_max-range_min);


	//initialize HMC opbject
	HMC<number_type> sampler(correlators, range_min, range_max, c_lengths, 1e-3, 90, 130, 1e-6);
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
	//sampler.walk(1e2, 10, 60*55, popt, 10);
	
	

	//sampler.walk_automatic();
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
	//sampler.walk_silently(1e4, "data", filenumber);
	//sampler.walk_silently_disregarding(1e4, 40.225, "data", filenumber);

	

	


	/* ACTUAL RUN */


	Vector<number_type> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	fill_from_region(popt, range_min, range_max);
	//sampler.walk(2e4, 60*30, popt, 10);
	sampler.walk(1e3, 60*30, popt, 10);

	

	

	return 0;
}