#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "auxiliary_files/vector.hpp"
#include "auxiliary_files/read_data.hpp"
#include "auxiliary_files/hmc.hpp"
#include "model.hpp"
#include "model.hpp"
#include <fstream>


int main ()
{
	Vector<double> t(0);
	Vector<double> A0l_A0l(0);
	Vector<double> A0l_A0s(0);
	Vector<double> A0l_Pl(0);
	Vector<double> A0l_Ps(0);
	Vector<double> A0s_A0l(0);
	Vector<double> A0s_A0s(0);
	Vector<double> A0s_Pl(0);
	Vector<double> A0s_Ps(0);
	Vector<double> Pl_A0l(0);
	Vector<double> Pl_A0s(0);
	Vector<double> Pl_Pl(0);
	Vector<double> Pl_Ps(0);
	Vector<double> Ps_A0l(0);
	Vector<double> Ps_A0s(0);
	Vector<double> Ps_Pl(0);
	Vector<double> Ps_Ps(0);
	Vector<double> d_A0l_A0l(0);
	Vector<double> d_A0l_A0s(0);
	Vector<double> d_A0l_Pl(0);
	Vector<double> d_A0l_Ps(0);
	Vector<double> d_A0s_A0l(0);
	Vector<double> d_A0s_A0s(0);
	Vector<double> d_A0s_Pl(0);
	Vector<double> d_A0s_Ps(0);
	Vector<double> d_Pl_A0l(0);
	Vector<double> d_Pl_A0s(0);
	Vector<double> d_Pl_Pl(0);
	Vector<double> d_Pl_Ps(0);
	Vector<double> d_Ps_A0l(0);
	Vector<double> d_Ps_A0s(0);
	Vector<double> d_Ps_Pl(0);
	Vector<double> d_Ps_Ps(0);

	// Read in measured data
	read_data("correlators/correlator_A0l-A0l", t, A0l_A0l, d_A0l_A0l);
	read_data("correlators/correlator_A0l-A0s", t, A0l_A0s, d_A0l_A0s);
	read_data("correlators/correlator_A0l-Pl", t, A0l_Pl, d_A0l_Pl);
	read_data("correlators/correlator_A0l-Ps", t, A0l_Ps, d_A0l_Ps);
	read_data("correlators/correlator_A0s-A0l", t, A0s_A0l, d_A0s_A0l);
	read_data("correlators/correlator_A0s-A0s", t, A0s_A0s, d_A0s_A0s);
	read_data("correlators/correlator_A0s-Pl", t, A0s_Pl, d_A0s_Pl);
	read_data("correlators/correlator_A0s-Ps", t, A0s_Ps, d_A0s_Ps);
	read_data("correlators/correlator_Pl-A0l", t, Pl_A0l, d_Pl_A0l);
	read_data("correlators/correlator_Pl-A0s", t, Pl_A0s, d_Pl_A0s);
	read_data("correlators/correlator_Pl-Pl", t, Pl_Pl, d_Pl_Pl);
	read_data("correlators/correlator_Pl-Ps", t, Pl_Ps, d_Pl_Ps);
	read_data("correlators/correlator_Ps-A0l", t, Ps_A0l, d_Ps_A0l);
	read_data("correlators/correlator_Ps-A0s", t, Ps_A0s, d_Ps_A0s);
	read_data("correlators/correlator_Ps-Pl", t, Ps_Pl, d_Ps_Pl);
	read_data("correlators/correlator_Ps-Ps", t, Ps_Ps, d_Ps_Ps);

	// Set up fitting model
	Model<double> correlators(
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
	Vector<double> popt(6);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	range_min[0] = 0.8;
	range_max[0] = 1.1;
	range_min[1] = 1.2;
	range_max[1] = 1.7;
	range_min[2] = 0.1;
	range_max[2] = 0.4;
	range_min[3] = 0.4;
	range_max[3] = 1.1;
	range_min[4] = 0.4;
	range_max[4] = 0.9;
	range_min[5] = 2.1;
	range_max[5] = 4;
	
	

	// Characteristic length scales for the parameters // default 1
	Vector<double> c_lengths(popt.size(), 1);
	// c_lengths[0] = 0.0015;
	// c_lengths[1] = 0.02;
	// c_lengths[2] = 0.003;
	// c_lengths[3] = 0.035;
	// c_lengths[4] = 0.015;
	// c_lengths[5] = 0.075;


	


	//initialize HMC opbject
	HMC<double> sampler(correlators, range_min, range_max, c_lengths, 9e-4, 4, 5, 1e0);
	sampler.bounds_fixed(false);
	sampler.do_analysis(true);
	//sampler.discard_from(50);

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	sampler.get_acceptance_rates(range_min, range_max, 200, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 300, 700, "preliminary_tools/correlation_times.txt");


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	
	fill_from_region(popt, range_min, range_max);

	Vector<double> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	//sampler.walk(1e4, 10, 60*30, popt, 10);
	//sampler.walk(1e5, 60*30, popt, 10);
	

	return 0;
}