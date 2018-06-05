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
	Vector<double> popt(20);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	range_min[0] = 0.986;
	range_max[0] = 0.9935;
	range_min[1] = 1.40;
	range_max[1] = 1.535;
	range_min[2] = -0.29;
	range_max[2] = 0.29;
	range_min[3] = 0.50;
	range_max[3] = 0.85;
	range_min[4] = -0.071;
	range_max[4] = 0.072;
	range_min[5] = 0.09;
	range_max[5] = 0.0185;
	range_min[6] = -0.20;
	range_max[6] = 0.20;
	range_min[7] = 0.67;
	range_max[7] = 1.1;
	range_min[8] = 0.46;
	range_max[8] = 0.9;
	range_min[9] = -0.85;
	range_max[9] = 0.8;
	range_min[10] = 0.72;
	range_max[10] = 0.88;
	range_min[11] = -0.74;
	range_max[11] = 0.94;
	range_min[12] = -0.95;
	range_max[12] = 1.34;
	range_min[13] = 1.13;
	range_max[13] = 1.38;
	range_min[14] = 2.0;
	range_max[14] = 3.15;
	range_min[15] = 6.3;
	range_max[15] = 17.8;
	range_min[16] = 2.05;
	range_max[16] = 2.6;
	range_min[17] = 6;
	range_max[17] = 26.5;
	range_min[18] = 5.8;
	range_max[18] = 17;
	range_min[19] = 2.07;
	range_max[19] = 2.43;


	
	

	// Characteristic length scales for the parameters // default 1
	Vector<double> c_lengths(popt.size(), 1);
	c_lengths = range_max-range_min;


	


	//initialize HMC opbject
	HMC<double> sampler(correlators, range_min, range_max, c_lengths, 5e-6, 40, 50, 1e-2);
	sampler.bounds_fixed(false);
	sampler.do_analysis(true);
	

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 50, 50, "preliminary_tools/acceptrates.txt");
	
	// returns autocorrelation lengths for random starting points (to adjust number of leapfrog steps)
	//sampler.get_optimal_number_of_steps(range_min, range_max, 300, 700, "preliminary_tools/correlation_times.txt");

	/*
	For the estimation of the lower and upper bounds of the parameters, values with a potential above
	the value below are discarded
	*/
	sampler.discard_from(150.);
	fill_from_region(popt, range_min, range_max);
	sampler.walk(3e4, 4, 60*55, popt, 10);

	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	
	



	Vector<double> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	//sampler.walk(5e4, 60*30, popt, 10);

	

	

	return 0;
}