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
	range_min[0] = 0.981644;
	range_max[0] = 0.994899;
	range_min[1] = 1.24214;
	range_max[1] = 1.67238;
	range_min[2] = -0.290866;
	range_max[2] = 0.291909;
	range_min[3] = 0.368444;
	range_max[3] = 1.08222;
	range_min[4] = -0.0704209;
	range_max[4] = 0.0704742;
	range_min[5] = 0.0077663;
	range_max[5] = 0.0280702;
	range_min[6] = -0.201429;
	range_max[6] = 0.202556;
	range_min[7] = 0.391631;
	range_max[7] = 1.45834;
	range_min[8] = 0.40843;
	range_max[8] = 2.18364;
	range_min[9] = -11.6915;
	range_max[9] = 8.77749;
	range_min[10] = 0.644657;
	range_max[10] = 0.876352;
	range_min[11] = -8.29216;
	range_max[11] = 7.30805;
	range_min[12] = -12.5464;
	range_max[12] = 9.72157;
	range_min[13] = 0.996292;
	range_max[13] = 1.40714;
	range_min[14] = 1.99583;
	range_max[14] = 4.57368;
	range_min[15] = 5.97847;
	range_max[15] = 110.062;
	range_min[16] = 2.00568;
	range_max[16] = 3.17356;
	range_min[17] = 6.79197;
	range_max[17] = 154.174;
	range_min[18] = 5.63385;
	range_max[18] = 191.832;
	range_min[19] = 1.98952;
	range_max[19] = 2.82087;


	
	

	// Characteristic length scales for the parameters // default 1
	Vector<double> c_lengths(popt.size(), 1);
	c_lengths = range_max-range_min;


	


	//initialize HMC opbject
	HMC<double> sampler(correlators, range_min, range_max, c_lengths, 5e-6, 40, 50, 1e1);
	//sampler.bounds_fixed(false);
	//sampler.do_analysis(true);
	

	/* PRELIMINARY RUN TOOLS */
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
	//sampler.walk(3e4, 4, 60*55, popt, 10);
	sampler.walk_automatic();

	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	
	



	Vector<double> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	//sampler.walk(5e4, 60*30, popt, 10);

	

	

	return 0;
}