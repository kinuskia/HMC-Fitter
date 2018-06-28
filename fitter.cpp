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
	Vector<number_type> popt(30);

	// Boolean vector for whether a fitting parameter is likely to be degenerate in terms of (+) <-> (-)
	Vector<bool> degenerate(popt.size());
	degenerate = false;
	// degenerate[2] = true;
	// degenerate[3] = true;
	// degenerate[4] = true;
	// degenerate[5] = true;

	// Estimated search region
	Vector<number_type> range_min(popt.size());
	Vector<number_type> range_max(popt.size());
	// Characteristic length scales for the parameters // default 1
	Vector<number_type> c_lengths(popt.size(), 1);
	// characteristic length scales are here relative to range_max-range-min ...
	range_min[0] = 0.9850;
	range_max[0] = 0.9857;
	c_lengths[0] = 1.;
	range_min[1] = 1.348;
	range_max[1] = 1.371;
	c_lengths[1] = 0.5;
	range_min[2] = 0.2644;
	range_max[2] = 0.2665;
	c_lengths[2] = 1.;
	range_min[3] = 0.516;
	range_max[3] = 0.549;
	c_lengths[3] = 0.5;
	range_min[4] = 0.06497;
	range_max[4] = 0.0655;
	c_lengths[4] = 1.0;
	range_min[5] = 0.0396;
	range_max[5] = 0.0439;
	c_lengths[5] = 0.5;
	range_min[6] = 0.1784;
	range_max[6] = 0.1801;
	c_lengths[6] = 1.;
	range_min[7] = 0.578;
	range_max[7] = 0.622;
	c_lengths[7] = 0.5;
	range_min[8] = 0.02241;
	range_max[8] = 0.02258;
	c_lengths[8] = 1.;
	range_min[9] = 0.0362;
	range_max[9] = 0.0385;
	c_lengths[9] = 0.3;
	range_min[10] = 0.538;//
	range_max[10] = 0.549;
	c_lengths[10] = 1.;
	range_min[11] = -0.0055;//
	range_max[11] = -0.00425;
	c_lengths[11] = 0.5;
	range_min[12] = 0.799;//
	range_max[12] = 0.808;
	c_lengths[12] = 1.;
	range_min[13] = -0.0039;
	range_max[13] = -0.0029;//
	c_lengths[13] = 0.5;
	range_min[14] = -0.125;
	range_max[14] = 0;
	c_lengths[14] = 0.5;
	range_min[15] = -0.0074;//
	range_max[15] = -0.0061;
	c_lengths[15] = 0.5;
	range_min[16] = -2.;
	range_max[16] = 0;
	c_lengths[16] = 1./21.;
	range_min[17] = 1.284;//
	range_max[17] = 1.297;
	c_lengths[17] = 1.0;
	range_min[18] = -0.0057;//
	range_max[18] = -0.0046;
	c_lengths[18] = 0.5;
	range_min[19] = -0.9;
	range_max[19] = 0.0;
	c_lengths[19] = 1./7.;
	range_min[20] = 2.175;
	range_max[20] = 2.219;//
	c_lengths[20] = 1.;
	range_min[21] = 1.515;
	range_max[21] = 1.592;//
	c_lengths[21] = 1.;
	range_min[22] = 2.124;
	range_max[22] = 2.152;//
	c_lengths[22] = 0.5;
	range_min[23] = 1.619;
	range_max[23] = 1.685;//
	c_lengths[23] = 1.;
	range_min[24] = 1.3;
	range_max[24] = 6;
	c_lengths[24] = 1./9.;
	range_min[25] = 1.552;
	range_max[25] = 1.611;//
	c_lengths[25] = 1.;
	range_min[26] = 1.3;
	range_max[26] = 12.5;
	c_lengths[26] = 1./12;
	range_min[27] = 2.10;
	range_max[27] = 2.1250;//
	c_lengths[27] = 0.5;
	range_min[28] = 1.68;
	range_max[28] = 1.74;//
	c_lengths[28] = 0.5;
	range_min[29] = 1.3;
	range_max[29] = 12;
	c_lengths[29] = 1./4;
	
	
	// ... and are now made absolute
	c_lengths *= (range_max-range_min);

	for (size_type j = 10; j < c_lengths.size(); ++j)
	{
		if (j == 14 || j == 16 || j == 19 || j == 24 || j == 26 || j == 29)
		{
			c_lengths[j] *= 6e-2;
		}
	}

	// artificially decrease length scale for the residual couplings constants and masses 
	// for (size_type j = 10; j < c_lengths.size(); ++j)
	// {
	// 	if (j == 10 || j == 11 || j == 12 || j == 13 || j == 15 || j == 17 || j == 18)
	// 		continue;
	// 	if (j == 20 || j == 21 || j == 22 || j == 23 || j == 25 || j == 27 || j == 28)
	// 		continue;
	// 	if (j >= 20)
	// 		c_lengths[j] /= 1.0e0; //4.4e1
	// 	if (j >= 10 && j < 20)
	// 		c_lengths[j] /= 1.0e0; //1.3e1
		
	// }
	

	


	


	//initialize HMC opbject
	HMC<number_type> sampler(correlators, range_min, range_max, c_lengths, 3e-3, 90, 130, 1e-3);
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
	sampler.walk_silently(1e4, "data", filenumber);

	/* ACTUAL RUN */


	Vector<number_type> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	//sampler.walk(5e4, 60*30, popt, 10);

	

	

	return 0;
}