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
	range_min[0] = 0.9852;
	range_max[0] = 0.9857;
	range_min[1] = 1.350;
	range_max[1] = 1.3675;
	range_min[2] = 0.2647;
	range_max[2] = 0.2660;
	range_min[3] = 0.518;
	range_max[3] = 0.54;
	range_min[4] = 0.06506;
	range_max[4] = 0.06539;
	range_min[5] = 0.040;
	range_max[5] = 0.0437;
	range_min[6] = 0.1755;
	range_max[6] = 0.1797;
	range_min[7] = 0.58;
	range_max[7] = 0.614;
	range_min[8] = 0.02245;
	range_max[8] = 0.022535;
	range_min[9] = 0.0363;
	range_max[9] = 0.0382;
	range_min[10] = 0.541;//
	range_max[10] = 0.5464;
	range_min[11] = -0.0053;//
	range_max[11] = -0.0044;
	range_min[12] = 0.8015;//
	range_max[12] = 0.807;
	range_min[13] = -0.0037;
	range_max[13] = -0.0030;//
	range_min[14] = -0.4;
	range_max[14] = 0;
	range_min[15] = -0.0073;//
	range_max[15] = -0.0062;
	range_min[16] = -1;
	range_max[16] = 0.0;
	range_min[17] = 1.288;//
	range_max[17] = 1.294;
	range_min[18] = -0.0055;//
	range_max[18] = -0.00465;
	range_min[19] = -0.8;
	range_max[19] = 0;
	range_min[20] = 2.18;
	range_max[20] = 2.21;//
	range_min[21] = 1.53;
	range_max[21] = 1.57;//
	range_min[22] = 2.126;
	range_max[22] = 2.147;//
	range_min[23] = 1.63;
	range_max[23] = 1.675;//
	range_min[24] = 1.3;
	range_max[24] = 3;
	range_min[25] = 1.56;
	range_max[25] = 1.60;//
	range_min[26] = 1.3;
	range_max[26] = 12.5;
	range_min[27] = 2.103;
	range_max[27] = 2.121;//
	range_min[28] = 1.692;
	range_max[28] = 1.731;//
	range_min[29] = 1.3;
	range_max[29] = 12;
	
	


	
	

	// Characteristic length scales for the parameters // default 1
	Vector<number_type> c_lengths(popt.size(), 1);
	c_lengths = range_max-range_min;

	// for the residual couplings constants and masses artificially decrease length scale
	for (size_type j = 10; j < c_lengths.size(); ++j)
	{
		if (j == 10 || j == 11 || j == 12 || j == 13 || j == 15 || j == 17 || j == 18)
			continue;
		if (j == 20 || j == 21 || j == 22 || j == 23 || j == 25 || j == 27 || j == 28)
			continue;
		if (j >= 20)
			c_lengths[j] /= 4.4e1;
		if (j >= 10 && j < 20)
			c_lengths[j] /= 0.9e1; 
		
	}
	

	


	


	//initialize HMC opbject
	HMC<number_type> sampler(correlators, range_min, range_max, c_lengths, 1e-3, 80, 120, 1e-3);
	//sampler.bounds_fixed(false);
	//sampler.do_analysis(true);
	

	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns the discretization error of H in a file (to adjust leapfrog step size)
	// sampler.get_H_errors(1e3, "preliminary_tools/H_errors.txt");
	// number_type mean;
	// number_type dev;
	// sampler.get_H(1e3, mean, dev);
	// std::cout << mean << " +/- " << dev <<"\n";


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
	
	std::string filenumber;
	if (argc > 1)
	{
		filenumber = argv[1];
	}

	//sampler.walk_automatic();
	//write_scripts(1, "script");
	
	sampler.walk_silently(1e4, "data", filenumber);

	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	// commented to avoid burn-in time (to be uncommented !!)
	
	// popt = 0;
	// popt[0] = 0.984527;
	// popt[1] = 1.35209;
	// popt[2] = 0.264117;
	// popt[3] = 0.523776;
	// popt[4] = 0.0644572;
	// popt[5] = 0.0440069;
	// popt[6] = 0.546241;
	// popt[7] = -0.00521682;
	// popt[8] = -0.000152412;
	// popt[9] = 2.18909;
	// popt[10] = 1.55419;
	// popt[11] = 1.38124;


	




	Vector<number_type> perr(popt.size());

	//sampler.intrinsic_err(popt, perr);

	
	//sampler.walk(5e4, 60*30, popt, 10);

	

	

	return 0;
}