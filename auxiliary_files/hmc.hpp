#ifndef HMC_HPP
#define HMC_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <assert.h>
#include <random>
#include <ctime>
#include "storage.hpp"
#include "../model.hpp"
#include <algorithm>
#include <fstream>

template<typename REAL>
class HMC
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	/* Standard constructor, fed with measured data */
	HMC(Model<number_type> model, Vector<number_type> range_min, Vector<number_type> range_max, Vector<number_type> c_lengths, number_type stepsize, size_type n_steps_min, size_type n_steps_max, number_type temperature)
	: range_min_(range_min)
	, range_max_(range_max) 
	, bounds_fixed_(false)
	, do_analysis_(true)
	, discard_data_(false)
	, chi2redmax_(0)
	, c_lengths_(c_lengths)
	, stepsize_(stepsize)
	, n_steps_min_(n_steps_min)
	, n_steps_max_(n_steps_max)
	, counter_(0)
	, counter_accepted_(0)
	, temperature_(temperature)
	, model_(model)
	{
	}


	/* reduced (!) chi2 sum divided by the global temperature used as potential energy */
	number_type potential(const Vector<number_type> & q)
	{
		return model_.potential(q)/temperature_;

	}

	/* setter for whether bounds are supposed to be fixed */
	void bounds_fixed(bool fixed)
	{
		bounds_fixed_ = fixed;
	}

	/* Setter of whether output includes analysis */
	void do_analysis (bool yes)
	{
		do_analysis_ = yes;
	}

	/* Setter for maximal chi2red value */
	void discard_from(number_type value)
	{
		chi2redmax_ = value;
		discard_data_ = true;
	}

	/* Setter for whether data above chi2redmax should be discarded */
	void discard_data(bool yes)
	{
		discard_data_ = yes;
	}

	/* check if given position is within bounds */
	bool is_in_range(Vector<number_type> position) const
	{
		bool inbounds = true;
		for (size_type i = 0; i < position.size(); ++i)
		{
			inbounds = (position[i] < range_max_[i]) && (position[i] > range_min_[i]);
			if (!inbounds)
			{
				break;
			}
		}
		return inbounds;
	}

	/* gradient of the potential */
	void grad_potential (Vector<number_type> & output, const Vector<number_type> & position)
	{
		assert(output.size() == position.size());
		for (size_type i = 0; i < output.size(); ++i)
		{
			// estimation of partial derivative with respect to q_i
			number_type h = 0.1*stepsize_*c_lengths_[i];
			Vector<number_type> position_forward = position;
			position_forward[i] += h;
			Vector<number_type> position_backward = position;
			position_backward[i] -= h;
			output[i] = (potential(position_forward)-potential(position_backward))/2.0/h;
		}

	}

	// TO DO : Move intrinsic error functions to the model class so that HMC doesn't need d_of_freedom()
	/* Calculate intrinsic uncertainty of parameter i */
	number_type intrinsic_err(const Vector<number_type> & position, size_type index)
	{
		size_type d_of_freedom = model_.d_of_freedom();
		number_type sec_deriv_chi2 = sec_deriv_potential(position, index) * d_of_freedom * temperature_;
		return sqrt(2.0/sec_deriv_chi2);
	}

	// TO DO : Move intrinsic error functions to the model class so that HMC doesn't need d_of_freedom()
	/* Calculate intrinsic uncertainty of fitting result */
	void intrinsic_err(const Vector<number_type> & position, Vector<number_type> & errors)
	{
		Matrix<number_type> Hessian(position.size(), position.size());
		Vector<number_type> jac(position.size());
		grad_potential(jac, position);
		for (size_type i = 0; i < Hessian.rowsize(); ++i)
		{
			for (size_type j = 0; j <= i; ++j) 
			{
				Hessian(i, j) = sec_deriv_potential(position, i, j);
				if (i != j)
				{
					Hessian(j, i) = Hessian(i, j); // Using symmetry of second derivatives
				}
			}
		}
		size_type d_of_freedom = model_.d_of_freedom();;
		Hessian *=  d_of_freedom * temperature_;
		
		for (size_type i = 0; i < Hessian.rowsize(); ++i)
		{
			for (size_type j = 0; j < Hessian.colsize(); ++j)
			{
				if (i == j)
					std::cout << Hessian(i, j) << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
		Hessian /= 2.0;

		Matrix<number_type> pcov(Hessian.rowsize(), Hessian.colsize());
		MatrixInverse(Hessian, pcov);
		for (size_type i = 0; i < Hessian.rowsize(); ++i)
		{
			for (size_type j = 0; j < Hessian.colsize(); ++j)
			{
				if (i == j)
					std::cout << pcov(i, j) << " ";
			}
			std::cout << "\n";
		}

		for (size_type i = 0; i < errors.size(); ++i)
		{
			errors[i] = sqrt(pcov[i][i]);
		}
	}

	/* second derivative of the potential with respect to parameter i */
	number_type sec_deriv_potential (const Vector<number_type> & position, size_type i)
	{
		assert( i >= 0 && i < position.size()); // Verify that index is not out of bounds 
		number_type h = 1e-1 * c_lengths_[i]*stepsize_;
		Vector<number_type> position_f = position;
		position_f[i] += h;
		number_type pot_f = potential(position_f);
		Vector<number_type> position_ff = position;
		position_ff[i] += 2.*h;
		number_type pot_ff = potential(position_ff);
		Vector<number_type> position_b = position;
		position_b[i] -= h;
		number_type pot_b = potential(position_b);
		Vector<number_type> position_bb = position;
		position_bb[i] -= 2.*h;
		number_type pot_bb = potential(position_bb);
		number_type pot = potential(position);

		return (-1.*pot_ff + 16.*pot_f -30.*pot + 16.*pot_b -1.*pot_bb)/12./h/h;
	}

	/* mixed second derivative of the potential with respect to parameters i ans j */
	number_type sec_deriv_potential (const Vector<number_type> & position, size_type i, size_type j)
	{
		if (i == j)
		{
			return sec_deriv_potential(position, i);
		}
		else
		{
			assert( i >= 0 && i < position.size()); // Verify that index is not out of bounds
			assert( j >= 0 && j < position.size());
			number_type hi = 1e-1 * c_lengths_[i]*stepsize_;
			number_type hj = 1e-1 * c_lengths_[j]*stepsize_;

			Vector<number_type> position_2b_2b = position; // f: forward, b: backward
			position_2b_2b[i] -= 2.*hi;
			position_2b_2b[j] -= 2.*hj;
			number_type pot_2b_2b = potential(position_2b_2b);

			Vector<number_type> position_2b_1b = position; 
			position_2b_1b[i] -= 2.*hi;
			position_2b_1b[j] -= 1.*hj;
			number_type pot_2b_1b = potential(position_2b_1b);

			Vector<number_type> position_2b_1f = position; 
			position_2b_1f[i] -= 2.*hi;
			position_2b_1f[j] += 1.*hj;
			number_type pot_2b_1f = potential(position_2b_1f);

			Vector<number_type> position_2b_2f = position; 
			position_2b_2f[i] -= 2.*hi;
			position_2b_2f[j] += 2.*hj;
			number_type pot_2b_2f = potential(position_2b_2f);

			Vector<number_type> position_1b_2b = position; 
			position_1b_2b[i] -= 1.*hi;
			position_1b_2b[j] -= 2.*hj;
			number_type pot_1b_2b = potential(position_1b_2b);

			Vector<number_type> position_1b_1b = position; 
			position_1b_1b[i] -= 1.*hi;
			position_1b_1b[j] -= 1.*hj;
			number_type pot_1b_1b = potential(position_1b_1b);

			Vector<number_type> position_1b_1f = position; 
			position_1b_1f[i] -= 1.*hi;
			position_1b_1f[j] += 1.*hj;
			number_type pot_1b_1f = potential(position_1b_1f);

			Vector<number_type> position_1b_2f = position; 
			position_1b_2f[i] -= 1.*hi;
			position_1b_2f[j] += 2.*hj;
			number_type pot_1b_2f = potential(position_1b_2f);

			Vector<number_type> position_1f_2b = position; 
			position_1f_2b[i] += 1.*hi;
			position_1f_2b[j] -= 2.*hj;
			number_type pot_1f_2b = potential(position_1f_2b);

			Vector<number_type> position_1f_1b = position; 
			position_1f_1b[i] += 1.*hi;
			position_1f_1b[j] -= 1.*hj;
			number_type pot_1f_1b = potential(position_1f_1b);

			Vector<number_type> position_1f_1f = position; 
			position_1f_1f[i] += 1.*hi;
			position_1f_1f[j] += 1.*hj;
			number_type pot_1f_1f = potential(position_1f_1f);

			Vector<number_type> position_1f_2f = position; 
			position_1f_2f[i] += 1.*hi;
			position_1f_2f[j] += 2.*hj;
			number_type pot_1f_2f = potential(position_1f_2f);

			Vector<number_type> position_2f_2b = position; 
			position_2f_2b[i] += 2.*hi;
			position_2f_2b[j] -= 2.*hj;
			number_type pot_2f_2b = potential(position_2f_2b);

			Vector<number_type> position_2f_1b = position; 
			position_2f_1b[i] += 2.*hi;
			position_2f_1b[j] -= 1.*hj;
			number_type pot_2f_1b = potential(position_2f_1b);

			Vector<number_type> position_2f_1f = position; 
			position_2f_1f[i] += 2.*hi;
			position_2f_1f[j] += 1.*hj;
			number_type pot_2f_1f = potential(position_2f_1f);

			Vector<number_type> position_2f_2f = position; 
			position_2f_2f[i] += 2.*hi;
			position_2f_2f[j] += 2.*hj;
			number_type pot_2f_2f = potential(position_2f_2f);

			number_type result = 
			(pot_2b_2b - 8. * pot_2b_1b + 8. * pot_2b_1f - pot_2b_2f) +
			(pot_1b_2b - 8. * pot_1b_1b + 8. * pot_1b_1f - pot_1b_2f) * (-8.) +
			(pot_1f_2b - 8. * pot_1f_1b + 8. * pot_1f_1f - pot_1f_2f) * 8. -
			(pot_2f_2b - 8. * pot_2f_1b + 8. * pot_2f_1f - pot_2f_2f);
			result /= (144. * hi * hi);
			

			return result;
		}
	}




	// void grad_potential_exact(Vector<number_type> & output, const Vector<number_type> & position)
	// {
	// 	output[0] = 0;
	// 	for (size_type i = 0; i < x_data_.size(); ++i)
	// 	{
	// 		output[0] += -2.0*(y_data_[i]-model_.f(x_data_[i], position))/dy_data_[i]/dy_data_[i]*1;
	// 	}
	// 	output[1] = 0;
	// 	for (size_type i = 0; i < x_data_.size(); ++i)
	// 	{
	// 		output[1] += -2.0*(y_data_[i]-model_.f(x_data_[i], position))/dy_data_[i]/dy_data_[i]*exp(-position[2]*x_data_[i]);
	// 	}
	// 	output[2] = 0;
	// 	for (size_type i = 0; i < x_data_.size(); ++i)
	// 	{
	// 		output[2] += -2.0*(y_data_[i]-model_.f(x_data_[i], position))/dy_data_[i]/dy_data_[i]*(-x_data_[i]*position[1])*exp(-position[2]*x_data_[i]);
	// 	}
	// }
	/* kinetic energy function */
	number_type kinetic (const Vector<number_type> & p)
	{
		number_type energy = 0;
		for (size_type i = 0; i < p.size(); ++i)
		{
			energy += p[i]*p[i]/2.0;
		}
		return energy;
	}

	/* Leapfrog integrator */
	void leapfrog(Vector<number_type> & q, Vector<number_type> & p)
	{
		// Introducing a step size vector: We scale the vector of characteristic length scales with the stepsize_ scalar
		Vector<number_type> stepsizes = stepsize_ * c_lengths_;
		// Make a half step for momentum at the beginning
		Vector<number_type> grad_U(q.size());
		grad_potential(grad_U, q);
		p -= stepsizes * grad_U / 2;

		// Draw random number of leapfrog steps
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_unif(n_steps_min_, n_steps_max_);
		size_type n_steps = dis_unif(gen);

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < n_steps; ++i)
		{
			// Make a full step for the position, update gradient of potential
			q += stepsizes * p;
			
			grad_potential(grad_U, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != n_steps - 1)
			{
				p -= stepsizes * grad_U;
			}
		}
		// make a half step for momentum at the end
		p -= stepsizes * grad_U / 2;
	}


	/* 
	Leapfrog integrator with a tempering parameter alpha:
	During the first half of the trajectory, there is a multiplication
	of the momenta by alpha. In the second half, correspondingly, there 
	is a division. The tempering scheme is supposed to be symmetrical
	*/
	// void leapfrog (Vector<number_type> & q, Vector<number_type> & p, number_type temperature_max)
	// {
	// 	// Draw random number of leapfrog steps
	// 	std::random_device rd;
	// 	std::mt19937 gen(rd());
	// 	std::uniform_int_distribution<> dis_unif(n_steps_min_, n_steps_max_);
	// 	size_type n_steps = dis_unif(gen);

	// 	// Adjust tempering parameter alpha according to maximal temperature to be reached
	// 	number_type alpha = pow(temperature_max, 2.0/n_steps);

	// 	// Make a half step for momentum at the beginning
	// 	Vector<number_type> grad_U(q.size());
	// 	grad_potential(grad_U, q);
	// 	p *= sqrt(alpha); //tempering
	// 	p -= stepsize_ * grad_U / 2;
		

	// 	// Alternate full steps for position and momentum
	// 	for (size_type i = 0; i < n_steps; ++i)
	// 	{
	// 		// Make a full step for the position, update gradient of potential
	// 		q += stepsize_ * p;

	// 		grad_potential(grad_U, q);
	// 		// Make a full step for the momentum, except at the end of the trajectory
	// 		if (i != n_steps - 1)
	// 		{
	// 			p -= stepsize_ * grad_U;

	// 			// tempering
	// 			if (i < n_steps/2)
	// 			{
	// 				p *= alpha;
	// 			}
	// 			else if ((n_steps%2 == 1) && (i == n_steps/2)) // special treatment of odd number of steps
	// 			{
	// 				// do nothing
	// 			}
	// 			else 
	// 			{
	// 				p /= alpha;
	// 			}
	// 		}
	// 	}
	// 	// make a half step for momentum at the end
	// 	p -= stepsize_ * grad_U / 2;
	// 	p /= sqrt(alpha); // tempering

	// }

	/* 
	Leapfrog integrator with a tempering parameter and adjustable number of leapfrog steps
	*/
	// void leapfrog (Vector<number_type> & q, Vector<number_type> & p, Storage<number_type> & q_copy, number_type temperature_max, size_type n_steps)
	// {
	// 	// Adjust tempering parameter alpha according to maximal temperature to be reached
	// 	number_type alpha = pow(temperature_max, 2.0/n_steps);

	// 	// Make a half step for momentum at the beginning
	// 	Vector<number_type> grad_U(q.size());
	// 	grad_potential(grad_U, q);
	// 	p *= sqrt(alpha); //tempering
	// 	p -= stepsize_ * grad_U / 2;
		

	// 	// Alternate full steps for position and momentum
	// 	for (size_type i = 0; i < n_steps; ++i)
	// 	{
	// 		// Make a full step for the position, update gradient of potential
	// 		q += stepsize_ * p;

	// 		// Read values of q into the storage vector
	// 		q_copy.read_in(q);

	// 		grad_potential(grad_U, q);
	// 		// Make a full step for the momentum, except at the end of the trajectory
	// 		if (i != n_steps - 1)
	// 		{
	// 			p -= stepsize_ * grad_U;

	// 			// tempering
	// 			if (i < n_steps/2)
	// 			{
	// 				p *= alpha;
	// 			}
	// 			else if ((n_steps%2 == 1) && (i == n_steps/2)) // special treatment of odd number of steps
	// 			{
	// 				// do nothing
	// 			}
	// 			else 
	// 			{
	// 				p /= alpha;
	// 			}
	// 		}
	// 	}
	// 	// make a half step for momentum at the end
	// 	p -= stepsize_ * grad_U / 2;
	// 	p /= sqrt(alpha); // tempering

	// }


	/* 
		Leapfrog integrator copying position to external vector after each step
		and offering adjustable number of iterations
	 */

	void leapfrog(Vector<number_type> & q, Vector<number_type> & p, Storage<number_type> & q_copy, size_type nb_iterations)
	{
		// Introducing a step size vector: We scale the vector of characteristic length scales with the stepsize_ scalar
		Vector<number_type> stepsizes = stepsize_ * c_lengths_;

		// Make a half step for momentum at the beginning
		Vector<number_type> grad_U(q.size());
		grad_potential(grad_U, q);
		p -= stepsizes * grad_U / 2;

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < nb_iterations; ++i)
		{
			// Make a full step for the position, update gradient of potential
			q += stepsizes * p;

			// Read values of q into the storage vector
			q_copy.read_in(q);
			

			grad_potential(grad_U, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != nb_iterations - 1)
			{
				p -= stepsizes * grad_U;
			}
		}
		// make a half step for momentum at the end
		p -= stepsizes * grad_U / 2;
	}



	/* do one Metropolis step */
	void step_forward(Vector<number_type> & current_q)
	{
		counter_++;
		Vector<number_type> q = current_q;

		// generate momenta from gaussian distribution
		std::normal_distribution<> dis_norm(0, 1); // mean 0, std deviation 1
		Vector<number_type> p(q.size());
		fill_random(p, dis_norm);
		Vector<number_type> current_p = p;

		// Compute trajectory using the Leapfrog method
		leapfrog(q, p);
		//leapfrog(q, p, temperature_);

		/* Negate momentum at the end of the trajectory to make the proposal symmetric
		(doesn't change the outcome of the algorithm, but is mathematically nicer)
		*/
		p = -p;

		// Evaluate potential (U) and kinetic (K) energies at start and end of trajectory
		number_type current_U = potential(current_q);
		number_type current_K = kinetic(current_p);
		number_type proposed_U = potential(q);
		number_type proposed_K = kinetic(p);


		/*
		Accept or reject the state at end of trajectory, returning either
		the position at the end of the trajectory or the initial position
		*/

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis_unif(0, 1);
		bool accepted = dis_unif(gen) < exp(current_U+current_K-proposed_U-proposed_K);
		if (bounds_fixed_)
		{
			accepted = accepted && is_in_range(q); // automatical rejection if proposition is out of bounds
		}
		if (accepted)
		{
			current_q = q; // accept
			counter_accepted_++;
		}
		// otherwise q is refused.



	}

	/* do n steps */
	void walk (size_type nb_steps, number_type max_duration, Vector<number_type> & initial, size_type progress_steps)
	{
		// reinitialize internal counters
		counter_ = 0;
		counter_accepted_ = 0;

		// Set timer
		std::time_t start = std::time(nullptr);
		number_type expected_duration;
		bool estimate_given = false;

		// Set up storage vector
		Storage<number_type> data(initial, 2); // records values of "initial" and two additional things

		// Start walking
		size_type counter = 0;
		while (counter < nb_steps)
		{
			step_forward(initial);

			// save updated data to storage
			data.read_in(initial); // save fitting parameters
			data.read_in(potential(initial)*temperature_); // save chi2_red
			data.read_in(acceptance_rate()); // save acceptance rate

			// After 10s: Estimate duration of the walk
			std::time_t end = std::time(nullptr);
			number_type diff = end - start;
			if ( (diff > 10) && !estimate_given)
			{
				expected_duration = nb_steps * diff / counter / 60;
				std::cout << "Computation time in min: " << expected_duration << "\n";
				estimate_given = true;
			}

			if (diff > max_duration) // Cancel calculation if time exceeds max_duration
			{
				std::cout << "Computation aborted as time ran out. " << "\n";
				break;
			}

			// Print progress
			for (size_type i = 1; i < progress_steps; ++i)
			{
				if (counter == size_type(i*nb_steps/progress_steps))
				{
					std::cout << "Progress: " << size_type(i*100./progress_steps) << "%" << "\n";
				}
			}

			counter++;
		}

		// write data to output file
		data.write("data.txt");
		
		if (discard_data_)
		{
			data.write("data_kept.txt", chi2redmax_);
		}

		// do analysis if requested
		if (do_analysis_)
		{
			// Calculate fitting result
			Vector<number_type> result(data.n_variables());
			Vector<number_type> result_err(result.size());
			data.mean(result, result_err);
			Vector<number_type> popt(initial.size());
			Vector<number_type> perr(initial.size());
			for (size_type i = 0; i < popt.size(); ++i)
			{
				popt[i] = result[i];
				perr[i] = result_err[i];
			}

			// Calculate minimal chi2_red and its error
			number_type chi2redmin = potential(popt)*temperature_;
			Vector<number_type> derivatives(popt.size());
			grad_potential(derivatives, popt);
			derivatives *= temperature_;
			 	// procedure to calculate integrated autocorrelation time
			number_type C_F = abs(data.gamma_chi2red(derivatives, 0));
			size_type W = 0;
			number_type gamma_old;
			number_type gamma_current;
			for (size_type t = 1; t <= data.entries_per_variable(); ++t)
			{
				gamma_current = abs(data.gamma_chi2red(derivatives, t));

				C_F += 2.0 * gamma_current;

				// finding optimal summation window W
				if (t>1)
				{
					if (gamma_old < gamma_current)
					{
						W = t;
						break;
					}
				}

				gamma_old = abs(gamma_current);
			}
			number_type v_F = abs(data.gamma_chi2red(derivatives, 0));

			number_type time_F = C_F/2.0/v_F;
			number_type time_F_err = time_F * sqrt(4./data.entries_per_variable()*(W+1./2- time_F )) + C_F*exp(-1.*W/time_F);
			
			std::cout << "Integrated autocorrelation time for chi2_red_min: "
			<< time_F << " + - " << time_F_err << "\n";
			number_type ESS = data.entries_per_variable()/(2.*time_F); // effective sample size
			number_type chi2redmin_err = sqrt(v_F/ESS);
		

			// report total calculation time
			std::time_t end = std::time(nullptr);
			number_type diff = (end - start)/60;
			std::cout << "Total calculation time in min : " << diff << "\n";

			// report burn-in length
			std::cout << "Burn-in length: " << data.get_burn_in_length() << "\n";

			// report fitting result
			std::cout << "FITTING RESULT: \n";
			for (size_type i = 0; i < initial.size(); ++i)
			{
				std::cout << "Parameter " << i << " : " << popt[i] << " + - " << perr[i]  << "\n";
			}

			
			number_type chi2redmean = result[popt.size()];
			number_type chi2reddiff = chi2redmean - chi2redmin;
			number_type chi2reddiff_theo = initial.size() * temperature_ / 2;
			std::cout << "chi2_red (best fit): " << chi2redmin << " + - " << chi2redmin_err << "\n";
			std::cout << "chi2_red (mean): " << chi2redmean << " + - " << result_err[popt.size()] << "\n";
			std::cout << "Ratio of measured difference to theoretical difference: " 
			<< chi2reddiff/chi2reddiff_theo << "\n";
			
			// report intrinsic uncertainties
			std::cout << "Intrinsic uncertainties: \n";
			Vector<number_type> err_intr(popt.size());
			intrinsic_err(popt, err_intr);
			for (size_type i = 0; i < err_intr.size(); ++i)
			{
				std::cout << "Parameter " << i << " : " << err_intr[i]  << "\n";
			}
		}
	}

	

	number_type acceptance_rate()
	{
		return 1.0 * counter_accepted_ / counter_;
	}

	/* preliminary run to estimate correct step size */
	void get_acceptance_rates(const Vector<number_type> & range_min, const Vector<number_type> & range_max, size_type nb_positions, size_type nb_leapfrog_steps, std::string filename)
	{
		std::vector<number_type> acceptance_rates(0);

		for (size_type i = 0; i<nb_positions; ++i)
		{
			// Draw random position from the search region
			Vector<number_type> popt(range_min.size());
			fill_from_region(popt, range_min, range_max);
			
			// do some leapfrog steps and save acceptance rate
			counter_ = 0;
			counter_accepted_ = 0;
			for (size_type j = 0; j < nb_leapfrog_steps; ++j)
			{
				step_forward(popt);
			}
			acceptance_rates.push_back(acceptance_rate());
		}

		// Save acceptance rate vector in a data file
		std::ofstream outfile(filename);
		for (size_type i = 0; i < acceptance_rates.size(); ++i)
		{
			outfile << acceptance_rates[i] << "\n";
		}

		
	}

	/* preliminary run to estimate number of leapfrog steps */
	void get_optimal_number_of_steps(const Vector<number_type> & range_min, const Vector<number_type> & range_max, size_type nb_positions, size_type nb_iterations, std::string filename)
	{
		std::vector<number_type> number_steps(0);

		for (size_type i = 0; i<nb_positions; ++i)
		{
			// Draw random position from the search region
			Vector<number_type> popt(range_min.size());
			fill_from_region(popt, range_min, range_max);
			
			// do one leapfrog step and save position values for each iteration
			Storage<number_type> q_values(popt);
			std::random_device rd;
			std::mt19937 gen(rd());
			std::normal_distribution<> dis_norm(0, 1);
			Vector<number_type> p(popt.size());
			fill_random(p, dis_norm);
			leapfrog(popt, p, q_values, nb_iterations);

			// calculate period of each position component using autocorrelation
			Vector<size_type> period_vector(popt.size());
			q_values.period(period_vector);
	

			// Half a period corresponds to the optimal number of steps
			for (size_type j = 0; j < period_vector.size(); ++j)
			{
				number_steps.push_back(period_vector[j]/2);
			}
			
		}

		// Save array of optimal number of steps in a data file
		std::ofstream outfile(filename);
		for (size_type i = 0; i < number_steps.size(); ++i)
		{
			outfile << number_steps[i] << "\n";
		}	
	}


private:
	/* measured data */
	Model<number_type> model_;
	Vector<number_type> range_min_;
	Vector<number_type> range_max_;
	bool bounds_fixed_;
	bool do_analysis_;
	bool discard_data_;
	number_type chi2redmax_;

	/* parameters for the leapfrog integrator */
	number_type stepsize_;
	size_type n_steps_min_;
	size_type n_steps_max_;
	number_type temperature_;
	Vector<number_type> c_lengths_; // characteristic length scales for parameters

	/* some statistics */
	size_type counter_;
	size_type counter_accepted_;

};

#endif