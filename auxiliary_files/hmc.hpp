#ifndef HMC_HPP
#define HMC_HPP

#include "vector.hpp"
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
	HMC(Vector<number_type> x_data, Vector<number_type> y_data, Vector<number_type> dy_data, number_type stepsize, size_type n_steps_min, size_type n_steps_max, number_type temperature)
	: x_data_(x_data)
	, y_data_(y_data)
	, dy_data_(dy_data) 
	, stepsize_(stepsize)
	, n_steps_min_(n_steps_min)
	, n_steps_max_(n_steps_max)
	, counter_(0)
	, counter_accepted_(0)
	, temperature_(temperature)
	, model_()
	{
	}


	/* reduced (!) chi2 sum divided by the global temperature used as potential energy */
	number_type potential(const Vector<number_type> & q)
	{
		/* checks if data arrays have equal lengths */
		assert(x_data_.size() == y_data_.size());
		assert(x_data_.size() == dy_data_.size());
		number_type chi2 = 0;
		for (size_type i = 0; i<x_data_.size(); ++i)
		{
			chi2 += (y_data_[i] - model_.f(x_data_[i], q))*(y_data_[i] - model_.f(x_data_[i], q))/dy_data_[i]/dy_data_[i];
		}

		size_type d_of_freedom = x_data_.size() - q.size();
		return chi2/d_of_freedom/temperature_;

	}

	/* gradient of the potential */
	void grad_potential (Vector<number_type> & output, const Vector<number_type> & position)
	{
		assert(output.size() == position.size());
		for (size_type i = 0; i < output.size(); ++i)
		{
			// estimation of partial derivative with respect to q_i
			Vector<number_type> position_forward = position;
			position_forward[i] += stepsize_;
			Vector<number_type> position_backward = position;
			position_backward[i] -= stepsize_;
			output[i] = (potential(position_forward)-potential(position_backward))/2.0/stepsize_;
		}

	}
	/* Calculate intrinsic uncertainty of parameter i */
	number_type intrinsic_err(const Vector<number_type> & position, size_type index)
	{
		size_type d_of_freedom = x_data_.size() - position.size();
		number_type sec_deriv_chi2 = sec_deriv_potential(position, index) * d_of_freedom * temperature_;
		return sqrt(2.0/sec_deriv_chi2);
	}

	/* second derivative of the potential with respect to parameter i */
	number_type sec_deriv_potential (const Vector<number_type> & position, size_type i)
	{
		number_type h = 0.1*stepsize_;
		Vector<number_type> position_forward = position;
		position_forward[i] += h;
		Vector<number_type> position_backward = position;
		position_backward[i] -= h;

		return (potential(position_backward)-2.0*potential(position)+potential(position_forward))/h/h;
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
		// Make a half step for momentum at the beginning
		Vector<number_type> grad_U(q.size());
		grad_potential(grad_U, q);
		p -= stepsize_ * grad_U / 2;

		// Draw random number of leapfrog steps
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_unif(n_steps_min_, n_steps_max_);
		size_type n_steps = dis_unif(gen);

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < n_steps; ++i)
		{
			// Make a full step for the position, update gradient of potential
			q += stepsize_ * p;
			
			grad_potential(grad_U, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != n_steps - 1)
			{
				p -= stepsize_ * grad_U;
			}
		}
		// make a half step for momentum at the end
		p -= stepsize_ * grad_U / 2;
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
		// Make a half step for momentum at the beginning
		Vector<number_type> grad_U(q.size());
		grad_potential(grad_U, q);
		p -= stepsize_ * grad_U / 2;

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < nb_iterations; ++i)
		{
			// Make a full step for the position, update gradient of potential
			q += stepsize_ * p;

			// Read values of q into the storage vector
			q_copy.read_in(q);
			

			grad_potential(grad_U, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != nb_iterations - 1)
			{
				p -= stepsize_ * grad_U;
			}
		}
		// make a half step for momentum at the end
		p -= stepsize_ * grad_U / 2;
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
		if (dis_unif(gen) < exp(current_U+current_K-proposed_U-proposed_K))
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

		// Calculate fitting result
		Vector<number_type> popt(initial.size());
		Vector<number_type> perr(initial.size());
		data.mean(popt, perr);
	

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

		number_type chi2redmin = potential(popt)*temperature_;
		number_type chi2redmean = data.mean(5);
		number_type chi2reddiff = chi2redmean - chi2redmin;
		number_type chi2reddiff_theo = initial.size() * temperature_ / 2;
		std::cout << "chi2_red (best fit): " << chi2redmin << "\n";
		std::cout << "chi2_red (mean): " << chi2redmean << "\n";
		std::cout << "Ratio of measured difference to theoretical difference: " 
		<< chi2reddiff/chi2reddiff_theo << "\n";
		// report intrinsic uncertainties
		std::cout << "Intrinsic uncertainties: \n";
		for (size_type i = 0; i < initial.size(); ++i)
		{
			std::cout << "Parameter " << i << " : " << intrinsic_err(popt, i)  << "\n";
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
	Vector<number_type> x_data_;
	Vector<number_type> y_data_;
	Vector<number_type> dy_data_;
	Model<number_type> model_;

	/* parameters for the leapfrog integrator */
	number_type stepsize_;
	size_type n_steps_min_;
	size_type n_steps_max_;
	number_type temperature_;

	/* some statistics */
	size_type counter_;
	size_type counter_accepted_;

};

#endif