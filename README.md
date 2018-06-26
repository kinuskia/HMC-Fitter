# HMC-Fitter

In the following, I would like to explain briefly how the methods of the Hybrid Monte Carlo (HMC) class are used to find the global minimum of a function. They are defined in the folder "auxiliary_files". For this project,
I wrote some useful classes to deal with Euklidian vectors and matrices. 
Especially the vector class is broadly used throughout this project. It is inherited from the C++ STL vector class. For instance, one
initializes a vector with 14 double entries as follows:
	Vector<double> q(14);

The central object is a HMC type which I would like to call "sampler". It provides methods to generate Markov chains following a Boltzmann distribution with the function to be minimized as energy function, divided by a (hypothetical) temperature which can be set manually. The starting point of a Markov chain is drawn randomly from a range that can be set as well 
It is initialized as follows:
	
	HMC<double> sampler(Model<double> model, Vector<double> range_min, Vector<double> range_max, Vector<double> c_lengths, double stepsize, std::size_t n_steps_min, std::size_t n_steps_max, double temperature);

It takes 8 arguments:
1. `Model<double> model`: A class which defining the function to be minimized, e.g. "model.hpp". The class has to be called `Model<type>` and contain at least the following methods: return the number of parameters of the problem (-> `size_type n_parameters() const`), the function to be minimized (-> number_type potential(`const Vector<number_type> &q)`) and a method to check if possible constraints on the parameters one has set are respected at position q (-> `bool constraints_respected(const Vector<number_type> &q)`).
2. `Vector<double> range_min`: Lower range for the components of the starting point.
3. `Vector<double> range_max`: Upper range for the components of the stating point.
4. `Vector<double> c_lengths`: characteristic length scales for each parameter, typically range_max - range_min, but if one has additional information, it can be put to use here.
5. `double stepsize`: The step size to be used by the Leapfrog integrator
6. `std::size_t n_steps_min`: The number of Leapfrog steps in each iteration is drawn randomly from values, the lower bound of which is defined here
7. `std::size_t n_steps_max`: Definition of the upper bound for the number of Leapfrog steps
8. `double temperature`: Set the temperature of the system

Before generating the first Markov chain, one has to tune the parameters of the algorithm, namely the Leapfrog step size and the number of Leapfrog steps. The following method can be used to find the optimal step size. It creates 70 Markov chains of length 50 with random starting points from the search region and saves the final acceptance rate for each chain in a file `acceptrates.txt`:

	sampler.get_acceptance_rates(range_min, range_max, 70, 50, "acceptrates.txt");

The next method is designed to find the optimal number of leapfrog steps. It draws randomly 300 starting points and lets evolve the system over 700 Leapfrog steps. After that, the values are analyzed for periodicities. Half of a period is considered an optimal number of step an is saved to the output file `correlation_times`.

	sampler.get_optimal_number_of_steps(range_min, range_max, 300, 700, "correlation_times.txt")

Once the HMC parameters have been tuned, a Markov chain of length 5e4 is generated as follows:

	sampler.walk(5e4, 60*30, popt, 5);

The starting point is passed by reference.
The calculation is aborted either after the Markov chain has grown to the set length or once `60*30` seconds have passed. The fitting result and uncertainties are calculated automatically and printed to the console. The last argument constitutes the number of progress steps: In this case there is an output once 20, 40, 60, 80 and 100 percent of the chain have been generated.

	
