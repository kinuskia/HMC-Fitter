# HMC-Fitter

In the following, I would like to explain briefly how the methods of the Hybrid Monte Carlo (HMC) class are used to fit a function to experimental data. 

The HMC method is an algorithm commonly used to sample from a Boltzmann distribution. The idea is therefore to define a Boltzmann distribution with the reduced chi2 as an energy function. We will keep the inverse temperature beta characterizing a canonical ensemble as a free parameter of the Boltzmann distribution to be tuned.

The folder is structured as follows: 
1. The file "fitter.cpp" contains the main function. Preliminary runs as well as the actual run to generate a Markov chain for the fitting result are controlled from here.
2. The fitting function is defined in "model.hpp". The user sets the number of fitting parameters, degrees of freedom, the fitting function, the reduced chi2 sum and possible constraints.
3. The folder "auxiliary_files" contains the code for the HMC class as well as for the auxiliary classes needed in this project. For instance,
I wrote some useful classes to deal with Euklidian vectors and matrices. 
Especially the vector class is broadly used throughout this project. It is inherited from the C++ STL vector class. For instance, one
initializes a vector with 14 double entries as follows:
	
	Vector<double> q(14);

4. The folder "preliminary tools", contains certain auxiliary text files originating from preliminary runs are saved.
5. The python script "plot-fit.py" can be used to plot histograms from Markov chain data for each parameter.

I would like to explain the fitting procedure with a concrete example. The folder "experimental_data" contains a text file with data points (time, y-value, y-uncertainty). We shall try to fit a gaussian model to the data:

	f(t) = N exp(-(x-mu)^2/2/sigma^2)

The fitting function thus has three parameters. This specific problem has already been defined in the "model.hpp" file, so we can proceed with the "fitter.cpp" file. First, the experimental data is read in and saved in specific vectors ("read_data"). Next, the fitting model is initialized with the experimental data ("gaussian"). After the initialization of the fitting parameter vector `popt`, a search region for each parameter as well as characteristic length scales are specified.

Next, a HMC object is initialized which I would like to call "sampler". 
Its declaration reads:
	
	HMC<double> sampler(Model<double> model, Vector<double> range_min, Vector<double> range_max, Vector<double> c_lengths, double stepsize, std::size_t n_steps_min, std::size_t n_steps_max, double beta);

It takes 8 arguments:
1. `Model<double> model`: A class which defines the chi2-function to be minimized.
2. `Vector<double> range_min`: Lower range for the components of the starting point.
3. `Vector<double> range_max`: Upper range for the components of the stating point.
4. `Vector<double> c_lengths`: characteristic length scales for each parameter, typically range_max - range_min, but if one has additional information, it can be put to use here.
5. `double stepsize`: The step size to be used by the Leapfrog integrator
6. `std::size_t n_steps_min`: The number of Leapfrog steps in each iteration is drawn randomly from values, the lower bound of which is defined here
7. `std::size_t n_steps_max`: Definition of the upper bound for the number of Leapfrog steps
8. `double beta`: Set the parameter beta of the system (inverse of a temperature)

Before generating the first Markov chain, one has to tune the parameters of the algorithm, namely the Leapfrog step size and the number of Leapfrog steps as well as the parameter beta.
The method

	sampler.tune_parameters();

executes a commmando-based routine which helps setting the HMC parameters. The routine first prints log(potential) evaluated in the search region to get an idea for beta. Potential times beta should amount to approximately 100. For the given search region, beta=1e0 should work fine.
Next, the routine asks for a trial step size and calls a method which creates 50 Markov chains of length 50 with random starting points from the search region and saves the final acceptance rate for each chain in a file `acceptrates.txt`. The data in this file may be plotted in a histogram with the python script "preliminary_tools/plot_accept.py". The overall acceptance rate should be around 90 percent. The choice of 3e-3 as step size should work. Next, one is asked to fix lower and upper bounds for the number of leapfrog steps. Usually Lmin=90 and Lmax=130 are a good start. A a later step, integrated autocorrelation times are computed. If they were to be unusually high (>>0.5), the bounds should be increased. After a certain number of chains (e.g. 15) has been specified, the method above again computes some acceptance rates, which may be plotted in a histogram to verify that the overall acceptance rate is still around 90 percent.
In the next step, the routine generates a Markov chain of length 200, which can be made longer as explained in the routine. The resulting Markov chain is saved to "data.txt" which has the following columns: chain index, parameters (3 in this case), reduced chi2 sum, momentary acceptance rate. The python script "plot-fit.py" may be used to plot histograms for each parameter and chi2red as well as to plot the variation of each parameter with respect to the chain index. One should verify that the autocorrelation times (they are also computed. If they show "nan", the chain might be too short) and the burn-in length of the chain are small. 



Once the HMC parameters have been tuned, the initialization of the HMC object can be updated with the correct parameter values. 
The tune_parameters method should be commented.

Now we can generate a Markov chain of a specific length, e.g. 1e4:

	fill_from_region(popt, range_min, range_max);
	sampler.walk(1e4, 60*30, popt, 5);

The starting point (random, generated by `fill_from_region`) is passed by reference.
The calculation is aborted either after the Markov chain has grown to the set length or once `60*30` seconds have passed. The fitting result and uncertainties are calculated automatically and printed to the console. The last argument constitutes the number of progress steps: In this case there is an output once 20, 40, 60, 80 and 100 percent of the chain have been generated.

The resulting Markov chain is again saved as "data.txt". The best-fit values of the fitting parameters as well as their uncertainties are printed to the console, together with a comparison of theoretical expected values and variance. If theory and experiment differ considerably from one another, one should use histogram plots to update the search region and tune the parameters again at a higher beta. 

Note: If the walk method ends with a `segmentation fault: 11`, it should be run once again. This bug is caused by an extremely long burn-in length which can happen occasionally.

	
