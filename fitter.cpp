#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "vector.hpp"
#include "hmc.hpp"
#include <fstream>

void insert(double number, std::size_t counter, Vector<double> &x_data, Vector<double> &y_data, Vector<double> &dy_data)
{
	if (counter%3 == 0)
		x_data.push_back(number);
	if (counter%3 == 1)
		y_data.push_back(number);
	if (counter%3 == 2)
		dy_data.push_back(number);
}

void read_data(std::string filename, Vector<double> &x_data, Vector<double> &y_data, Vector<double> &dy_data)
{
	std::ifstream infile(filename);
	std::string line;

	std::size_t counter_numbers = 0;
	while (infile)
	{
		std::getline(infile, line); // Read in current line
		if (line == "")
		{
			continue;  // ignore empty lines
		}
		bool end_of_number = false;
		std::string numberstring = "";
		for (int i = 0; i < line.length(); ++i)
		{
			if (line[i] != ' ') // Read in until spaces
			{
				end_of_number = false;
				numberstring += line[i];
			}
			else if (!end_of_number)
			{
				end_of_number = true;
				double number = std::stod(numberstring);
				insert(number, counter_numbers, x_data, y_data, dy_data); // save number in the correct vector
				counter_numbers++;
				numberstring = "";
			}
			else
			{
				continue; // ignore second, third ... space in a series of spaces
			}
		}
		double number = std::stod(numberstring); // end of line finishes number
		insert(number, counter_numbers, x_data, y_data, dy_data);
		counter_numbers++;
		numberstring = "";
	}

}

int main ()
{
	Vector<double> x_data(0);
	Vector<double> y_data(0);
	Vector<double> dy_data(0);

	// Read in measured data
	read_data("Measurement_data_generator/measurements.txt", x_data, y_data, dy_data);

	// Vector for fitting parameters
	Vector<double> popt(3);

	// Estimated search region
	Vector<double> range_min(popt.size());
	Vector<double> range_max(popt.size());
	range_min[0] = 0;
	range_max[0] = 6;
	range_min[1] = -1;
	range_max[1] = 5;
	range_min[2] = 0;
	range_max[2] = 2;
	// range_min[3] = -0.5;
	// range_max[3] = 1.5;

	//initialize HMC opbject
	HMC<double> sampler(x_data, y_data, dy_data, 9e-3, 15);
	
	/* PRELIMINARY RUN TOOLS */
	// draws positions and returns acceptance rates in a file (to adjust leapfrog step size)
	//sampler.get_acceptance_rates(range_min, range_max, 500, 50, "acceptrates.txt");
	


	/* ACTUAL RUN */

	// initial guess for fitting variables : random pick from region above
	fill_from_region(popt, range_min, range_max);
	sampler.walk(3e6, 60*8, popt, 10);
	

	return 0;
}