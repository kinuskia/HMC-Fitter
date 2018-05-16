#ifndef READ_DATA_HPP
#define READ_DATA_HPP

#include <fstream>
#include "vector.hpp"

/* free functions to read in a data file and save it in a vector */

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

#endif