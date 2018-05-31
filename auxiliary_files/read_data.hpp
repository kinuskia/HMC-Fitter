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

void read_data(std::string filename, Vector<double> &x_data, Vector<double> &y_data, Vector<double> &dy_data, std::size_t skip_rows = 0)
{
	x_data.clear();
	y_data.clear();
	dy_data.clear();
	std::ifstream infile(filename);
	std::string line;

	std::size_t counter_numbers = 0;
	std::size_t counter_lines = 0;
	while (infile)
	{
		std::getline(infile, line); // Read in current line
		if (line == "")
		{
			continue;  // ignore empty lines
		}

		if (counter_lines < skip_rows) // skip first rows if one wished so
		{
			counter_lines++;
			continue;
		}

		bool end_of_number = true;
		std::string numberstring = "";
		for (int i = 0; i < line.length(); ++i)
		{
			// Detect beginning of a number
			bool is_number = std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-');
			bool is_seperator = (line[i] == ' ' && line[i] == '\t');
			if (is_number) // Detect beginning of a number
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
		counter_lines++;
	}

}

#endif