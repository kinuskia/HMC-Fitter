#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include <random>
#include <assert.h>
#include <fstream>

typedef double number_type;
typedef std::size_t size_type;


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
	}

}







int main ()
{
	Vector<double> t(0);
	Vector<double> Pl_Pl(0);
	Vector<double> d_Pl_Pl(0);

	read_data("../correlators/correlator_Pl-Pl", t, Pl_Pl, d_Pl_Pl);

	for (size_type i = 0; i < t.size(); ++i)
	{
		std::cout << t[i] << " " << Pl_Pl[i] << " " << d_Pl_Pl[i] << "\n";
	}
	

	return 0;
}