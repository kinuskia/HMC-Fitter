#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include <random>
#include <assert.h>

typedef double number_type;
typedef std::size_t size_type;

template<typename REAL>
class Playground
{
public:
	typedef double number_type;
	typedef std::size_t size_type;
	Playground(Vector<number_type> c_lengths, number_type stepsize)
	: c_lengths_(c_lengths)
	, stepsize_(stepsize)
	{}

	number_type potential (Vector<number_type> x)
{
	number_type result = -cos(sqrt((2*x[0]-x[1])*(2*x[0]-x[1])+(2*x[0]+x[1])*(2*x[0]+x[1])));
	return result;
}

/* second derivative of the potential with respect to parameter i */
	number_type sec_deriv_potential (const Vector<number_type> & position, size_type i)
	{
		assert( i >= 0 && i < position.size()); // Verify that index is not out of bounds 
		number_type h = 1e0 * c_lengths_[i]*stepsize_;
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

		number_type result = (-1.*pot_ff + 16.*pot_f -30.*pot + 16.*pot_b -1.*pot_bb)/12./h/h;
		std::cout << result << "\n";
		return result;
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
			number_type hi = 1e-7 * c_lengths_[i]*stepsize_;
			number_type hj = 1e-7 * c_lengths_[j]*stepsize_;

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
			std::cout << result << "\n";

			return result;
		}
	}

	/* Calculate intrinsic uncertainty of fitting result */
	void intrinsic_err(const Vector<number_type> & position, Vector<number_type> & errors)
	{
		Matrix<number_type> Hessian(position.size(), position.size());
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
		for (size_type i = 0; i < Hessian.rowsize(); ++i)
		{
			for (size_type j = 0; j < Hessian.colsize(); ++j)
			{
				std::cout << Hessian(i, j) << " ";
			}
			std::cout << "\n";
		}
		
		//size_type d_of_freedom = x_data_.size() - position.size();
		//Hessian *=  d_of_freedom * temperature_;
		Hessian /= 2.0;

		Matrix<number_type> pcov(Hessian.rowsize(), Hessian.colsize());
		MatrixInverse(Hessian, pcov);
		for (size_type i = 0; i < Hessian.rowsize(); ++i)
		{
			for (size_type j = 0; j < Hessian.colsize(); ++j)
			{
				std::cout << pcov(i, j) << " ";
			}
			std::cout << "\n";
		}

		for (size_type i = 0; i < errors.size(); ++i)
		{
			errors[i] = sqrt(pcov[i][i]);
		}
	}



private:
	Vector<number_type> c_lengths_;
	number_type stepsize_;
};








int main ()
{
	Vector<number_type> popt(2,0);
	Vector<number_type> c_lengths(2, 1);
	number_type stepsize = 1e-3;

	Playground<number_type> playground(c_lengths, stepsize);

	
	Vector<number_type> perr(2);
	//playground.intrinsic_err(popt, perr);

	Matrix<number_type> H(3,3);
	H(0,0) = 1;
	H(0,1) = 1;
	H(0,2) = 2;
	H(1,0) = 1;
	H(1,1) = 2;
	H(1,2) = 1;
	H(2,0) = 2;
	H(2,1) = 1;
	H(2,2) = 6;

	Matrix<number_type> H_inv = H;

	MatrixInverse_cholesky(H, H_inv);
	
	
	// for (size_type i = 0; i < popt.size(); ++i)
	// 	std::cout << popt[i] << "\n";

	// for (size_type i = 0; i < perr.size(); ++i)
	// 	std::cout << perr[i] << "\n";

	for (size_type i = 0; i < H.colsize(); ++i)
	{
		for (size_type j = 0; j < H.rowsize(); ++j)
		{
			std::cout << H_inv(i,j) << " ";
		}
		std::cout << "\n";
	}



	return 0;
}