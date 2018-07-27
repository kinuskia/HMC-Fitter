#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include <assert.h>
#include "auxiliary_files/vector.hpp"

template<typename REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Model(
	Vector<number_type> t,
	Vector<number_type> Vl_Vl,
	Vector<number_type> Vl_Vs1,
	Vector<number_type> Vl_Vs2,
	Vector<number_type> Vl_Vs3,
	Vector<number_type> Vs1_Vl,
	Vector<number_type> Vs1_Vs1,
	Vector<number_type> Vs1_Vs2,
	Vector<number_type> Vs1_Vs3,
	Vector<number_type> Vs2_Vl,
	Vector<number_type> Vs2_Vs1,
	Vector<number_type> Vs2_Vs2,
	Vector<number_type> Vs2_Vs3,
	Vector<number_type> Vs3_Vl,
	Vector<number_type> Vs3_Vs1,
	Vector<number_type> Vs3_Vs2,
	Vector<number_type> Vs3_Vs3,
	Vector<number_type> d_Vl_Vl,
	Vector<number_type> d_Vl_Vs1,
	Vector<number_type> d_Vl_Vs2,
	Vector<number_type> d_Vl_Vs3,
	Vector<number_type> d_Vs1_Vl,
	Vector<number_type> d_Vs1_Vs1,
	Vector<number_type> d_Vs1_Vs2,
	Vector<number_type> d_Vs1_Vs3,
	Vector<number_type> d_Vs2_Vl,
	Vector<number_type> d_Vs2_Vs1,
	Vector<number_type> d_Vs2_Vs2,
	Vector<number_type> d_Vs2_Vs3,
	Vector<number_type> d_Vs3_Vl,
	Vector<number_type> d_Vs3_Vs1,
	Vector<number_type> d_Vs3_Vs2,
	Vector<number_type> d_Vs3_Vs3
	) // default constructor
	: t_(t)
	, Vl_Vl_(Vl_Vl)
	, Vl_Vs1_(Vl_Vs1)
	, Vl_Vs2_(Vl_Vs2)
	, Vl_Vs3_(Vl_Vs3)
	, Vs1_Vl_(Vs1_Vl)
	, Vs1_Vs1_(Vs1_Vs1)
	, Vs1_Vs2_(Vs1_Vs2)
	, Vs1_Vs3_(Vs1_Vs3)
	, Vs2_Vl_(Vs2_Vl)
	, Vs2_Vs1_(Vs2_Vs1)
	, Vs2_Vs2_(Vs2_Vs2)
	, Vs2_Vs3_(Vs2_Vs3)
	, Vs3_Vl_(Vs3_Vl)
	, Vs3_Vs1_(Vs3_Vs1)
	, Vs3_Vs2_(Vs3_Vs2)
	, Vs3_Vs3_(Vs3_Vs3)
	, d_Vl_Vl_(d_Vl_Vl)
	, d_Vl_Vs1_(d_Vl_Vs1)
	, d_Vl_Vs2_(d_Vl_Vs2)
	, d_Vl_Vs3_(d_Vl_Vs3)
	, d_Vs1_Vl_(d_Vs1_Vl)
	, d_Vs1_Vs1_(d_Vs1_Vs1)
	, d_Vs1_Vs2_(d_Vs1_Vs2)
	, d_Vs1_Vs3_(d_Vs1_Vs3)
	, d_Vs2_Vl_(d_Vs2_Vl)
	, d_Vs2_Vs1_(d_Vs2_Vs1)
	, d_Vs2_Vs2_(d_Vs2_Vs2)
	, d_Vs2_Vs3_(d_Vs2_Vs3)
	, d_Vs3_Vl_(d_Vs3_Vl)
	, d_Vs3_Vs1_(d_Vs3_Vs1)
	, d_Vs3_Vs2_(d_Vs3_Vs2)
	, d_Vs3_Vs3_(d_Vs3_Vs3)
	{}

	// /* print content of three vectors to console */
	// void print_content() const
	// {
	// 	Vector<number_type> x = Vs3_Vs3_;
	// 	Vector<number_type> dx = d_Vs3_Vs3_;
	// 	assert(t_.size() == x.size());
	// 	assert(x.size() == dx.size());
	// 	for (size_type i = 0; i < t_.size(); ++i)
	// 	{
	// 		std::cout << t_[i] << "\t" << x[i] << " +/- " << dx[i] <<"\n"; 
	// 	}

	// }

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 12;
	}
public: // needs to become private once I focus on intrinsic errors
	/* degrees of freedom */
	size_type d_of_freedom() const
	{
		return t_.size() * 2 * 2 - n_parameters();
	}
private:
	/* fitting model for correlator C_ij(t) */
	/*
	nomenclature for correlators
	index	correlator
	0		Vl
	1		Vs1

	
	
	nomenclature for fitting parameters
	index	parameter
	0		mass ground state
	1 		mass of first excited state
	2		Z0 of m1
	3		Z0 of m2
	4		Z1 of m1
	5		Z1 of m2
	6		A_00
	7		A_01
	8		A_11
	9		m_00
	10		m_01
	11		m_11
	*/
	number_type C(number_type t, size_type i, size_type j, const Vector<number_type> & popt)
	{
		assert(n_parameters() == popt.size());
		number_type T = 64.;
		number_type m1 = popt[0];
		number_type m2 = popt[1];
		number_type Zl1;
		number_type Zl2;
		number_type Zr1;
		number_type Zr2;
		number_type Ar;
		number_type mr;
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			// Coupling constants of mass 1
			Zl1 = popt[2];
			Zr1 = popt[2];
			// Coupling constants of mass 2
			Zl2 = popt[3];
			Zr2 = popt[3];
			// Residual coupling constant
			Ar = popt[6];
			//Residual mass
			mr = popt[9];

			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;
		}

		if ((i == 0 && j == 1) || (i == 1 && j == 0))
		{
			// Coupling constants of mass 1
			Zl1 = popt[2];
			Zr1 = popt[4];
			// Coupling constants of mass 2
			Zl2 = popt[3];
			Zr2 = popt[5];
			// Residual coupling constant
			Ar = popt[7];
			//Residual mass
			mr = popt[10];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;	
		}

		// if ((i == 0 && j == 2) || (i == 2 && j == 0))
		// {
		// 	// Coupling constants of mass 1
		// 	Zl1 = popt[2];
		// 	Zr1 = popt[6];
		// 	// Coupling constants of mass 2
		// 	Zl2 = popt[3];
		// 	Zr2 = popt[7];
		// 	// Residual coupling constant
		// 	Ar = popt[12];
		// 	//Residual mass
		// 	mr = popt[22];

		// 	// result
		// 	result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
		// 	result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
		// 	result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
		// 	return result;	
		// }

		// if ((i == 0 && j == 3) || (i == 3 && j == 0))
		// {
		// 	// Coupling constants of mass 1
		// 	Zl1 = popt[2];
		// 	Zr1 = popt[8];
		// 	// Coupling constants of mass 2
		// 	Zl2 = popt[3];
		// 	Zr2 = popt[9];
		// 	// Residual coupling constant
		// 	Ar = popt[13];
		// 	//Residual mass
		// 	mr = popt[23];
		// 	// result
		// 	result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
		// 	result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
		// 	result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
		// 	return result;	
		// }

		
		if (i == 1 && j == 1)
		{
			// Coupling constants of mass 1
			Zl1 = popt[4];
			Zr1 = popt[4];
			// Coupling constants of mass 2
			Zl2 = popt[5];
			Zr2 = popt[5];
			// Residual coupling constant
			Ar = popt[8];
			//Residual mass
			mr = popt[11];	
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;
		}

		// if ((i == 1 && j == 2) || (i == 2 && j == 1))
		// {
		// 	// Coupling constants of mass 1
		// 	Zl1 = popt[4];
		// 	Zr1 = popt[6];
		// 	// Coupling constants of mass 2
		// 	Zl2 = popt[5];
		// 	Zr2 = popt[7];
		// 	// Residual coupling constant
		// 	Ar = popt[15];
		// 	//Residual mass
		// 	mr = popt[25];
		// 	// result
		// 	result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
		// 	result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
		// 	result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
		// 	return result;	
		// }

	// 	if ((i == 1 && j == 3) || (i == 3 && j == 1))
	// 	{
	// 		// Coupling constants of mass 1
	// 		Zl1 = popt[4];
	// 		Zr1 = popt[8];
	// 		// Coupling constants of mass 2
	// 		Zl2 = popt[5];
	// 		Zr2 = popt[9];
	// 		// Residual coupling constant
	// 		Ar = popt[16];
	// 		//Residual mass
	// 		mr = popt[26];
	// 		// result
	// 		result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
	// 		result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
	// 		result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
	// 		return result;	
	// 	}

	// 	if (i == 2 && j == 2)
	// 	{
	// 		// Coupling constants of mass 1
	// 		Zl1 = popt[6];
	// 		Zr1 = popt[6];
	// 		// Coupling constants of mass 2
	// 		Zl2 = popt[7];
	// 		Zr2 = popt[7];
	// 		// Residual coupling constant
	// 		Ar = popt[17];
	// 		//Residual mass
	// 		mr = popt[27];	
		
	// 		// result
	// 		result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
	// 		result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
	// 		result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
	// 		return result;
	// 	}

	// 	if ((i == 2 && j == 3) || (i == 3 && j == 2))
	// 	{
	// 		// Coupling constants of mass 1
	// 		Zl1 = popt[6];
	// 		Zr1 = popt[8];
	// 		// Coupling constants of mass 2
	// 		Zl2 = popt[7];
	// 		Zr2 = popt[9];
	// 		// Residual coupling constant
	// 		Ar = popt[18];
	// 		//Residual mass
	// 		mr = popt[28];
	// 		// result
	// 		result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
	// 		result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
	// 		result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
	// 		return result;	
	// 	}

	// 	if (i == 3 && j == 3)
	// 	{
	// 		// Coupling constants of mass 1
	// 		Zl1 = popt[8];
	// 		Zr1 = popt[8];
	// 		// Coupling constants of mass 2
	// 		Zl2 = popt[9];
	// 		Zr2 = popt[9];
	// 		// Residual coupling constant
	// 		Ar = popt[19];
	// 		//Residual mass
	// 		mr = popt[29];	
	// 		// result
	// 		result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
	// 		result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
	// 		result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
	// 		return result;
	// 	}


	 	return result;
	}

	number_type C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = Vl_Vl_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = Vl_Vs1_[t_index];
			return result;
		}
		if (i == 0 && j == 2)
		{
			result = Vl_Vs2_[t_index];
			return result;
		}
		if (i == 0 && j == 3)
		{
			result = Vl_Vs3_[t_index];
			return result;
		}
		if (i == 1 && j == 0)
		{
			result = Vs1_Vl_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = Vs1_Vs1_[t_index];
			return result;
		}
		if (i == 1 && j == 2)
		{
			result = Vs1_Vs2_[t_index];
			return result;
		}
		if (i == 1 && j == 3)
		{
			result = Vs1_Vs3_[t_index];
			return result;
		}
		if (i == 2 && j == 0)
		{
			result = Vs2_Vl_[t_index];
			return result;
		}
		if (i == 2 && j == 1)
		{
			result = Vs2_Vs1_[t_index];
			return result;
		}
		if (i == 2 && j == 2)
		{
			result = Vs2_Vs2_[t_index];
			return result;
		}
		if (i == 2 && j == 3)
		{
			result = Vs2_Vs3_[t_index];
			return result;
		}
		if (i == 3 && j == 0)
		{
			result = Vs3_Vl_[t_index];
			return result;
		}
		if (i == 3 && j == 1)
		{
			result = Vs3_Vs1_[t_index];
			return result;
		}
		if (i == 3 && j == 2)
		{
			result = Vs3_Vs2_[t_index];
			return result;
		}
		if (i == 3 && j == 3)
		{
			result = Vs3_Vs3_[t_index];
			return result;
		}
		
		

		return result;
	}

	number_type d_C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = d_Vl_Vl_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = d_Vl_Vs1_[t_index];
			return result;
		}
		if (i == 0 && j == 2)
		{
			result = d_Vl_Vs2_[t_index];
			return result;
		}
		if (i == 0 && j == 3)
		{
			result = d_Vl_Vs3_[t_index];
			return result;
		}
		if (i == 1 && j == 0)
		{
			result = d_Vs1_Vl_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = d_Vs1_Vs1_[t_index];
			return result;
		}
		if (i == 1 && j == 2)
		{
			result = d_Vs1_Vs2_[t_index];
			return result;
		}
		if (i == 1 && j == 3)
		{
			result = d_Vs1_Vs3_[t_index];
			return result;
		}
		if (i == 2 && j == 0)
		{
			result = d_Vs2_Vl_[t_index];
			return result;
		}
		if (i == 2 && j == 1)
		{
			result = d_Vs2_Vs1_[t_index];
			return result;
		}
		if (i == 2 && j == 2)
		{
			result = d_Vs2_Vs2_[t_index];
			return result;
		}
		if (i == 2 && j == 3)
		{
			result = d_Vs2_Vs3_[t_index];
			return result;
		}
		if (i == 3 && j == 0)
		{
			result = d_Vs3_Vl_[t_index];
			return result;
		}
		if (i == 3 && j == 1)
		{
			result = d_Vs3_Vs1_[t_index];
			return result;
		}
		if (i == 3 && j == 2)
		{
			result = d_Vs3_Vs2_[t_index];
			return result;
		}
		if (i == 3 && j == 3)
		{
			result = d_Vs3_Vs3_[t_index];
			return result;
		}
		
		

		return result;
	}

public:
	/* reduced chi2 sum as potential function to be minimized */
	number_type potential(const Vector<number_type> & q)
	{
		/* checks if data arrays have equal lengths */
		size_type N = t_.size();
		assert(Vl_Vl_.size() == N);
		assert(Vl_Vs1_.size() == N);
		assert(Vl_Vs2_.size() == N);
		assert(Vl_Vs3_.size() == N);
		assert(Vs1_Vl_.size() == N);
		assert(Vs1_Vs1_.size() == N);
		assert(Vs1_Vs2_.size() == N);
		assert(Vs1_Vs3_.size() == N);
		assert(Vs2_Vl_.size() == N);
		assert(Vs2_Vs1_.size() == N);
		assert(Vs2_Vs2_.size() == N);
		assert(Vs2_Vs3_.size() == N);
		assert(Vs3_Vl_.size() == N);
		assert(Vs3_Vs1_.size() == N);
		assert(Vs3_Vs2_.size() == N);
		assert(Vs3_Vs3_.size() == N);
		assert(d_Vl_Vl_.size() == N);
		assert(d_Vl_Vs1_.size() == N);
		assert(d_Vl_Vs2_.size() == N);
		assert(d_Vl_Vs3_.size() == N);
		assert(d_Vs1_Vl_.size() == N);
		assert(d_Vs1_Vs1_.size() == N);
		assert(d_Vs1_Vs2_.size() == N);
		assert(d_Vs1_Vs3_.size() == N);
		assert(d_Vs2_Vl_.size() == N);
		assert(d_Vs2_Vs1_.size() == N);
		assert(d_Vs2_Vs2_.size() == N);
		assert(d_Vs2_Vs3_.size() == N);
		assert(d_Vs3_Vl_.size() == N);
		assert(d_Vs3_Vs1_.size() == N);
		assert(d_Vs3_Vs2_.size() == N);
		assert(d_Vs3_Vs3_.size() == N);
		number_type chi2 = 0;
		for (size_type k = 0; k<t_.size(); ++k)
		{
			for (size_type i = 0; i < 2; ++i)
			{
				for (size_type j = 0; j < 2; ++j)
				{
					chi2 += pow(C(t_[k], i, j, q) - C_exp(k, i, j), 2)/pow(d_C_exp(k, i, j), 2);
				}
			}
		}

		return chi2/d_of_freedom();

	}

	/* Constraints to the parameters */
	bool constraints_respected(const Vector<number_type> &q)
	{
		bool respected = true;
		if (q[1] < q[0] )
		{
			respected = false; // mass of excited state smaller than ground state
			return respected;
		}

		/* masses are positive */
		if (q[0] < 0)
		{
			respected = false; 
			return respected;
		}

		if ((q[9] < q[1]) || (q[10] < q[1]) || (q[11] < q[1]) ) // residual masses have to be bigger than m1 and m2
		{
			respected = false; 
			return respected;
		}


		// restrictions to coupling constants that one can demand without loss of generality
		if (q[2] < 0 || q[3] < 0 )
		{
			respected = false;
			return respected;
		}

		// /* coupling constants */
		// if (q[2] * q[3] < 0) // diagonal coupling constants are positive
		// {
		// 	respected = false;
		// }

		return respected;
	}

private:
	// correlators
	Vector<number_type> t_;
	Vector<number_type> Vl_Vl_;
	Vector<number_type> Vl_Vs1_;
	Vector<number_type> Vl_Vs2_;
	Vector<number_type> Vl_Vs3_;
	Vector<number_type> Vs1_Vl_;
	Vector<number_type> Vs1_Vs1_;
	Vector<number_type> Vs1_Vs2_;
	Vector<number_type> Vs1_Vs3_;
	Vector<number_type> Vs2_Vl_;
	Vector<number_type> Vs2_Vs1_;
	Vector<number_type> Vs2_Vs2_;
	Vector<number_type> Vs2_Vs3_;
	Vector<number_type> Vs3_Vl_;
	Vector<number_type> Vs3_Vs1_;
	Vector<number_type> Vs3_Vs2_;
	Vector<number_type> Vs3_Vs3_;
	Vector<number_type> d_Vl_Vl_;
	Vector<number_type> d_Vl_Vs1_;
	Vector<number_type> d_Vl_Vs2_;
	Vector<number_type> d_Vl_Vs3_;
	Vector<number_type> d_Vs1_Vl_;
	Vector<number_type> d_Vs1_Vs1_;
	Vector<number_type> d_Vs1_Vs2_;
	Vector<number_type> d_Vs1_Vs3_;
	Vector<number_type> d_Vs2_Vl_;
	Vector<number_type> d_Vs2_Vs1_;
	Vector<number_type> d_Vs2_Vs2_;
	Vector<number_type> d_Vs2_Vs3_;
	Vector<number_type> d_Vs3_Vl_;
	Vector<number_type> d_Vs3_Vs1_;
	Vector<number_type> d_Vs3_Vs2_;
	Vector<number_type> d_Vs3_Vs3_;




};


#endif