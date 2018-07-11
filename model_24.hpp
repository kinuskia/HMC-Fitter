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
	Vector<number_type> A0l_A0l,
	Vector<number_type> A0l_A0s,
	Vector<number_type> A0l_Pl,
	Vector<number_type> A0l_Ps,
	Vector<number_type> A0s_A0l,
	Vector<number_type> A0s_A0s,
	Vector<number_type> A0s_Pl,
	Vector<number_type> A0s_Ps,
	Vector<number_type> Pl_A0l,
	Vector<number_type> Pl_A0s,
	Vector<number_type> Pl_Pl,
	Vector<number_type> Pl_Ps,
	Vector<number_type> Ps_A0l,
	Vector<number_type> Ps_A0s,
	Vector<number_type> Ps_Pl,
	Vector<number_type> Ps_Ps,
	Vector<number_type> d_A0l_A0l,
	Vector<number_type> d_A0l_A0s,
	Vector<number_type> d_A0l_Pl,
	Vector<number_type> d_A0l_Ps,
	Vector<number_type> d_A0s_A0l,
	Vector<number_type> d_A0s_A0s,
	Vector<number_type> d_A0s_Pl,
	Vector<number_type> d_A0s_Ps,
	Vector<number_type> d_Pl_A0l,
	Vector<number_type> d_Pl_A0s,
	Vector<number_type> d_Pl_Pl,
	Vector<number_type> d_Pl_Ps,
	Vector<number_type> d_Ps_A0l,
	Vector<number_type> d_Ps_A0s,
	Vector<number_type> d_Ps_Pl,
	Vector<number_type> d_Ps_Ps
	) // default constructor
	: t_(t)
	, A0l_A0l_(A0l_A0l)
	, A0l_A0s_(A0l_A0s)
	, A0l_Pl_(A0l_Pl)
	, A0l_Ps_(A0l_Ps)
	, A0s_A0l_(A0s_A0l)
	, A0s_A0s_(A0s_A0s)
	, A0s_Pl_(A0s_Pl)
	, A0s_Ps_(A0s_Ps)
	, Pl_A0l_(Pl_A0l)
	, Pl_A0s_(Pl_A0s)
	, Pl_Pl_(Pl_Pl)
	, Pl_Ps_(Pl_Ps)
	, Ps_A0l_(Ps_A0l)
	, Ps_A0s_(Ps_A0s)
	, Ps_Pl_(Ps_Pl)
	, Ps_Ps_(Ps_Ps)
	, d_A0l_A0l_(d_A0l_A0l)
	, d_A0l_A0s_(d_A0l_A0s)
	, d_A0l_Pl_(d_A0l_Pl)
	, d_A0l_Ps_(d_A0l_Ps)
	, d_A0s_A0l_(d_A0s_A0l)
	, d_A0s_A0s_(d_A0s_A0s)
	, d_A0s_Pl_(d_A0s_Pl)
	, d_A0s_Ps_(d_A0s_Ps)
	, d_Pl_A0l_(d_Pl_A0l)
	, d_Pl_A0s_(d_Pl_A0s)
	, d_Pl_Pl_(d_Pl_Pl)
	, d_Pl_Ps_(d_Pl_Ps)
	, d_Ps_A0l_(d_Ps_A0l)
	, d_Ps_A0s_(d_Ps_A0s)
	, d_Ps_Pl_(d_Ps_Pl)
	, d_Ps_Ps_(d_Ps_Ps)
	{}

	// /* print content of three vectors to console */
	// void print_content() const
	// {
	// 	Vector<number_type> x = Ps_Ps_;
	// 	Vector<number_type> dx = d_Ps_Ps_;
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
		return 24;
	}
public: // needs to become private once I focus on intrinsic errors
	/* degrees of freedom */
	size_type d_of_freedom() const
	{
		return t_.size() * 4 * 4 - n_parameters();
	}
private:
	/* fitting model for correlator C_ij(t) */
	/*
	nomenclature for correlators
	index	correlator
	0		Pl
	1		Ps
	2		A0l
	3		A0s
	
	
	nomenclature for fitting parameters
	index	parameter
	0		mass ground state
	1 		mass of first excited state
	2		Z0 of m1
	3		Z0 of m2
	4		Z1 of m1
	5		Z1 of m2
	6		Z2 of m1
	7		Z2 of m2
	8		Z3 of m1
	9		Z3 of m2
	10		A_00
	11		A_01
	12		A_02
	13		A_03
	14		A_12
	15		A_22
	16		A_23
	17		m_00
	18		m_01
	19		m_02
	20		m_03
	21		m_12
	22		m_22
	23		m_23
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
			Ar = popt[10];
			//Residual mass
			mr = popt[17];
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
			Ar = popt[11];
			//Residual mass
			mr = popt[18];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;	
		}

		if ((i == 0 && j == 2) || (i == 2 && j == 0))
		{
			// Coupling constants of mass 1
			Zl1 = popt[2];
			Zr1 = popt[6];
			// Coupling constants of mass 2
			Zl2 = popt[3];
			Zr2 = popt[7];
			// Residual coupling constant
			Ar = popt[12];
			//Residual mass
			mr = popt[19];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
			return result;	
		}

		if ((i == 0 && j == 3) || (i == 3 && j == 0))
		{
			// Coupling constants of mass 1
			Zl1 = popt[2];
			Zr1 = popt[8];
			// Coupling constants of mass 2
			Zl2 = popt[3];
			Zr2 = popt[9];
			// Residual coupling constant
			Ar = popt[13];
			//Residual mass
			mr = popt[20];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
			return result;	
		}

		
		if (i == 1 && j == 1)
		{
			// Coupling constants of mass 1
			Zl1 = popt[4];
			Zr1 = popt[4];
			// Coupling constants of mass 2
			Zl2 = popt[5];
			Zr2 = popt[5];
			// Residual coupling constant
			//Ar = popt[14];
			//Residual mass
			//mr = popt[22];	
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			//result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;
		}

		if ((i == 1 && j == 2) || (i == 2 && j == 1))
		{
			// Coupling constants of mass 1
			Zl1 = popt[4];
			Zr1 = popt[6];
			// Coupling constants of mass 2
			Zl2 = popt[5];
			Zr2 = popt[7];
			// Residual coupling constant
			Ar = popt[14];
			//Residual mass
			mr = popt[21];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
			return result;	
		}

		if ((i == 1 && j == 3) || (i == 3 && j == 1))
		{
			// Coupling constants of mass 1
			Zl1 = popt[4];
			Zr1 = popt[8];
			// Coupling constants of mass 2
			Zl2 = popt[5];
			Zr2 = popt[9];
			// Residual coupling constant
			//Ar = popt[16];
			//Residual mass
			//mr = popt[26];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			//result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
			return result;	
		}

		if (i == 2 && j == 2)
		{
			// Coupling constants of mass 1
			Zl1 = popt[6];
			Zr1 = popt[6];
			// Coupling constants of mass 2
			Zl2 = popt[7];
			Zr2 = popt[7];
			// Residual coupling constant
			Ar = popt[15];
			//Residual mass
			mr = popt[22];	
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;
		}

		if ((i == 2 && j == 3) || (i == 3 && j == 2))
		{
			// Coupling constants of mass 1
			Zl1 = popt[6];
			Zr1 = popt[8];
			// Coupling constants of mass 2
			Zl2 = popt[7];
			Zr2 = popt[9];
			// Residual coupling constant
			Ar = popt[16];
			//Residual mass
			mr = popt[23];
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;	
		}

		if (i == 3 && j == 3)
		{
			// Coupling constants of mass 1
			Zl1 = popt[8];
			Zr1 = popt[8];
			// Coupling constants of mass 2
			Zl2 = popt[9];
			Zr2 = popt[9];
			// Residual coupling constant
			//Ar = popt[19];
			//Residual mass
			//mr = popt[29];	
			// result
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * cosh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * cosh(m2*(T/2.-t));
			//result += 2.*Ar * exp(-mr*T/2.) * cosh(mr*(T/2.-t));
			return result;
		}


		return result;
	}

	number_type C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = Pl_Pl_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = Pl_Ps_[t_index];
			return result;
		}
		if (i == 0 && j == 2)
		{
			result = Pl_A0l_[t_index];
			return result;
		}
		if (i == 0 && j == 3)
		{
			result = Pl_A0s_[t_index];
			return result;
		}
		if (i == 1 && j == 0)
		{
			result = Ps_Pl_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = Ps_Ps_[t_index];
			return result;
		}
		if (i == 1 && j == 2)
		{
			result = Ps_A0l_[t_index];
			return result;
		}
		if (i == 1 && j == 3)
		{
			result = Ps_A0s_[t_index];
			return result;
		}
		if (i == 2 && j == 0)
		{
			result = A0l_Pl_[t_index];
			return result;
		}
		if (i == 2 && j == 1)
		{
			result = A0l_Ps_[t_index];
			return result;
		}
		if (i == 2 && j == 2)
		{
			result = A0l_A0l_[t_index];
			return result;
		}
		if (i == 2 && j == 3)
		{
			result = A0l_A0s_[t_index];
			return result;
		}
		if (i == 3 && j == 0)
		{
			result = A0s_Pl_[t_index];
			return result;
		}
		if (i == 3 && j == 1)
		{
			result = A0s_Ps_[t_index];
			return result;
		}
		if (i == 3 && j == 2)
		{
			result = A0s_A0l_[t_index];
			return result;
		}
		if (i == 3 && j == 3)
		{
			result = A0s_A0s_[t_index];
			return result;
		}
		
		

		return result;
	}

	number_type d_C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = d_Pl_Pl_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = d_Pl_Ps_[t_index];
			return result;
		}
		if (i == 0 && j == 2)
		{
			result = d_Pl_A0l_[t_index];
			return result;
		}
		if (i == 0 && j == 3)
		{
			result = d_Pl_A0s_[t_index];
			return result;
		}
		if (i == 1 && j == 0)
		{
			result = d_Ps_Pl_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = d_Ps_Ps_[t_index];
			return result;
		}
		if (i == 1 && j == 2)
		{
			result = d_Ps_A0l_[t_index];
			return result;
		}
		if (i == 1 && j == 3)
		{
			result = d_Ps_A0s_[t_index];
			return result;
		}
		if (i == 2 && j == 0)
		{
			result = d_A0l_Pl_[t_index];
			return result;
		}
		if (i == 2 && j == 1)
		{
			result = d_A0l_Ps_[t_index];
			return result;
		}
		if (i == 2 && j == 2)
		{
			result = d_A0l_A0l_[t_index];
			return result;
		}
		if (i == 2 && j == 3)
		{
			result = d_A0l_A0s_[t_index];
			return result;
		}
		if (i == 3 && j == 0)
		{
			result = d_A0s_Pl_[t_index];
			return result;
		}
		if (i == 3 && j == 1)
		{
			result = d_A0s_Ps_[t_index];
			return result;
		}
		if (i == 3 && j == 2)
		{
			result = d_A0s_A0l_[t_index];
			return result;
		}
		if (i == 3 && j == 3)
		{
			result = d_A0s_A0s_[t_index];
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
		assert(A0l_A0l_.size() == N);
		assert(A0l_A0s_.size() == N);
		assert(A0l_Pl_.size() == N);
		assert(A0l_Ps_.size() == N);
		assert(A0s_A0l_.size() == N);
		assert(A0s_A0s_.size() == N);
		assert(A0s_Pl_.size() == N);
		assert(A0s_Ps_.size() == N);
		assert(Pl_A0l_.size() == N);
		assert(Pl_A0s_.size() == N);
		assert(Pl_Pl_.size() == N);
		assert(Pl_Ps_.size() == N);
		assert(Ps_A0l_.size() == N);
		assert(Ps_A0s_.size() == N);
		assert(Ps_Pl_.size() == N);
		assert(Ps_Ps_.size() == N);
		assert(d_A0l_A0l_.size() == N);
		assert(d_A0l_A0s_.size() == N);
		assert(d_A0l_Pl_.size() == N);
		assert(d_A0l_Ps_.size() == N);
		assert(d_A0s_A0l_.size() == N);
		assert(d_A0s_A0s_.size() == N);
		assert(d_A0s_Pl_.size() == N);
		assert(d_A0s_Ps_.size() == N);
		assert(d_Pl_A0l_.size() == N);
		assert(d_Pl_A0s_.size() == N);
		assert(d_Pl_Pl_.size() == N);
		assert(d_Pl_Ps_.size() == N);
		assert(d_Ps_A0l_.size() == N);
		assert(d_Ps_A0s_.size() == N);
		assert(d_Ps_Pl_.size() == N);
		assert(d_Ps_Ps_.size() == N);
		number_type chi2 = 0;
		for (size_type k = 0; k<t_.size(); ++k)
		{
			for (size_type i = 0; i < 4; ++i)
			{
				for (size_type j = 0; j < 4; ++j)
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

		if ((q[17] < q[1]) || (q[18] < q[1]) || (q[19] < q[1]) || (q[20] < q[1]) || (q[21] < q[1]) || (q[22] < q[1]) || (q[23] < q[1])) // residual masses have to be bigger than m1 and m2
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
	Vector<number_type> A0l_A0l_;
	Vector<number_type> A0l_A0s_;
	Vector<number_type> A0l_Pl_;
	Vector<number_type> A0l_Ps_;
	Vector<number_type> A0s_A0l_;
	Vector<number_type> A0s_A0s_;
	Vector<number_type> A0s_Pl_;
	Vector<number_type> A0s_Ps_;
	Vector<number_type> Pl_A0l_;
	Vector<number_type> Pl_A0s_;
	Vector<number_type> Pl_Pl_;
	Vector<number_type> Pl_Ps_;
	Vector<number_type> Ps_A0l_;
	Vector<number_type> Ps_A0s_;
	Vector<number_type> Ps_Pl_;
	Vector<number_type> Ps_Ps_;
	Vector<number_type> d_A0l_A0l_;
	Vector<number_type> d_A0l_A0s_;
	Vector<number_type> d_A0l_Pl_;
	Vector<number_type> d_A0l_Ps_;
	Vector<number_type> d_A0s_A0l_;
	Vector<number_type> d_A0s_A0s_;
	Vector<number_type> d_A0s_Pl_;
	Vector<number_type> d_A0s_Ps_;
	Vector<number_type> d_Pl_A0l_;
	Vector<number_type> d_Pl_A0s_;
	Vector<number_type> d_Pl_Pl_;
	Vector<number_type> d_Pl_Ps_;
	Vector<number_type> d_Ps_A0l_;
	Vector<number_type> d_Ps_A0s_;
	Vector<number_type> d_Ps_Pl_;
	Vector<number_type> d_Ps_Ps_;




};


#endif