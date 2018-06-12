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
	//0		Pl
	//1		Ps
	//2		A0l
	//3		A0s
	0		A0l
	1		A0s
	


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
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
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
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
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
		// 	Ar = popt[10];
		// 	//Residual mass
		// 	mr = popt[16];
		// 	// result
		// 	result += Zl1*Zr1/2./m1 * (exp(-m1*t) - exp(-m1*(T-t)));
		// 	result += Zl2*Zr2/2./m2 * (exp(-m2*t) - exp(-m2*(T-t)));
		// 	result += Ar * (exp(-mr*t) - exp(-mr*(T-t)));
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
			result += Zl1*Zr1/m1 * exp(-m1*T/2.) * sinh(m1*(T/2.-t));
			result += Zl2*Zr2/m2 * exp(-m2*T/2.) * sinh(m2*(T/2.-t));
			result += 2.*Ar * exp(-mr*T/2.) * sinh(mr*(T/2.-t));
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
		// 	Ar = popt[12];
		// 	//Residual mass
		// 	mr = popt[18];
		// 	// result
		// 	result += Zl1*Zr1/2./m1 * (exp(-m1*t) - exp(-m1*(T-t)));
		// 	result += Zl2*Zr2/2./m2 * (exp(-m2*t) - exp(-m2*(T-t)));
		// 	result += Ar * (exp(-mr*t) - exp(-mr*(T-t)));
		// 	return result;	
		// }

		// if ((i == 2 && j == 2))
		// {
		// 	// Coupling constants of mass 1
		// 	Zl1 = popt[6];
		// 	Zr1 = popt[6];
		// 	// Coupling constants of mass 2
		// 	Zl2 = popt[7];
		// 	Zr2 = popt[7];
		// 	// Residual coupling constant
		// 	Ar = popt[13];
		// 	//Residual mass
		// 	mr = popt[19];
		// 	// result
		// 	result += Zl1*Zr1/2./m1 * (exp(-m1*t) + exp(-m1*(T-t)));
		// 	result += Zl2*Zr2/2./m2 * (exp(-m2*t) + exp(-m2*(T-t)));
		// 	result += Ar * (exp(-mr*t) + exp(-mr*(T-t)));
		// 	return result;	
		// }

		return result;
	}

	number_type C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = A0l_A0l_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = A0l_A0s_[t_index];
			return result;
		}
		// if (i == 0 && j == 2)
		// {
		// 	result = Pl_A0l_[t_index];
		// 	return result;
		// }
		if (i == 1 && j == 0)
		{
			result = A0s_A0l_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = A0s_A0s_[t_index];
			return result;
		}
		// if (i == 1 && j == 2)
		// {
		// 	result = Ps_A0l_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 0)
		// {
		// 	result = A0l_Pl_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 1)
		// {
		// 	result = A0l_Ps_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 2)
		// {
		// 	result = A0l_A0l_[t_index];
		// 	return result;
		// }

		return result;
	}

	number_type d_C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result = 0;
		if (i == 0 && j == 0)
		{
			result = d_A0l_A0l_[t_index];
			return result;
		}
		if (i == 0 && j == 1)
		{
			result = d_A0l_A0s_[t_index];
			return result;
		}
		// if (i == 0 && j == 2)
		// {
		// 	result = d_Pl_A0l_[t_index];
		// 	return result;
		// }
		if (i == 1 && j == 0)
		{
			result = d_A0s_A0l_[t_index];
			return result;
		}
		if (i == 1 && j == 1)
		{
			result = d_A0s_A0s_[t_index];
			return result;
		}
		// if (i == 1 && j == 2)
		// {
		// 	result = d_Ps_A0l_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 0)
		// {
		// 	result = d_A0l_Pl_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 1)
		// {
		// 	result = d_A0l_Ps_[t_index];
		// 	return result;
		// }
		// if (i == 2 && j == 2)
		// {
		// 	result = d_A0l_A0l_[t_index];
		// 	return result;
		// }


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
		if (q[2] < 0 || q[4] < 0)
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