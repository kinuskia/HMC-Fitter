#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include <assert.h>

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
	, Pl_A0l_(A0s_A0l)
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
	, d_Pl_A0l_(d_A0s_A0l)
	, d_Pl_A0s_(d_Pl_A0s)
	, d_Pl_Pl_(d_Pl_Pl)
	, d_Pl_Ps_(d_Pl_Ps)
	, d_Ps_A0l_(d_Ps_A0l)
	, d_Ps_A0s_(d_Ps_A0s)
	, d_Ps_Pl_(d_Ps_Pl)
	, d_Ps_Ps_(d_Ps_Ps)
	{}

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 6;
	}
public: // needs to become private once I focus on intrinsic errors
	/* degrees of freedom */
	size_type d_of_freedom() const
	{
		return t_.size() * 1 * 1 - n_parameters();
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
	4		A_00
	5		m_00

	*/
	number_type C(number_type t, size_type i, size_type j, const Vector<number_type> & popt)
	{
		assert(n_parameters() == popt.size());
		number_type T = 64.;
		number_type m1 = popt[0];
		number_type m2 = popt[1];
		number_type Z11;
		number_type Z12;
		number_type Z21;
		number_type Z22;
		number_type Ar;
		number_type mr;
		number_type result;
		if (i == 0 && j == 0)
		{
			Z11 = popt[2];
			Z21 = popt[2];
			Z12 = popt[3];
			Z22 = popt[3];
			Ar = popt[4];
			mr = popt[5];
			result += Z11*Z21/2./m1 * (exp(-m1*t) + exp(-m1*(T-t)));
			result += Z12*Z22/2./m2 * (exp(-m2*t) + exp(-m2*(T-t)));
			result += Ar * (exp(-mr*t) + exp(-mr*(T-t)));
		}

		return result;
	}
	number_type C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result;
		if (i == 0 && j == 0)
		{
			result = Pl_Pl_[t_index];
		}
		return result;
	}

	number_type d_C_exp(size_type t_index, size_type i, size_type j)
	{
		number_type result;
		if (i == 0 && j == 0)
		{
			result = d_Pl_Pl_[t_index];
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
			for (size_type i = 0; i < 1; ++i)
			{
				for (size_type j = 0; j < 1; ++j)
				{
					chi2 += pow(C(t_[k], i, j, q) - C_exp(k, i, j), 2)/pow(d_C_exp(k, i, j), 2);
				}
			}
		}

		return chi2/d_of_freedom();

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