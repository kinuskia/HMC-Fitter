#ifndef CURVE_FIT_HPP
#define CURVE_FIT_HPP


template <typename REAL>
class Curve_fit
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	/* Default constructor */
	Curve_fit(
	HMC<number_type> sampler
	)
	: sampler_(sampler)
	{}

	size_type d_of_freedom() const
	{
		return sampler_.model_.t_.size() - sampler_.model_.n_parameters();
	}


private:
	// HMC sampler which also includes fitting model
	HMC<number_type> sampler_;

	
};

#endif