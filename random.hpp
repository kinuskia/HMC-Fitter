#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>


template<typename REAL>
class Uniform_real_distribution
{
public:
	typedef REAL number_type;
	typedef std::size_t size_type;
	Uniform_real_distribution( number_type min, number_type max)
	: rd_()
	, gen_(rd_())
	, min_(min)
	, max_(max)
	, dis_unif_(min, max)
	{}
	Uniform_real_distribution (const Uniform_real_distribution<REAL> & other)  // copy constructor
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->min_ = other.min_;
		this->max_ = other.max_;
		dis_unif_ = std::uniform_real_distribution<>(other.min_, other.max_);
	};
	Uniform_real_distribution & operator= (const Uniform_real_distribution<REAL> & other)  // copy assignment operator
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->min_ = other.min_;
		this->max_ = other.max_;
		dis_unif_ = std::uniform_real_distribution<>(other.min_, other.max_);
		return *this;
	};

	
	number_type draw()
	{
		number_type result = dis_unif_(gen_);
		return result;
	}


private:
	std::random_device rd_;
	std::mt19937 gen_;
	number_type min_;
	number_type max_;
	std::uniform_real_distribution<> dis_unif_;
};


template<typename REAL>
class Uniform_int_distribution
{
public:
	typedef REAL number_type;
	typedef std::size_t size_type;
	Uniform_int_distribution( number_type min, number_type max)
	: rd_()
	, gen_(rd_())
	, min_(min)
	, max_(max)
	, dis_unif_(min, max)
	{}
	Uniform_int_distribution (const Uniform_int_distribution<REAL> & other)  // copy constructor
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->min_ = other.min_;
		this->max_ = other.max_;
		dis_unif_ = std::uniform_int_distribution<>(other.min_, other.max_);
	};
	Uniform_int_distribution & operator= (const Uniform_int_distribution<REAL> & other)  // copy assignment operator
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->min_ = other.min_;
		this->max_ = other.max_;
		dis_unif_ = std::uniform_int_distribution<>(other.min_, other.max_);
		return *this;
	};

	
	number_type draw()
	{
		number_type result = dis_unif_(gen_);
		return result;
	}


private:
	std::random_device rd_;
	std::mt19937 gen_;
	number_type min_;
	number_type max_;
	std::uniform_int_distribution<> dis_unif_;
};


template<typename REAL>
class Normal_distribution
{
public:
	typedef REAL number_type;
	typedef std::size_t size_type;
	Normal_distribution( number_type mu, number_type sigma)
	: rd_()
	, gen_(rd_())
	, mu_(mu)
	, sigma_(sigma)
	, dis_norm_(mu, sigma)
	{}
	Normal_distribution (const Normal_distribution<REAL> & other)  // copy constructor
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->mu_ = other.mu_;
		this->sigma_ = other.sigma_;
		dis_norm_ = std::normal_distribution<>(other.mu_, other.sigma_);
	};
	Normal_distribution & operator= (const Normal_distribution<REAL> & other)  // copy assignment operator
	{
		(this->rd_)();
		(this->gen_)=std::mt19937(rd_());
		this->mu_ = other.mu_;
		this->sigma_ = other.sigma_;
		dis_norm_ = std::normal_distribution<>(other.mu_, other.sigma_);
		return *this;
	};

	
	number_type draw()
	{
		number_type result = dis_norm_(gen_);
		return result;
	}


private:
	std::random_device rd_;
	std::mt19937 gen_;
	number_type mu_;
	number_type sigma_;
	std::normal_distribution<> dis_norm_;
};



#endif