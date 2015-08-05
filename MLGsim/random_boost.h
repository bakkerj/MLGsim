#ifndef RANDOM_BOOST_H
#define RANDOM_BOOST_H

#include <cmath>
#include <limits>
#include <boost/random.hpp>

namespace boost_rng
{
	boost::mt19937 rng_i;
	boost::lagged_fibonacci607 rng_f;
};

// generate random double out of [0,1)
double uniform()
{
	return boost_rng::rng_f();
}

void set_seed(unsigned seed)
{
	boost_rng::rng_i.seed(boost::mt19937::result_type(seed));
	boost_rng::rng_f.seed(boost_rng::rng_i);
}
	
double normal(double mean, double sig)
{
	static bool cached = false;
	static double cache;
	double x, y, r2;
	
	if (cached) {
		cached = false;
		return mean + sig * cache;
	}
		
	do {
		// choose x,y in uniform square (-1,-1) to (+1,+1)
		x = -1 + 2 * uniform();
		y = -1 + 2 * uniform();
		
		// see if it is in the unit circle
		r2 = x * x + y * y;
	} while (r2 > 1.0 || r2 == 0);
	
	const double factor = sqrt (-2.0 * log (r2) / r2);
	cache = x * factor;
	cached = true;
	
	// Box-Muller transform
	return mean + sig * y * factor;
}

// generate random int out of [0, maximum)
int random_number(int maximum)
{
	return boost_rng::rng_i() % maximum;  //should be ok, since MT generates good lower bits
}

#endif	//RANDOM_BOOST_H