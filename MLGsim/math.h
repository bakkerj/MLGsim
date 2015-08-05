#pragma once
#ifndef __INC_MGL_MATH_HH__
#define __INC_MGL_MATH_HH__

#include "random_boost.h"

namespace mlg {

  class Math {
  public:

// changed to replace unix-specific function random(): uniform() is based on Boost libraries
		//static double fRandom() { return (double)random() / (RAND_MAX + 1.0); }
    static double fRandom() { return uniform(); }	

    static double binom(size_t n, size_t k)
    {
			//std::cout << "n = " << n << " k = " << k << std::endl;
      assert(n >= k);
      double r = 1.0;
      for (; k; --k, --n)
				r *= (double)n / (double)k;
      return r;
    }

// added to replace unix-specific function nearbyint()
		static double round(double num)
		{
			if ((num + 0.5) >= (int(num) + 1))
				return double(int(num) + 1);
			else
				return double(int(num));
		}
  };
};

#endif // __INC_MGL_MATH_HH__
