// ================================================================
// Copyright (c) 2004 John Kerl.
// kerl.john.r@gmail.com
//
// This code and information is provided as is without warranty of
// any kind, either expressed or implied, including but not limited to
// the implied warranties of merchantability and/or fitness for a
// particular purpose.
//
// No restrictions are placed on copy or reuse of this code, as long
// as these paragraphs are included in the code.
// ================================================================

#include <math.h>
#include "rcmrand.h"

// Box-Muller, from _Numerical Recipes in C_.
// Returns a normally distributed deviate with zero mean and unit variance.

// ----------------------------------------------------------------
double gasdev(void)
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if (iset == 0) {
		// We don't have an extra deviate handy, so pick two uniform
		// numbers in the unit square, repeating until we get a pair
		// in the unit circle.
		do {
			v1 = 2.0 * URANDOM() - 1.0;
			v2 = 2.0 * URANDOM() - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		// Now make the Box-Muller transformation to get two normal
		// deviates.  Return one and save the other for next time.
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else {
		iset = 0;
		return gset;
	}
}

// ----------------------------------------------------------------
double gasdev_ms(double mean, double stddev)
{
	return mean + stddev * gasdev();
}
