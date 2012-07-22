// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include "rho.h"
#include "util.h"
#include "vector_misc.h"

// ----------------------------------------------------------------
// rho_k(pi) = 1/|X| #{x in X: ell_x(pi) \le k}
//
// Output:  rho must be an array of length N+1.
// Scratch: counts must be an array of length N+1.
//
// Example:
// Cycle decomposition is (1) (2) (3) (4) (5 6) (7 8) (9 10 11)
// Cycle type is [3 2 2 1 1 1 1]
// Cycle counts are [0 4 2 1 0 0 0 0 0 0 0 0]
//   meaning that there are 4 1-cycles, 2 2-cycles, and 1 3-cycle,
//   and no cycles of any other length.

void get_rho_L_pi(points_t* ppoints, double* rho, int* counts, int *plmax)
{
	int N = ppoints->N; // N = |X|
	double recip_N = 1.0 / N;
	int j;
	double rhoval;

	pmt_get_cycle_counts_and_lmax(ppoints, counts, plmax);

	rhoval = 0.0;
	for (j = 0; j <= N; j++) {
		rhoval += j * counts[j];
		rho[j] = rhoval * recip_N;
	}
}

// ----------------------------------------------------------------
// Horizontal axis runs from 0 to N.
// Vertical   axis runs from 0 to 1.
//
// +------------------------+
// |         *              |
// |      *                 |
// |    *                   |
// |  *                     |
// |                        |
// |*                       |
// |                        |
// |                        |
// |                        |
// |                        |
// +------------------------+

double compute_rho_infty(double* rho_fit, double* rho_avg, double* rho,
	int N, int verbose)
{
	int k;
	int maxk = 0;
	double dev;
	double maxdev = 0.0;
	double rho_finite_T;
	double rho_infty_T;

	for (k = 1; k <= N; k++) {
		dev = rho_avg[k] - (double)k/(double)N;
		if (dev > maxdev) {
			maxdev = dev;
			maxk = k;
		}
	}

	rho_finite_T = rho_avg[maxk] - (double)maxk/(double)N;
	rho_infty_T  = 1.0 - rho_finite_T;

	printf("# rho_infty max dev = %11.7lf at k = %d\n", maxdev, maxk);

	for (k = 0; k <= N; k++) {
		rho_fit[k] = (double)k/(double)N + rho_finite_T;
		if (rho_fit[k] > 1.0)
			rho_fit[k] = 1.0;
	}

	return rho_infty_T;
}
