// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef RHO_H
#define RHO_H
#include "points.h"
#include "pmt.h"

// ----------------------------------------------------------------
// Output:  rho must be an array of length N+1.
// plmax is a reference argument which contains pi's lmax upon return.
// Scratch: counts must be an array of length N+1.
void get_rho_L_pi(points_t* ppoints, double* rho, int* counts, int *plmax);

// See the paper for details.  rho_avg and rho_fit are 0:N arrays (of
// length N+1).
double compute_rho_infty(double* rho_fit, double* rho_avg, double* rho,
	int N, int verbose);

#endif // RHO_H
