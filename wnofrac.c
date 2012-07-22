// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rcmrand.h"
#include "pmt.h"

int has_only_short_jumps(points_t* ppoints);

// ----------------------------------------------------------------
int main(int argc, char **argv)
{
	int do_seed = 0;
	int seed    = 0;

	points_t* ppoints = 0;
	int L = 40;
	int d = 3;
	int Wx, Wy, Wz;
	int W2;

	//int reps = 1000000;
	int reps = 1000;
	int repi;

	int num_only_short_jumps = 0;

	// Seed the pseudorandom-number generator.
	if (do_seed)
		SRANDOM(seed);
	else
		STRANDOM();

	// Lattice points and permutations:
	ppoints = get_cubic_lattice_points(L, d);
	set_unif_rand_pmt(ppoints);
	set_up_cycinfo_list(ppoints);

	for (repi = 0; repi < reps; repi++) {
		set_unif_rand_pmt(ppoints);
		recompute_cycinfo_list(ppoints);

		// Skip if not short jump lengths.
		//
		if (!has_only_short_jumps(ppoints))
			continue;
		num_only_short_jumps++;

		// Compute Wx, Wy, Wz
		get_pmt_winding_numbers(ppoints, &Wx, &Wy, &Wz);
		printf("%5d %5d %5d\n", Wx, Wy, Wz);
		W2 = Wx*Wx + Wy*Wy + Wz*Wz;

		// Compute W2
	}

	return 0;
}

// ----------------------------------------------------------------
int has_only_short_jumps(points_t* ppoints)
{
	int i, j, k;
	int L = ppoints->dims.nx;
	point_t*** level1;
	point_t**  level2;
	point_t*   level3;
	point_t*   px;

	level1 = ppoints->lattice;
	for (i = 0; i < L; i++) {
		level2 = level1[i];
		for (j = 0; j < L; j++) {
			level3 = level2[j];
			for (k = 0; k < L; k++) {
				px = &level3[k];
				if (px->fwd_d >= 5) {
					printf("%lf\n", px->fwd_d);
					return 0;
				}
			}
		}
	}
	return 1;
}
