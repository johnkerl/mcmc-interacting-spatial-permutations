// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef DOTPLOT_H
#define DOTPLOT_H

#include "points.h"
#include "pmt.h"

void print_points_and_cycles(points_t* ppoints);
void print_points_and_cycles_min_max(points_t* ppoints,
	int min_ell, int max_ell);

#endif // DOTPLOT_H
