// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef METROPOLIS_H
#define METROPOLIS_H
#include "points.h"
#include "pmt.h"
#include "metro_stats.h"
#include "mcmc_params.h"

void SO_sweep(points_t* ppoints, mcmc_params_t* pparams,
	double* alphas, double* pH, double* pD, double* pV,
	metro_stats_t* pmetro_stats);

void R_sweep(points_t* ppoints, mcmc_params_t* pparams);

void worm_sweep(points_t* ppoints, mcmc_params_t* pparams,
	int do_print_worm, int do_print_wormspan,
	double* alphas, double* pH, double* pD, double* pV,
	metro_stats_t* pmetro_stats);

#endif // METROPOLIS_H
