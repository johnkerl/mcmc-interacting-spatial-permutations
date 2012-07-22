// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcmc_params.h"
#include "util.h"
#include "rcmrand.h"
#include "points.h"
#include "checks.h"

static void check_magic_number(unsigned file_value,
	char* var_name, char* func_name);

// ----------------------------------------------------------------
void print_mcmc_params(mcmc_params_t* p)
{
	printf("# RNG = %s\n", get_rng_type_desc(p->rng_type));
	if (p->do_seed)
		printf("# RNG seed = 0x%08x\n", p->seed);
	printf("# L = %d d = %d N = %d\n", p->L, p->d, ipow(p->L, p->d));
#ifdef PERIODIC_BOUNDARY_CONDITIONS
	printf("# Boundary conditions: periodic.\n");
#else
	printf("# Boundary conditions: Dirichlet.\n");
#endif
	printf("# T = %11.7lf beta = %11.7lf alpha0 = %11.7lf\n",
		1.0 / p->beta, p->beta, p->alpha0);

	printf("# gamma = %11.7lf\n", p->gamma);
	printf("# Interactions:  %s.\n",
		get_interaction_type_desc(p->interaction_type));
	printf("# Initial permutation: %s.\n",
		p->rand_init_pi ? "uniform random" :
		p->wx1_init_pi  ? "Wx=1" :
		p->wy1_init_pi  ? "Wy=1" :
		p->wz1_init_pi  ? "Wz=1" :
		"identity");
	printf("# Site selection for Metropolis sweeps: %s.\n",
		p->do_random_sweep ? "random" : "sequential");

#ifdef CHECK_H
	printf("# H checking enabled.\n");
#endif
#ifdef CHECK_CYCINFO
	printf("# Cycinfo checking enabled.\n");
#endif

	if (p->force_therm == -1) {
		printf("# Thermalization is detected by %d turning points of system "
			"energy\n", p->therm_ntp);
		printf("# smoothed over a %d-point window.\n",
			p->smoothing_window_size);
	}
	else {
		printf("# %d pre-accumulation sweeps having occurred.\n",
			p->force_therm);
	}

	printf("# Terminate after %d accumulation%s:\n",
		p->termination_num_acc,
		p->termination_num_acc == 1 ? "" : "s");
	printf("# * %d %s  sweep%s per accumulation.\n",
		p->gsweeps_per_acc,
		p->allow_cycle_reversals ? "SAR" : "SO",
		p->gsweeps_per_acc == 1 ? "" : "s");

	printf("# * %d worm sweep%s per accumulation.\n",
		p->wsweeps_per_acc, p->wsweeps_per_acc == 1 ? "" : "s");
}

// ----------------------------------------------------------------
void set_default_mcmc_params(mcmc_params_t* p)
{
	memset((void*)p, 0xee, sizeof(*p));

	p->magic_number            = RCM_MAGIC_NUMBER;

	p->rng_type                = RCM_WHICH;
	p->seed                    = DEFAULT_SEED;
	p->do_seed                 = DEFAULT_DO_SEED;

	p->L                       = DEFAULT_L;
	p->d                       = DEFAULT_D;

	p->beta                    = DEFAULT_BETA;
	p->alpha0                  = DEFAULT_ALPHA0;
	p->interaction_type        = DEFAULT_INTERACTION_TYPE;
	p->gamma                   = DEFAULT_GAMMA;

	p->rand_init_pi            = DEFAULT_RAND_INIT_PI;
	p->wx1_init_pi             = DEFAULT_WX1_INIT_PI;
	p->wy1_init_pi             = DEFAULT_WY1_INIT_PI;
	p->wz1_init_pi             = DEFAULT_WZ1_INIT_PI;
	p->do_random_sweep         = DEFAULT_DO_RANDOM_SWEEP;
	p->allow_cycle_reversals   = DEFAULT_ALLOW_CYCLE_REVERSALS;
	p->gsweeps_per_acc         = DEFAULT_GSWEEPS_PER_ACC;
	p->wsweeps_per_acc         = DEFAULT_WSWEEPS_PER_ACC;
	p->termination_num_acc     = DEFAULT_TERMINATION_NUM_ACC;

	p->force_therm             = DEFAULT_FORCE_THERM;
	p->smoothing_window_size   = DEFAULT_SMOOTHING_WINDOW_SIZE;
	p->therm_ntp               = DEFAULT_THERM_NTP;

	p->magic_number2           = RCM_MAGIC_NUMBER;
}

// ----------------------------------------------------------------
void get_mcmc_params_from_file(mcmc_params_t* pparams, FILE* fp)
{
	read_from_file_or_die((void *)pparams, sizeof(*pparams), 1, fp);
	check_magic_number(pparams->magic_number,  "magic number",
		"get_mcmc_params_from_file");
	check_magic_number(pparams->magic_number2, "magic number 2",
		"get_mcmc_params_from_file");
}

// ----------------------------------------------------------------
void populate_footer(mcmc_footer_t* pfooter,
	double start_time, double therm_time, double end_time,
	metro_stats_t* pmetro_stats)
{
	memset((void*)pfooter, 0, sizeof(*pfooter));
	pfooter->magic_number  = RCM_MAGIC_NUMBER;
	pfooter->start_time    = start_time;
	pfooter->therm_time    = therm_time;
	pfooter->end_time      = end_time;
	pfooter->metro_stats   = *pmetro_stats;
	pfooter->magic_number2 = RCM_MAGIC_NUMBER;
}

// ----------------------------------------------------------------
void read_footer(mcmc_footer_t* pfooter, FILE* fp)
{
	read_from_file_or_die((void*)pfooter, sizeof(*pfooter), 1, fp);
	check_magic_number(pfooter->magic_number,  "magic number",   "read_footer");
	check_magic_number(pfooter->magic_number2, "magic number 2", "read_footer");
}

// ----------------------------------------------------------------
void write_footer(mcmc_footer_t* pfooter, FILE* fp)
{
	write_to_file_or_die((void*)pfooter, sizeof(*pfooter), 1, fp);
}

// ----------------------------------------------------------------
void display_footer(mcmc_footer_t* pfooter)
{
	// Display Metropolis stats.
	report_metro_stats(&pfooter->metro_stats);

	// Display elapsed time.
	printf("# Elapsed thermalization seconds: %.6lf\n",
		pfooter->therm_time - pfooter->start_time);
	printf("# Elapsed accumulation   seconds: %.6lf\n",
		pfooter->end_time - pfooter->therm_time);
	printf("# Elapsed total seconds:          %.6lf\n",
		pfooter->end_time - pfooter->start_time);
}

// ----------------------------------------------------------------
static void check_magic_number(unsigned file_value,
	char* var_name, char* func_name)
{
	if (file_value != RCM_MAGIC_NUMBER) {
		fprintf(stderr,
		"%s:  expected %s 0x%08x; got 0x%08x.\n",
			func_name, var_name,
			RCM_MAGIC_NUMBER, file_value);
		fprintf(stderr, "Incorrect file format, or file format has changed?\n");
		exit(1);
	}
}
