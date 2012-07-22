// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef MCMC_PARAMS_H
#define MCMC_PARAMS_H

#include <stdio.h>
#include "interactions.h"
#include "metro_stats.h"

// ----------------------------------------------------------------
#define RCM_MAGIC_NUMBER              0x214d4352 // "RCM!"

#define DEFAULT_SEED                  0
#define DEFAULT_DO_SEED               0

#define DEFAULT_L                     10
#define DEFAULT_D                     3

#define DEFAULT_BETA                  (1.0/6.8)
#define DEFAULT_ALPHA0                0.0
#define DEFAULT_INTERACTION_TYPE      NO_INTERACTIONS
#define DEFAULT_GAMMA                 0.05

#define DEFAULT_RAND_INIT_PI          0
#define DEFAULT_WX1_INIT_PI           0
#define DEFAULT_WY1_INIT_PI           0
#define DEFAULT_WZ1_INIT_PI           0
#define DEFAULT_DO_RANDOM_SWEEP       0
#define DEFAULT_TERMINATION_NUM_ACC  10000

// Default:  SAR rather than SO, worm, or hybrid.
#define DEFAULT_ALLOW_CYCLE_REVERSALS 1
#define DEFAULT_GSWEEPS_PER_ACC       1
#define DEFAULT_WSWEEPS_PER_ACC       0

#define DEFAULT_FORCE_THERM           -1
#define DEFAULT_SMOOTHING_WINDOW_SIZE 200
#define DEFAULT_THERM_NTP             40

// ----------------------------------------------------------------
typedef struct _mcmc_params_t {

	unsigned magic_number; // Indicates file type.

	int      rng_type;
	unsigned seed;
	int      do_seed;

	// Data dimensions.
	int     L;
	int     d;

	// Temperature and interaction parameters.
	double  beta;
	double  alpha0;
	int     interaction_type;
	// Pad is for identical struct packing on 32-bit & 64-bit processors.
	int     pad;
	double  gamma;

	int     tbd1;
	int     tbd2;
	int     tbd3;
	int     tbd4;

	// Metropolis parameters
	int     rand_init_pi;
	int     wx1_init_pi;
	int     wy1_init_pi;
	int     wz1_init_pi;
	int     do_random_sweep;
	int     allow_cycle_reversals;
	int     gsweeps_per_acc;
	int     wsweeps_per_acc;
	int     termination_num_acc;

	// Thermalization parameters
	int     force_therm;
	int     smoothing_window_size;
	int     therm_ntp;

	// Helps to detect changes in this header structure:  E.g. if
	// a new program is reading an old file.
	unsigned magic_number2;

} mcmc_params_t;

void print_mcmc_params(mcmc_params_t* p);
void set_default_mcmc_params(mcmc_params_t* p);
void get_mcmc_params_from_file(mcmc_params_t* pparams, FILE* fp);

// ----------------------------------------------------------------
typedef struct _mcmc_footer_t {

	// Helps to detect changes in this header structure:  E.g. if
	// a new program is reading an old file.
	unsigned magic_number;

	// Pad is for identical struct packing on 32-bit & 64-bit processors.
	int     pad;

	double start_time;
	double therm_time;
	double end_time;
	metro_stats_t metro_stats;

	int     tbd[4];

	// Helps to detect changes in this header structure:  E.g. if
	// a new program is reading an old file.
	unsigned magic_number2;

} mcmc_footer_t;

void populate_footer(mcmc_footer_t* pfooter,
	double start_time, double therm_time, double end_time,
	metro_stats_t* pmetro_stats);
void write_footer(mcmc_footer_t* pfooter, FILE* fp);
void read_footer(mcmc_footer_t* pfooter, FILE* fp);
void display_footer(mcmc_footer_t* pfooter);

#endif // MCMC_PARAMS_H
