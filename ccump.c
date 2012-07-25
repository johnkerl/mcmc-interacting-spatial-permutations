// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//
// Code for experimenting with the _correlated uniform Markov process_  as
// defined in Appendix B of http://johnkerl.org/rcm/kerl-dis-gc-ssp.pdf
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cump_params.h"
#include "rcmrand.h"
#include "util.h"

// ----------------------------------------------------------------
// Parameter defaults.
#define DEFAULT_SEED                  0
#define DEFAULT_DO_SEED               0

#define DEFAULT_ETA                   0.99

#define DEFAULT_NY                    10

#define DEFAULT_RAND_Y0               0
#define DEFAULT_TERMINATION_NUM_ACC   10000
#define DEFAULT_THERM_NUM_ITER        10000

// ----------------------------------------------------------------
// Parameters.
typedef struct _cump_params_t {
	int      rng_type;
	unsigned seed;
	int      do_seed;

	double  eta;

	double  etac;
	double  a;
	double  b;
	double  s;
	double  width;

	int     NY;

	int     rand_Y0;
	int     termination_num_acc;
	int     therm_num_iter;

} cump_params_t;

// ----------------------------------------------------------------
// Prototypes for local functions.

static void initialize_Ys(double* Ys, cump_params_t* pparams);
static void update_Ys    (double* Ys, cump_params_t* pparams);
static void print_Ys     (double* Ys, cump_params_t* pparams);

static void usage(char * argv0);
static int parse_command_line(int argc, char** argv, cump_params_t* pparams);
static void set_default_cump_params(cump_params_t* p);
static void print_cump_params(cump_params_t* p);

// ================================================================
int main(int argc, char **argv)
{
	double* Ys = 0;
	cump_params_t params;
	int i;

	// Set parameters.
	set_default_cump_params(&params);
	if (!parse_command_line(argc, argv, &params))
		usage(argv[0]);

	// Seed the pseudorandom-number generator.
	if (params.do_seed)
		SRANDOM(params.seed);
	else
		STRANDOM();

	// Allocate and initialize memory.
	Ys = double_malloc_or_die(params.NY);

	// Print some information.
	print_cump_params(&params);
	fflush(stdout);

	// Choose the initial values.
	initialize_Ys(Ys, &params);

	// Thermalization loop
	for (i = 0; i < params.therm_num_iter; i++) {
		update_Ys(Ys, &params);
	}

	// Accumulation loop
	for (i = 0; i < params.termination_num_acc; i++) {
		update_Ys(Ys, &params);
		print_Ys (Ys, &params);
	}

	// Display means of random variables.
	//if (lparams.compute_rvs) {
	//	compute_rv_means(&rvs, ppoints);
	//	display_rv_means(&rvs, &params, ppoints, &verbosity);
	//}

	// Free memory.  This is done at exit, of course, but it's good hygiene in
	// case this code is copied and pasted into something larger:  one wants
	// one's code, as well as its offspring, to be free of memory leaks.
	// Also, explicitly freeing memory allows the valgrind tool to tell us
	// what (if anything) we missed freeing.
	free(Ys);

	return 0;
}

// ----------------------------------------------------------------
static void initialize_Ys(double* Ys, cump_params_t* pparams)
{
	int j;
	double width = pparams->b - pparams->a;

	if (pparams->rand_Y0) {
		for (j = 0; j < pparams->NY; j++)
			Ys[j] = pparams->a + width * URANDOM();
	}
	else {
		for (j = 0; j < pparams->NY; j++)
			Ys[j] = pparams->a;
	}
}

// ----------------------------------------------------------------
static void update_Ys(double* Ys, cump_params_t* pparams)
{
	double U;
	int j;

	for (j = 0; j < pparams->NY; j++) {
		U = pparams->a + pparams->width * URANDOM();
		Ys[j] = pparams->eta * Ys[j] + pparams->etac * U;
	}
}

// ----------------------------------------------------------------
static void print_Ys(double* Ys, cump_params_t* pparams)
{
	int j;

	for (j = 0; j < pparams->NY; j++) {
		if (j > 0)
			printf(" ");
		printf("%11.7f", Ys[j]);
	}
	printf("\n");
}

// ----------------------------------------------------------------
static void usage(char * argv0)
{
	fprintf(stderr, "Usage: %s [options]\n", argv0);
	fprintf(stderr, "Options:\n");

	fprintf(stderr,

	"-h, --help:     Print this message.\n"
	"seed=[number]:  Specify pseudorandom seed (default: use time/PID value).\n"
	"\n"

	"randinit=[0,1]: Y0 is a (0) or uniform random (1).\n"
	"\n"

	"eta=[number]:   Specify eta.  Must be >= 0 and < 1.\n"
	"\n"

	"nthm=[number]:  Specify number of thermalization iterations.\n"
	"nacc=[number]:  Specify number of batches.\n"
	"\n");

	exit(1);
}

// ----------------------------------------------------------------
static int parse_command_line(int argc, char** argv, cump_params_t* pparams)
{
	int argi;

	for (argi = 1; argi < argc; argi++) {
		if (strcmp(argv[argi], "-h") == 0)
			return 0;
		else if (strcmp(argv[argi], "--help") == 0)
			return 0;

		else if (sscanf(argv[argi], "seed=%i", &pparams->seed) == 1)
			pparams->do_seed = 1;

		else if (sscanf(argv[argi], "ri=%d", &pparams->rand_Y0) == 1)
			;

		else if (sscanf(argv[argi], "eta=%lf", &pparams->eta) == 1)
			;

		else if (sscanf(argv[argi], "nthm=%d", &pparams->therm_num_iter) == 1)
			;
		else if (sscanf(argv[argi], "nacc=%d", &pparams->termination_num_acc) == 1)
			;

		else
			return 0;
	}

	if ((pparams->eta < 0.0) || (pparams->eta >= 1.0))
		return 0;

	pparams->etac = 1.0 - pparams->eta;
	pparams->s = sqrt((1 + pparams->eta) / (1 - pparams->eta));
	pparams->b = 0.5 * (1 + pparams->s);
	pparams->a = 0.5 * (1 - pparams->s);
	pparams->width = pparams->b - pparams->a;

	return 1;
}

// ----------------------------------------------------------------
static void set_default_cump_params(cump_params_t* p)
{
	memset((void*)p, 0xee, sizeof(*p));

	p->rng_type                = RCM_WHICH;
	p->seed                    = DEFAULT_SEED;
	p->do_seed                 = DEFAULT_DO_SEED;

	p->eta                     = DEFAULT_ETA;

	p->NY                      = DEFAULT_NY;

	p->rand_Y0                 = DEFAULT_RAND_Y0;
	p->termination_num_acc     = DEFAULT_TERMINATION_NUM_ACC;
	p->therm_num_iter          = DEFAULT_THERM_NUM_ITER;
}

// ----------------------------------------------------------------
static void print_cump_params(cump_params_t* p)
{
	printf("# RNG = %s\n", get_rng_type_desc(p->rng_type));
	if (p->do_seed)
		printf("# RNG seed = 0x%08x\n", p->seed);
	printf("# eta = %11.7lf\n", p->eta);

	printf("# %d independent Y_t sequence%s.\n", p->NY, p->NY == 1 ? "" : "s");
	printf("# Y0: %s.\n", p->rand_Y0 ? "uniform random" : "a");

	printf("# Thermalization is declared after %d iteration%s.\n",
		p->therm_num_iter,
		p->therm_num_iter == 1 ? "" : "s");

	printf("# Terminate after %d accumulation%s.\n",
		p->termination_num_acc,
		p->termination_num_acc == 1 ? "" : "s");
}
