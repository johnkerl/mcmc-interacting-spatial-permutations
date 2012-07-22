// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mcmc_params.h"
#include "mcmc_rvs.h"
#include "rcmrand.h"
#include "util.h"
#include "interactions.h"
#include "energy.h"
#include "rho.h"
#include "vector_misc.h"
#include "stats.h"
#include "times.h"
#include "seps.h"

// ----------------------------------------------------------------
typedef struct _rv_local_params_t {
	char* pmts_fn;
	FILE* pmts_fp;
	int   file_info_only;
} rv_local_params_t;

static int parse_command_line(int argc, char** argv,
	rv_local_params_t* plparams, mcmc_rvs_t* prvs, rv_verbosity_t* pv);
static void usage(char * argv0);

// ----------------------------------------------------------------
int main(int argc, char **argv)
{
	// Lattice points and permutations:
	points_t*     ppoints = 0;

	// MCMC parameters:
	mcmc_params_t  params;
	mcmc_footer_t  footer;
	double*        alphas;
	rv_local_params_t lparams;
	int            sweeps_per_acc = 0;

	// Random variables:
	mcmc_rvs_t     rvs;
	rv_verbosity_t verbosity;
	double         H = 0.0;
	double         D = 0.0;
	double         V = 0.0;

	// Loop variable:
	int acci, sweepi = 0;


	// Set parameters.
	lparams.pmts_fn = 0;
	lparams.pmts_fp = 0;
	lparams.file_info_only = 0;
	memset(&verbosity, 0, sizeof(verbosity));
	if (!parse_command_line(argc, argv, &lparams, &rvs, &verbosity))
		usage(argv[0]);

	// Open the realization file.  Read MCMC parameters from the header.
	lparams.pmts_fp = open_file_or_die(lparams.pmts_fn, "rb");
	get_mcmc_params_from_file(&params, lparams.pmts_fp);
	sweeps_per_acc = params.gsweeps_per_acc + params.wsweeps_per_acc;
	printf("# Using realization file \"%s\".\n", lparams.pmts_fn);

	// Display the MCMC parameters.
	print_mcmc_params(&params);
	if (lparams.file_info_only) {
		close_file_or_die(lparams.pmts_fp);
		exit(0);
	}

	// Allocate memory.
	ppoints = get_cubic_lattice_points(params.L, params.d);
	set_up_cycinfo_list(ppoints);
#ifdef CHECK_CYCINFO
	sanity_check_cycinfo_list(ppoints);
#endif
	alphas = get_alphas(params.interaction_type, ppoints->N, params.alpha0);
	allocate_rvs(&rvs, &params, ppoints);


	for (acci = 0; acci < params.termination_num_acc; acci++) {
		if (verbosity.sweepi) {
			printf("\n\n%s", LARGE_SEP);
			printf("# realization %d\n\n", sweepi);
		}

		// Read the next permutation from the file.
		read_pmt_from_file(ppoints, lparams.pmts_fp);
		recompute_cycinfo_list(ppoints);

#ifdef CHECK_CYCINFO
		sanity_check_cycinfo_list(ppoints);
#endif

		H = get_H_of_pi(ppoints, params.beta, alphas,
			params.interaction_type, &D, &V);

		if (verbosity.H)
			printf("%5d %11.7lf %11.7lf %11.7lf # H\n", sweepi, H, D, V);
		display_realization_rvs(H, D, V, &rvs, ppoints, &verbosity);
		accumulate_rvs(&rvs, ppoints, &params, H, D, V);
	}

	compute_rv_means(&rvs, ppoints);
	display_rv_means(&rvs, &params, ppoints, &verbosity);

	read_footer(&footer, lparams.pmts_fp);
	display_footer(&footer);

	// Free memory.  This is done at exit, of course, but it's good hygiene in
	// case this code is copied and pasted into something larger:  one wants
	// one's code, as well as its offspring, to be free of memory leaks.
	free_rvs(&rvs);
	free_points(ppoints);
	free(alphas);
	close_file_or_die(lparams.pmts_fp);

	return 0;
}

// ----------------------------------------------------------------
static int parse_command_line(int argc, char** argv,
	rv_local_params_t* plparams, mcmc_rvs_t* prvs, rv_verbosity_t* pv)
{
	int argi;

	for (argi = 1; argi < argc; argi++) {
		if (strcmp(argv[argi], "-h") == 0)
			return 0;
		else if (strcmp(argv[argi], "--help") == 0)
			return 0;

		else if (strncmp(argv[argi], "rzn=", strlen("rzn=")) == 0)
			plparams->pmts_fn = &argv[argi][strlen("rzn=")];
		else if (strcmp(argv[argi], "info") == 0)
			plparams->file_info_only = 1;
		else if (sscanf(argv[argi], "info=%d", &plparams->file_info_only) == 1)
			;

		else if (sscanf(argv[argi], "siv=%d",    &pv->sweepi)         == 1)
			;
		else if (sscanf(argv[argi], "rhov=%d",   &pv->rho)            == 1)
			;
		else if (sscanf(argv[argi], "pmtv=%d",   &pv->pmt)            == 1)
			;
		else if (sscanf(argv[argi], "cv=%d",     &pv->counts)         == 1)
			;
		else if (sscanf(argv[argi], "cdv=%d",    &pv->cycdec)         == 1)
			;
		else if (sscanf(argv[argi], "hv=%d",     &pv->H)              == 1)
			;
		else if (sscanf(argv[argi], "tiv=%d",    &pv->tauint_Hs)      == 1)
			;
		else if (sscanf(argv[argi], "r2v=%d",    &pv->r2)             == 1)
			;
		else if (sscanf(argv[argi], "lmaxv=%d",  &pv->lmax)           == 1)
			;
		else if (sscanf(argv[argi], "ellv=%d",   &pv->ell)            == 1)
			;
		else if (sscanf(argv[argi], "ellbarv=%d",&pv->ellbar)         == 1)
			;
		else if (sscanf(argv[argi], "c0sv=%d",   &pv->cycle_0_spatlen)== 1)
			;
		else if (sscanf(argv[argi], "mslv=%d",   &pv->meanspatlen)    ==1)
			;
		else if (sscanf(argv[argi], "wcv=%d",    &pv->winding_cycle)  == 1)
			;
		else if (sscanf(argv[argi], "wnv=%d",    &pv->winding_number) == 1)
			;
		else if (sscanf(argv[argi], "jlv=%d",    &pv->jump_length)    == 1)
			;
		else if (sscanf(argv[argi], "jtv=%d",    &pv->jump_target)    == 1)
			;
		else if (sscanf(argv[argi], "jvv=%d",    &pv->jump_vector)    == 1)
			;
		else if (sscanf(argv[argi], "kinkdv=%d", &pv->kink_dot)       == 1)
			;
		else if (sscanf(argv[argi], "kinkav=%d", &pv->kink_cos_angle) == 1)
			;
		else if (sscanf(argv[argi], "anglev=%d", &pv->cos_angle)      == 1)
			;

		else
			return 0;
	}
	if (plparams->pmts_fn == 0)
		return 0;
	return 1;
}

// ----------------------------------------------------------------
static void usage(char * argv0)
{
	fprintf(stderr, "Usage: %s [options] rzn={realization file name}\n", argv0);
	fprintf(stderr, "Options:\n");

	fprintf(stderr,

	"-h, --help:     Print this message.\n"
	"seed=[number]:  Specify pseudorandom seed (default: use time/PID value).\n"
	"\n"

	"siv=[0,1,2]:    Print sweep number.\n"
	"rhov=[0,1]:     Specify rho verbosity.\n"
	"pmtv=[0,1]:     Specify permutation verbosity.\n"
	"cv=[0,1]:       Specify cycle-count (r_ell) verbosity.\n"
	"cdv=[0,1,2]:    Specify cycle-decomposition verbosity.\n"
	"hv=[0,1]:       Specify energy verbosity.\n"
	"tiv=[0,1]:      Specify tauint verbosity.\n"
	"r2v=[0,1]:      Specify r2 verbosity.\n"
	"lmaxv=[0,1]:    Specify ell_max verbosity.\n"
	"ellv=[0,1]:     Specify ell(x) verbosity.\n"
	"ellbarv=[0,1]:  Specify ellbar verbosity.\n"
	"c0sv=[number]:  Specify verbosity for cycle 0 spatial length.\n"
	"mslv=[number]:  Specify verbosity for mean spatial length.\n"
	"wcv=[0,1]:      Specify winding-cycle verbosity.\n"
	"wnv=[0,1,2]:    Specify winding-number verbosity.\n"
	"jlv=[0,1]:      Specify verbosity for jump lengths.\n"
	"jtv=[0,1]:      Specify verbosity for jump targets.\n"
	"jvv=[0,1]:      Specify verbosity for jump vectors.\n"
	"kinkdv=[0,1]:   Specify verbosity for kink dot products.\n"
	"kinkav=[0,1]:   Specify verbosity for kink angle cosines.\n"
	"anglev=[0,1]:   Specify verbosity for angle cosines.\n");

	exit(1);
}
