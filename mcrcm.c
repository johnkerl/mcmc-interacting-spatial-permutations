// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mcmc_params.h"
#include "mcmc_rvs.h"
#include "thermalization.h"
#include "rcmrand.h"
#include "util.h"
#include "points.h"
#include "energy.h"
#include "metropolis.h"
#include "times.h"
#include "seps.h"
#include "dotplot.h"

// ----------------------------------------------------------------
typedef struct _local_params_t {
	char* pmts_fn;
	FILE* pmts_fp;
	int compute_rvs;
	int dry_run;
	int print_pre_therm_RVs;
	int H_verbosity;
	int do_print_worm;
	int do_print_wormspan;
} local_params_t;

static void do_one_sweep(points_t* ppoints, mcmc_params_t* pparams,
	local_params_t* plparams, int sweepi, int thermalized, therm_ctl_t* ptherm,
	mcmc_rvs_t* prvs, rv_verbosity_t* pverbosity, metro_stats_t* pmetro_stats,
	double* alphas, double* pH, double* pD, double* pV);

static int parse_command_line(int argc, char** argv,
	local_params_t* plparams, mcmc_params_t* pparams,
	mcmc_rvs_t* prvs, rv_verbosity_t* pverbosity);

static void usage(char * argv0);

// ----------------------------------------------------------------
int main(int argc, char **argv)
{
	// Lattice points and permutations:
	points_t* ppoints = 0;
	// MCMC parameters:
	mcmc_params_t params;
	mcmc_footer_t footer;
	double* alphas;

	local_params_t lparams;
	double  start_time, therm_time, end_time;

	// Random variables:
	double         H = 0.0, D = 0.0, V = 0.0;
	mcmc_rvs_t     rvs;
	rv_verbosity_t verbosity;

	// Metropolis control variables:
	int            acci, sweepi = 0;

	// Metropolis statistics:
	metro_stats_t  metro_stats;

	// Thermalization-detection variables:
	therm_ctl_t therm;


	// Set parameters.
	set_default_mcmc_params(&params);
	memset(&verbosity, 0, sizeof(verbosity));
	memset(&lparams,   0, sizeof(lparams));
	lparams.compute_rvs = 1;

	if (!parse_command_line(argc, argv, &lparams, &params, &rvs, &verbosity)) {
		usage(argv[0]);
	}

	// Seed the pseudorandom-number generator.
	if (params.do_seed)
		SRANDOM(params.seed);
	else
		STRANDOM();
	start_time = get_sys_time_double();

	// Allocate and initialize memory.
	ppoints = get_cubic_lattice_points(params.L, params.d);
	alphas = get_alphas(params.interaction_type, ppoints->N, params.alpha0);
	allocate_rvs(&rvs, &params, ppoints);

	// Initialize the smoothing window.
	therm_ctl_init(&therm, params.smoothing_window_size,
		params.therm_ntp, params.force_therm);

	// Print some information.
	print_mcmc_params(&params);
	fflush(stdout);
	if (lparams.dry_run) {
		printf("# Dry run only.  Exiting.\n");
		exit(0);
	}

	// Write MCMC parameters (file-header information) to the realization file.
	if (lparams.pmts_fn) {
		lparams.pmts_fp = open_file_or_die(lparams.pmts_fn, "wb");
		write_to_file_or_die((void*)&params, sizeof(params), 1,
			lparams.pmts_fp);
	}

	// Choose the initial permutation.
	if (params.rand_init_pi)
		set_unif_rand_pmt(ppoints);
	else {
		if (params.wx1_init_pi)
			set_wx_pmt(ppoints);
		if (params.wy1_init_pi)
			set_wy_pmt(ppoints);
		if (params.wz1_init_pi)
			set_wz_pmt(ppoints);
	}
	// else, leave it at the identity permutation.

	set_up_cycinfo_list(ppoints);
#ifdef CHECK_CYCINFO
	sanity_check_cycinfo_list(ppoints);
#endif

	H = get_H_of_pi(ppoints, params.beta, alphas,
		params.interaction_type, &D, &V);

	// Make one sweep to exclude the initial zero-energy state from the
	// smoothing window.
	printf("#\n");
	printf("# Initial HDV = %11.7lf %11.7lf %11.7lf\n", H, D, V);
	do_one_sweep(ppoints, &params, &lparams, sweepi++,
		therm.thermalized, &therm,
		&rvs, &verbosity, &metro_stats, alphas, &H, &D, &V);
	printf("# Sweep 0 HDV = %11.7lf %11.7lf %11.7lf\n", H, D, V);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Thermalization loop
	clear_metro_stats(&metro_stats);
	for (; !therm.thermalized; sweepi++) {
		do_one_sweep(ppoints, &params, &lparams, sweepi,
			therm.thermalized, &therm,
			&rvs, &verbosity, &metro_stats, alphas, &H, &D, &V);
		if (lparams.print_pre_therm_RVs)
			display_realization_rvs(H, D, V, &rvs, ppoints, &verbosity);
	}
	therm_time = get_sys_time_double();
	printf("# Thermalization complete:  sweep %6d H = %11.7lf\n", sweepi, H);
	printf("# Ntherm = %6d\n", sweepi);
	fflush(stdout);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Accumulation loop
	clear_metro_stats(&metro_stats);
	for (acci = 0; acci < params.termination_num_acc; acci++, sweepi++) {

		//if ((acci % 1000) == 0) {
		//	// ctime() includes a trailing carriage return.
		//	printf("# acci = %7d  H = %11.7lf  %s",
		//		acci, H, get_sys_time_string());
		//	fflush(stdout);
		//}

		do_one_sweep(ppoints, &params, &lparams, sweepi,
			therm.thermalized, &therm,
			&rvs, &verbosity, &metro_stats, alphas, &H, &D, &V);
		if (lparams.compute_rvs)
			accumulate_rvs(&rvs, ppoints, &params, H, D, V);
		display_realization_rvs(H, D, V, &rvs, ppoints, &verbosity);
	}
	end_time = get_sys_time_double();
	printf("# Accumulation complete:  acc %6d H = %11.7lf\n",
		params.termination_num_acc, H);
	printf("#\n");
	fflush(stdout);

	//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Display means of random variables.
	if (lparams.compute_rvs) {
		compute_rv_means(&rvs, ppoints);
		display_rv_means(&rvs, &params, ppoints, &verbosity);
	}

	// Write MCMC footer to standard output, and to the realization file.
	populate_footer(&footer, start_time, therm_time, end_time, &metro_stats);
	display_footer(&footer);
	if (lparams.pmts_fn)
		write_footer(&footer, lparams.pmts_fp);

	// Close the realization file.
	if (lparams.pmts_fp) {
		close_file_or_die(lparams.pmts_fp);
		lparams.pmts_fp = 0;
	}

	if (lparams.pmts_fn)
		printf("# Wrote realization file \"%s\".\n", lparams.pmts_fn);

	// Free memory.  This is done at exit, of course, but it's good hygiene in
	// case this code is copied and pasted into something larger:  one wants
	// one's code, as well as its offspring, to be free of memory leaks.
	// Also, explicitly freeing memory allows the valgrind tool to tell us
	// what (if anything) we missed freeing.
	free_rvs(&rvs);
	free_points(ppoints);
	free(alphas);
	therm_ctl_free(&therm);

	return 0;
}

// ----------------------------------------------------------------
static void do_one_sweep(points_t* ppoints, mcmc_params_t* pparams,
	local_params_t* plparams, int sweepi, int thermalized, therm_ctl_t* ptherm,
	mcmc_rvs_t* prvs, rv_verbosity_t* pverbosity, metro_stats_t* pmetro_stats,
	double* alphas, double* pH, double* pD, double* pV)
{
	char* tadesc = thermalized ? "accum" : "therm";
	int i;

	if (pverbosity->sweepi == 2) {
		printf("\n\n%s", LARGE_SEP);
		printf("# %s phase sweep %d\n\n", tadesc, sweepi);
		fflush(stdout);
	}
	if (pverbosity->sweepi == 1) {
		printf("# %s phase sweep %d\n", tadesc, sweepi);
		fflush(stdout);
	}

	// Quit if the file "__stop__" is present in the current directory.
	// (This makes it easy to kill a multiprocessor job.)
	check_stop();

	for (i = 0; i < pparams->gsweeps_per_acc; i++)
		SO_sweep(ppoints, pparams, alphas, pH, pD, pV, pmetro_stats);
	if (pparams->allow_cycle_reversals)
		R_sweep(ppoints, pparams);
	for (i = 0; i < pparams->wsweeps_per_acc; i++)
		worm_sweep(ppoints, pparams, plparams->do_print_worm,
			plparams->do_print_wormspan, alphas, pH, pD, pV, pmetro_stats);

	// Write the current permutation to the realization file.
	if (plparams->pmts_fp && thermalized)
		write_pmt_to_file(ppoints, plparams->pmts_fp);

	check_H(ppoints, pparams->beta, alphas, pparams->interaction_type,
		*pH, *pD, *pV);
	sanity_check_cycinfo_list(ppoints);

	// Update the smoothing window for thermalization detection.
	therm_ctl_update(ptherm, *pH);

	// Display system energy, unsmoothed as well as smoothed.
	if (plparams->H_verbosity >= 1) {
		if (!thermalized) {
			if (plparams->print_pre_therm_RVs) {
				printf("%5d %11.7lf %11.7lf %11.7lf %11.7lf %11.7lf # %d H D V\n",
					sweepi, *pH, *pD, *pV,
					ptherm->ssH, thermalized * ptherm->ssH, thermalized);
			}
		}
		else {
			printf("%5d %11.7lf %11.7lf %11.7lf %11.7lf %11.7lf # %d H\n",
				sweepi, *pH, *pD, *pV, *pH, thermalized * *pH,
				thermalized);
		}
		fflush(stdout);
	}
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

	"L=[number]:     Specify extent of cubic lattice.\n"
	"d=[1,2,3]:      Specify dimension.\n"
	"\n"

	"randinit=[0,1]: Initial pi is id (0) or uniform random (1).\n"
	"\n"

	"T=[number]:     Specify T (temperature).\n"
	"beta=[number]:  Specify beta (inverse temperature).\n"
	"\n"

	"alpha=[number]: Specify alpha (for cycle interactions).\n"
	"gamma=[number]: Specify gamma (open-worm weight).\n"
	"none:           No interactions.\n"
	"r2:             r_2 interactions.\n"
	"rell:           r_ell interactions.\n"
	"vzwl:           V(z,w) lower interactions.\n"
	"vzwu:           V(z,w) upper interactions.\n"
	"\n"

	"gspa=[number]:  Specify number of SO   sweeps per acc.\n"
	"wspa=[number]:  Specify number of worm sweeps per acc.\n"
	"cycrev=[0,1]:   Allow cycle reversals.\n"
	"\n"

	"nsm=[number]:   Specify size of H-smoothing window.\n"
	"ntp=[number]:   The system is deemed thermalized when the smoothed energy\n"
	"                reaches this many turning points.\n"
	"ft=[number]:    Force thermalization after at most [number] sweeps.\n"
	"\n"

	"nacc=[number]:  Specify number of accumulations.\n"
	"rndswp=[0,1]:   Randomize sites for Metropolis sweeps.\n"
	"\n"

	"pw=[0,1]:       Print worm sites on each MCMC step.\n"
	"pws=[0,1]:      Print worm span on each MCMC step.\n"
	"\n"

	"siv=[0,1,2]:    Print sweep number.\n"
	"pt=[0,1]:       Print pre-thermalization random-variable values.\n"
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
	"anglev=[0,1]:   Specify verbosity for angle cosines.\n"
	"\n"
	"rzn={realization file name}: Write out a realization file for use by rvrcm.\n"
	"crv=[0,1]:      Compute random variables or not, e.g. crv=0 along with rzn=...\n"
	"\n"
	"dryrun=[0,1]:   Do a dry run:  print header and footer with no MCMC steps.\n");

	exit(1);
}

// ----------------------------------------------------------------
static int parse_command_line(int argc, char** argv,
	local_params_t* plparams, mcmc_params_t* pparams,
	mcmc_rvs_t* prvs, rv_verbosity_t* pv)
{
	int argi;

	for (argi = 1; argi < argc; argi++) {
		if (strcmp(argv[argi], "-h") == 0)
			return 0;
		else if (strcmp(argv[argi], "--help") == 0)
			return 0;

		else if (sscanf(argv[argi], "dryrun=%d", &plparams->dry_run) == 1)
			;
		else if (strcmp(argv[argi], "dryrun") == 0)
			plparams->dry_run = 1;
		else if (sscanf(argv[argi], "pt=%d", &plparams->print_pre_therm_RVs) == 1)
			;
		else if (sscanf(argv[argi], "hv=%d", &plparams->H_verbosity) == 1)
			;
		else if (sscanf(argv[argi], "pw=%d", &plparams->do_print_worm) == 1)
			;
		else if (sscanf(argv[argi], "pws=%d", &plparams->do_print_wormspan)==1)
			;
		else if (strcmp(argv[argi], "pw") == 0)
			plparams->do_print_worm = 1;
		else if (strcmp(argv[argi], "pws") == 0)
			plparams->do_print_wormspan = 1;
		else if (strncmp(argv[argi], "rzn=", strlen("rzn=")) == 0)
			plparams->pmts_fn = &argv[argi][strlen("rzn=")];
		else if (sscanf(argv[argi], "crv=%d", &plparams->compute_rvs) == 1)
			;

		else if (sscanf(argv[argi], "seed=%i", &pparams->seed) == 1)
			pparams->do_seed = 1;

		else if (sscanf(argv[argi], "L=%d", &pparams->L) == 1)
			;
		else if (sscanf(argv[argi], "d=%d", &pparams->d) == 1)
			;

		else if (sscanf(argv[argi], "ri=%d", &pparams->rand_init_pi) == 1)
			;
		else if (sscanf(argv[argi], "wxi=%d", &pparams->wx1_init_pi) == 1)
			;
		else if (sscanf(argv[argi], "wyi=%d", &pparams->wy1_init_pi) == 1)
			;
		else if (sscanf(argv[argi], "wzi=%d", &pparams->wz1_init_pi) == 1)
			;

		else if (sscanf(argv[argi], "T=%lf", &pparams->beta) == 1)
			pparams->beta = 1.0 / pparams->beta;
		else if (sscanf(argv[argi], "beta=%lf", &pparams->beta) == 1)
			;
		else if (sscanf(argv[argi], "alpha0=%lf", &pparams->alpha0) == 1)
			;

		else if (sscanf(argv[argi], "none_alpha0=%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = NO_INTERACTIONS;
		else if (sscanf(argv[argi], "r2_alpha0=%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = R2_INTERACTIONS;
		else if (sscanf(argv[argi], "rell_alpha0=%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = RELL_INTERACTIONS;

		else if (sscanf(argv[argi], "none_alpha0_%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = NO_INTERACTIONS;
		else if (sscanf(argv[argi], "r2_alpha0_%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = R2_INTERACTIONS;
		else if (sscanf(argv[argi], "rell_alpha0_%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = RELL_INTERACTIONS;

		else if (sscanf(argv[argi], "alpha=%lf", &pparams->alpha0) == 1)
			pparams->interaction_type = R2_INTERACTIONS;
		else if (sscanf(argv[argi], "gamma=%lf", &pparams->gamma) == 1)
			;


		else if (sscanf(argv[argi], "cycrev=%d", &pparams->allow_cycle_reversals) == 1)
			;
		else if (sscanf(argv[argi], "gspa=%d", &pparams->gsweeps_per_acc) == 1)
			;
		else if (sscanf(argv[argi], "wspa=%d", &pparams->wsweeps_per_acc) == 1)
			;

		else if (sscanf(argv[argi], "hyb=%d,%d",
			&pparams->gsweeps_per_acc, &pparams->wsweeps_per_acc) == 2)
		{
			pparams->allow_cycle_reversals = 1;
		}

		else if (strcmp(argv[argi], "none") == 0)
			pparams->interaction_type = NO_INTERACTIONS;
		else if (strcmp(argv[argi], "r2") == 0)
			pparams->interaction_type = R2_INTERACTIONS;
		else if (strcmp(argv[argi], "rell") == 0)
			pparams->interaction_type = RELL_INTERACTIONS;
		else if (strcmp(argv[argi], "vzwl") == 0)
			pparams->interaction_type = VZWL_INTERACTIONS;
		else if (strcmp(argv[argi], "vzwu") == 0)
			pparams->interaction_type = VZWU_INTERACTIONS;

		else if (sscanf(argv[argi], "nsm=%d",
			&pparams->smoothing_window_size) == 1)
			;
		else if (sscanf(argv[argi], "ntp=%d", &pparams->therm_ntp) == 1)
			;

		else if (sscanf(argv[argi], "ft=%d", &pparams->force_therm) == 1)
			;

		else if (sscanf(argv[argi], "nacc=%d", &pparams->termination_num_acc) == 1)
			;
		else if (sscanf(argv[argi], "rndswp=%d", &pparams->do_random_sweep)==1)
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
	return 1;
}
