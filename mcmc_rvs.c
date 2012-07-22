// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "mcmc_rvs.h"
#include "rho.h"
#include "pmt.h"
#include "stats.h"
#include "vector_misc.h"

#define SMALL_SEP "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"

// ----------------------------------------------------------------
void allocate_rvs(mcmc_rvs_t* prvs, mcmc_params_t* pparams,
	points_t* ppoints)
{
	int N = ppoints->N;
	int k;
	int midx    = ppoints->dims.nx/2;
	int midy    = ppoints->dims.ny/2;
	int midz    = ppoints->dims.nz/2;
	int num_acc = pparams->termination_num_acc;

	prvs->unitk         = double_malloc_or_die(N+1);
	prvs->rho           = double_malloc_or_die(N+1);
	prvs->rho_sum       = double_malloc_or_die(N+1);
	prvs->rho_avg       = double_malloc_or_die(N+1);
	prvs->rho_fit       = double_malloc_or_die(N+1);
	prvs->rho_dev       = double_malloc_or_die(N+1);
	prvs->delta_rho     = double_malloc_or_die(N+1);
	prvs->delta_rho_avg = double_malloc_or_die(N+1);
	prvs->counts        = int_malloc_or_die   (N+1);
	prvs->counts_sum    = double_malloc_or_die(N+1);
	prvs->counts_avg    = double_malloc_or_die(N+1);
	prvs->num_acc       = pparams->termination_num_acc;

	prvs->rho_sum_count = 0;
	for (k = 0; k <= ppoints->N; k++) {
		prvs->rho_sum[k] = 0.0;
		prvs->rho_avg[k] = 0.0;
		prvs->rho_fit[k] = 0.0;
		prvs->rho_dev[k] = 0.0;
		prvs->unitk[k] = (double)k/(double)ppoints->N;
		prvs->counts_sum[k] = 0.0;
		prvs->counts_avg[k] = 0.0;
	}
	prvs->maxjumplen = 0.0;

	init_second_order_cssm_stats(&prvs->H_2cos, num_acc);
	init_first_order_stats      (&prvs->D_1os, num_acc);
	init_first_order_stats      (&prvs->V_1os, num_acc);

	init_second_order_cssm_stats(&prvs->r2_2cos,   num_acc);
	init_second_order_cssm_stats(&prvs->lmax_2cos, num_acc);

	init_second_order_stats(&prvs->jumplenbar_2os,        num_acc);
	init_second_order_stats(&prvs->ell0_2os,              num_acc);
	init_second_order_cssm_stats(&prvs->ellbar_2cos,      num_acc);
	init_second_order_stats(&prvs->cycle_0_spatlen_2os,   num_acc);
	init_second_order_cssm_stats(&prvs->meanspatlen_2cos, num_acc);
	init_second_order_cssm_stats(&prvs->recipmeanspatlen_2cos, num_acc);
	init_second_order_stats(&prvs->ellfbar_2os,           num_acc);
	init_second_order_stats(&prvs->meanspatlenf_2os,      num_acc);
	init_second_order_stats(&prvs->winding_number_2os,    num_acc);
	init_second_order_cssm_stats(&prvs->fS_2cos,          num_acc);
	init_second_order_cssm_stats(&prvs->fW_2cos,          num_acc);

	prvs->ellxy_stats.px = &ppoints->lattice[0][0][0];
	prvs->ellxy_stats.py = &ppoints->lattice[midx][midy][midz];

	init_second_order_pair_stats(&prvs->ellxy_stats.pair_stats, num_acc);

	prvs->jump0_stats.ppoint = &ppoints->lattice[0][0][0];
	init_second_order_stats(&prvs->jump0_stats.sostats, num_acc);

	prvs->kink_stats.ppoint = &ppoints->lattice[0][0][0];
	init_second_order_stats(&prvs->kink_stats.dot_sostats, num_acc);
	init_second_order_stats(&prvs->kink_stats.cos_angle_sostats, num_acc);

	prvs->angle_stats[0].px = &ppoints->lattice[0][0][0];
	prvs->angle_stats[0].py = &ppoints->lattice[0][0][1];

	prvs->angle_stats[1].px = &ppoints->lattice[0][0][0];
	prvs->angle_stats[1].py = &ppoints->lattice[0][1][1];

	prvs->angle_stats[2].px = &ppoints->lattice[0][0][0];
	prvs->angle_stats[2].py = &ppoints->lattice[1][1][1];

	prvs->angle_stats[3].px = &ppoints->lattice[0][0][0];
	prvs->angle_stats[3].py = &ppoints->lattice[midx][midy][midz];

	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		init_second_order_stats(&prvs->angle_stats[k].dot_sostats, num_acc);
		init_second_order_stats(&prvs->angle_stats[k].cos_angle_sostats,
			num_acc);
	}

	ihisto_init(&prvs->Wx_histo, WX_HISTO_MIN, WX_HISTO_MAX);
}

// ----------------------------------------------------------------
// Find per-realization values.  Accumulate sums in a separate routine, below.
void compute_rvs(mcmc_rvs_t* prvs, points_t* ppoints, mcmc_params_t* pparams,
	double H, double D, double V)
{
	int k;
	kink_stats_t* pk = &prvs->kink_stats; // Keystroke saver

	prvs->H = H;
	prvs->D = D;
	prvs->V = V;

	get_rho_L_pi(ppoints, prvs->rho, prvs->counts, &prvs->lmax);
	prvs->r2 = prvs->counts[2];

	prvs->cycle_0_spatlen = get_cycle_0_spatlen(ppoints);
	get_mean_cycle_length(ppoints, &prvs->ellbar, &prvs->meanspatlen,
		&prvs->ellfbar, &prvs->meanspatlenf);
	if (prvs->meanspatlen == 0.0)
		prvs->recipmeanspatlen = 0.0;
	else
		prvs->recipmeanspatlen = 1.0 / prvs->meanspatlen;

	get_pmt_winding_numbers(ppoints, &prvs->Wx, &prvs->Wy, &prvs->Wz);
	prvs->winding_number =
		(double)prvs->Wx*prvs->Wx +
		(double)prvs->Wy*prvs->Wy +
		(double)prvs->Wz*prvs->Wz;
	prvs->fS = prvs->winding_number // f_S = <W^2> L^2 / 3 beta N.
		* (double)pparams->L*(double)pparams->L
		/ 3.0 / pparams->beta / (double)ppoints->N;
	prvs->fW = get_fraction_of_cycles_which_wind(ppoints);

	prvs->ellx = prvs->ellxy_stats.px->pcycinfo->cyclen;
	prvs->elly = prvs->ellxy_stats.py->pcycinfo->cyclen;
	prvs->jumplenbar = get_mean_jump_length(ppoints, &prvs->maxjumplen);

	pk->dot = get_jump_dot(pk->ppoint, pk->ppoint->pfwd, ppoints);
	pk->cos_angle_defined = get_jump_cos_angle(pk->ppoint, pk->ppoint->pfwd,
		ppoints, &pk->cos_angle);

	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		angle_stats_t* pcurr = &prvs->angle_stats[k];
		pcurr->dot = get_jump_dot(pcurr->px, pcurr->py, ppoints);
		pcurr->cos_angle_defined = get_jump_cos_angle(pcurr->px, pcurr->py,
			ppoints, &pcurr->cos_angle);
	}
}

// ----------------------------------------------------------------
void accumulate_rvs(mcmc_rvs_t* prvs, points_t* ppoints,
	mcmc_params_t* pparams, double H, double D, double V)
{
	int k;
	kink_stats_t* pk = &prvs->kink_stats; // Keystroke saver

	compute_rvs(prvs, ppoints, pparams, H, D, V);

	prvs->rho_sum_count++;
	vector_accumulate(prvs->rho_sum, prvs->rho, ppoints->N+1);
	vector_iaccumulate(prvs->counts_sum, prvs->counts, ppoints->N+1);

	if (pk->cos_angle_defined)
		accum_second_order_stats(&pk->cos_angle_sostats, pk->cos_angle);

	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		angle_stats_t* pcurr = &prvs->angle_stats[k];
		accum_second_order_stats(&pcurr->dot_sostats, pcurr->dot);
		if (pcurr->cos_angle_defined)
			accum_second_order_stats(&pcurr->cos_angle_sostats,
				pcurr->cos_angle);
	}

	// Winding-number histogram.
	ihisto_incr(&prvs->Wx_histo, prvs->Wx);

	accum_second_order_cssm_stats(&prvs->H_2cos, H);
	accum_first_order_stats      (&prvs->D_1os,  D);
	accum_first_order_stats      (&prvs->V_1os,  V);

	accum_second_order_cssm_stats(&prvs->r2_2cos, (double)prvs->r2);
	accum_second_order_cssm_stats(&prvs->lmax_2cos, (double)prvs->lmax);

	accum_second_order_stats(&prvs->jumplenbar_2os, prvs->jumplenbar);
	accum_second_order_stats(&prvs->ell0_2os, prvs->ellx);
	accum_second_order_cssm_stats(&prvs->ellbar_2cos, prvs->ellbar);
	accum_second_order_stats(&prvs->cycle_0_spatlen_2os, prvs->cycle_0_spatlen);
	accum_second_order_cssm_stats(&prvs->meanspatlen_2cos, prvs->meanspatlen);
	accum_second_order_cssm_stats(&prvs->recipmeanspatlen_2cos,
		prvs->recipmeanspatlen);
	accum_second_order_stats(&prvs->ellfbar_2os, prvs->ellfbar);
	accum_second_order_stats(&prvs->meanspatlenf_2os, prvs->meanspatlenf);
	accum_second_order_stats(&prvs->winding_number_2os, prvs->winding_number);
	accum_second_order_cssm_stats(&prvs->fS_2cos, prvs->fS);
	accum_second_order_cssm_stats(&prvs->fW_2cos, prvs->fW);

	accum_second_order_pair_stats(&prvs->ellxy_stats.pair_stats,
		(double)prvs->ellx, (double)prvs->elly);

	accum_second_order_stats(&prvs->jump0_stats.sostats,
		prvs->jump0_stats.ppoint->fwd_d);

	accum_second_order_stats(&pk->dot_sostats, pk->dot);
}

// ----------------------------------------------------------------
void display_realization_rvs(double H, double D, double V,
	mcmc_rvs_t* prvs, points_t* ppoints, rv_verbosity_t* pvb)
{
	int k;

	if (pvb->pmt)
		pmt_print_images(ppoints);
	if (pvb->cycdec == 1)
		pmt_print_cyclens(ppoints);
	if (pvb->cycdec == 2)
		pmt_print_cycdec(ppoints);

	// xxx opt to constrain by len: only print for ello to elhi.
	if (pvb->counts >= 2) {
		for (k = 1; k <= ppoints->N; k++) {
			printf(" %d", prvs->counts[k]);
		}
		printf("\n");
	}

	// H, D, & V are printed from within the calling routine.
	if (pvb->r2 >= 1) {
		printf("%d # r2\n", prvs->r2);
	}

	if (pvb->lmax)
		printf("%d # lmax\n", prvs->lmax);
	if (pvb->ell)
		printf("%d %d # ellx elly\n", prvs->ellx, prvs->elly);
	if (pvb->ellbar) {
		printf("%11.7lf # ellbar\n",   prvs->ellbar);
		printf("%11.7lf # ellfbar\n",  prvs->ellfbar);
	}
	if (pvb->cycle_0_spatlen)
		printf("%.7lf # cyc0 spatlen\n", prvs->cycle_0_spatlen);
	if (pvb->meanspatlen) {
		printf("%.7lf # meanspatlen\n",  prvs->meanspatlen);
		printf("%.7lf # recipmeanspatlen\n",  prvs->recipmeanspatlen);
		printf("%.7lf # meanspatlenf\n", prvs->meanspatlenf);
	}

	if (pvb->winding_cycle) {
		if (prvs->winding_number != 0.0) {
			printf("\n%s", SMALL_SEP);
			printf("# Winding cycles (accumulation phase):\n");
			print_winding_cycles(ppoints);
			printf("\n%s", SMALL_SEP);
			printf("\n");
		}
	}
	if (pvb->winding_number== 1) {
		printf("%.7lf # Winding number\n", prvs->winding_number);
	}
	if (pvb->winding_number >= 2) {
		printf("%6d %6d %6d # Wx Wy Wz\n",
			prvs->Wx, prvs->Wy, prvs->Wz);
	}

	if (pvb->jump_length == 1) {
		printf("%11.7lf # fwd_d\n", prvs->jump0_stats.ppoint->fwd_d);
	}
	if (pvb->jump_target == 1) {
		printf("(%d,%d,%d)%d # jump target\n",
			prvs->jump0_stats.ppoint->pfwd->c.selfi,
			prvs->jump0_stats.ppoint->pfwd->c.selfj,
			prvs->jump0_stats.ppoint->pfwd->c.selfk,
			prvs->jump0_stats.ppoint->fwd_dsq);
	}
	if (pvb->jump_vector == 1) {
		point_t* x = &ppoints->lattice[0][0][0];
		point_t* y = x->pfwd;
		coord_t diff = get_coord_diff(&x->c, &y->c, &ppoints->dims);
		printf("%3d %3d %3d # jump vector\n",
			diff.selfi, diff.selfj, diff.selfk);
	}

	if (pvb->kink_dot) {
		printf("%11.7lf # kink dot\n", prvs->kink_stats.dot);
	}
	if (pvb->kink_cos_angle) {
		if (prvs->kink_stats.cos_angle_defined)
			printf("%11.7lf # kink cos(angle) \n", prvs->kink_stats.cos_angle);
		else
			printf("%11.7lf # kink cos(angle) \n", -2.0);
	}

	if (pvb->cos_angle) {
		for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
			angle_stats_t* pcurr = &prvs->angle_stats[k];
			if (pcurr->cos_angle_defined)
				printf("%11.7lf ", pcurr->cos_angle);
			else
				printf("%11.7lf ", -2.0);
		}
		printf("# cos_angle\n");
	}
	if (pvb->dot) {
		for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
			angle_stats_t* pcurr = &prvs->angle_stats[k];
			printf("%11.7lf ", pcurr->dot);
		}
		printf("# dot\n");
	}

	fflush(stdout);
}

// ----------------------------------------------------------------
void compute_rv_means(mcmc_rvs_t* prvs, points_t* ppoints)
{
	double recip_count = 1.0 / prvs->rho_sum_count;
	int k;
	kink_stats_t* pk = &prvs->kink_stats; // Keystroke saver

	vector_scale(prvs->rho_avg, prvs->rho_sum, recip_count, ppoints->N+1);
	vector_scale(prvs->counts_avg, prvs->counts_sum, recip_count, ppoints->N+1);

	// Compute first differences of final rho.
	compute_vector_delta(prvs->delta_rho,     prvs->rho,     ppoints->N+1);
	compute_vector_delta(prvs->delta_rho_avg, prvs->rho_avg, ppoints->N+1);

	// Compute rho_infty(T).
	prvs->rho_infty_T = compute_rho_infty(prvs->rho_fit, prvs->rho_avg,
		prvs->rho, ppoints->N, 2);
	for (k = 0; k <= ppoints->N; k++)
		prvs->rho_dev[k] = prvs->rho_avg[k] - (double)k/(double)ppoints->N;

	finish_second_order_cssm_stats(&prvs->H_2cos);
	finish_first_order_stats      (&prvs->D_1os);
	finish_first_order_stats      (&prvs->V_1os);

	finish_second_order_cssm_stats(&prvs->r2_2cos);
	finish_second_order_cssm_stats(&prvs->lmax_2cos);

	finish_second_order_stats(&prvs->jumplenbar_2os);
	finish_second_order_stats(&prvs->ell0_2os);
	finish_second_order_cssm_stats(&prvs->ellbar_2cos);
	finish_second_order_stats(&prvs->cycle_0_spatlen_2os);
	finish_second_order_cssm_stats(&prvs->meanspatlen_2cos);
	finish_second_order_cssm_stats(&prvs->recipmeanspatlen_2cos);
	finish_second_order_stats(&prvs->ellfbar_2os);
	finish_second_order_stats(&prvs->meanspatlenf_2os);
	finish_second_order_stats(&prvs->winding_number_2os);
	finish_second_order_cssm_stats(&prvs->fS_2cos);
	finish_second_order_cssm_stats(&prvs->fW_2cos);

	finish_second_order_pair_stats(&prvs->ellxy_stats.pair_stats);

	finish_second_order_stats(&prvs->jump0_stats.sostats);

	finish_second_order_stats(&pk->dot_sostats);
	if (pk->cos_angle_defined)
		finish_second_order_stats(&pk->cos_angle_sostats);

	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		angle_stats_t* pcurr = &prvs->angle_stats[k];
		finish_second_order_stats(&pcurr->dot_sostats);
		finish_second_order_stats(&pcurr->cos_angle_sostats);
	}
}

// ----------------------------------------------------------------
void display_rv_means(mcmc_rvs_t* prvs, mcmc_params_t* pparams,
	points_t* ppoints, rv_verbosity_t* pvb)
{
	int k;
	const int max_count = 10;
	int cutoff_count = ppoints->N < max_count ? ppoints->N : max_count;
	kink_stats_t* pk = &prvs->kink_stats; // Keystroke saver

	if (pvb->rho == 1) {
		printf("#k/N rho[0,k] <rho[0,k]> <rho[0,k]>-k/N rhofit "
			"rho[k,k] <rho[k,k]>\n");
		printf("#--- -------- ---------- -------------- ------ "
			"-------- ----------\n");
		vector_col7_printf_decimated(
			prvs->unitk,
			prvs->rho,
			prvs->rho_avg,
			prvs->rho_dev,
			prvs->rho_fit,
			prvs->delta_rho,
			prvs->delta_rho_avg,
			ppoints->N+1, "%.6lf");
	}

	print_name_and_dvalue("fI", prvs->rho_infty_T);
	printf("#\n");

	if (pvb->counts == 1) {
		printf("#k <counts[k]>\n");
		for (k = 1; k <= ppoints->N; k++) {
			printf("%4d %11.7lf\n", k, prvs->counts_avg[k]);
		}
	}
	else {
		double sum = 0.0;
		for (k = 1; k <= cutoff_count; k++)
			printf("# k= %2d <counts> %11.7lf\n", k, prvs->counts_avg[k]);
		for (k = 1; k <= ppoints->N; k++)
			sum += k * prvs->counts_avg[k];
		print_name_and_dvalue("k*counts sum", sum);
	}
	printf("#\n");

	// xxx per RV ... tauexp = tauint_to_tauexp(prvs->tauint);
	// xxx per RV ... eta    = tauint_to_eta   (prvs->tauint);

	// xxx per RV ... print_name_and_dvalue("tauint", prvs->tauint);
	// xxx per RV ... print_name_and_dvalue("tauexp", tauexp);
	// xxx per RV ... print_name_and_dvalue("eta",    eta);
	// xxx per RV ... printf("#\n");

	printf("# GRU quantities:\n");
	print_name_and_dvalue("alpha GRU", 0.25 / pparams->beta);
	print_name_and_dvalue("<lmax>/|X|",
		(double)prvs->lmax_2cos.avg/(double)ppoints->N);
	print_name_and_dvalue("phi(alpha)", prvs->rho_infty_T);
	print_name_and_dvalue("macroscopic-cycle quotient",
		(double)prvs->lmax_2cos.avg / (double)ppoints->N / prvs->rho_infty_T);
	printf("#\n");

	print_second_order_cssm_stats(&prvs->H_2cos, "H");
	printf("#\n");

	print_first_order_stats      (&prvs->D_1os, "D");
	print_first_order_stats      (&prvs->V_1os, "V");
	print_name_and_dvalue("mean_h", prvs->H_2cos.avg / ppoints->N);
	print_name_and_dvalue("mean_d", prvs->D_1os.avg / ppoints->N);
	print_name_and_dvalue("mean_v", prvs->V_1os.avg / ppoints->N);
	printf("#\n");

	print_second_order_cssm_stats(&prvs->r2_2cos,        "r2");
	printf("#\n");
	print_second_order_cssm_stats(&prvs->lmax_2cos,      "lmax");
	printf("#\n");

	print_second_order_stats(&prvs->jumplenbar_2os,      "jumplenbar");
	print_name_and_dvalue("maxjumplen", prvs->maxjumplen);
	printf("#\n");
	print_second_order_stats(&prvs->ell0_2os,            "ell0");
	printf("#\n");
	print_second_order_cssm_stats(&prvs->ellbar_2cos,    "ellbar");
	printf("#\n");
	print_second_order_stats(&prvs->ellfbar_2os,         "ellfbar");
	printf("#\n");
	print_second_order_stats(&prvs->cycle_0_spatlen_2os, "cyc0spatlen");
	printf("#\n");
	print_second_order_cssm_stats(&prvs->meanspatlen_2cos, "meanspatlen");
	printf("#\n");
	print_second_order_cssm_stats(&prvs->recipmeanspatlen_2cos,
		"recipmeanspatlen");
	printf("#\n");
	print_second_order_stats(&prvs->meanspatlenf_2os,    "meanspatlenf");
	printf("#\n");

	print_name_and_dvalue("<ell0>/L",
		prvs->ell0_2os.avg / ppoints->dims.nx);
	print_name_and_dvalue("<ellbar>/L",
		prvs->ellbar_2cos.avg / ppoints->dims.nx);
	print_name_and_dvalue("<ellfbar>/L",
		prvs->ellfbar_2os.avg / ppoints->dims.nx);
	print_name_and_dvalue("<cyc0spatlen>/L",
		prvs->cycle_0_spatlen_2os.avg / ppoints->dims.nx);
	print_name_and_dvalue("<meanspatlen>/L",
		prvs->meanspatlen_2cos.avg / ppoints->dims.nx);
	print_name_and_dvalue("<meanspatlenf>/L",
		prvs->meanspatlenf_2os.avg / ppoints->dims.nx);
	printf("#\n");

	print_name_and_dvalue("1/<ell0>",
		1.0 / prvs->ell0_2os.avg);
	print_name_and_dvalue("1/<ellbar>",
		1.0 / prvs->ellbar_2cos.avg);
	print_name_and_dvalue("1/<cyc0spatlen>",
		1.0 / prvs->cycle_0_spatlen_2os.avg);
	printf("#\n");

	print_second_order_stats(&prvs->winding_number_2os, "wno");
	printf("#\n");
	print_second_order_cssm_stats(&prvs->fS_2cos,       "fS");
	print_name_and_dvalue("recip_fS",
		(prvs->fS_2cos.avg == 0.0) ? 0.0 : 1.0 / prvs->fS_2cos.avg);
	printf("#\n");
	print_second_order_cssm_stats(&prvs->fW_2cos,       "fW");
	print_name_and_dvalue("recip_fW",
		(prvs->fW_2cos.avg == 0.0) ? 0.0 : 1.0 / prvs->fW_2cos.avg);
	printf("#\n");

	print_second_order_pair_stats(&prvs->ellxy_stats.pair_stats,
		"ellx", "elly");
	printf("#\n");

	print_second_order_stats(&prvs->jump0_stats.sostats, "jump0len");
	printf("#\n");

	print_second_order_stats(&pk->dot_sostats, "kinkdot");
	print_second_order_stats(&pk->cos_angle_sostats, "kinkcosangle");
	print_name_and_ivalue("kinkcosangle num defined",
		pk->cos_angle_sostats.time_series.nused);
	printf("#\n");

	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		angle_stats_t* pcurr = &prvs->angle_stats[k];

		printf("# * Angle stats for sites %d,%d,%d and %d,%d,%d:\n",
			pcurr->px->c.selfi, pcurr->px->c.selfj, pcurr->px->c.selfk,
			pcurr->py->c.selfi, pcurr->py->c.selfj, pcurr->py->c.selfk);

		print_second_order_stats(&pcurr->dot_sostats, "dot");
		print_second_order_stats(&pcurr->cos_angle_sostats, "cos_angle");
		print_name_and_ivalue("num defined",
			pcurr->cos_angle_sostats.time_series.nused);
	}
	printf("#\n");

	ihisto_report(&prvs->Wx_histo, "Wx");
	printf("#\n");
}

// ----------------------------------------------------------------
void free_rvs(mcmc_rvs_t* prvs)
{
	int k;

	free(prvs->unitk);
	free(prvs->rho);
	free(prvs->rho_sum);
	free(prvs->rho_avg);
	free(prvs->rho_fit);
	free(prvs->rho_dev);
	free(prvs->delta_rho);
	free(prvs->delta_rho_avg);
	free(prvs->counts);
	free(prvs->counts_sum);
	free(prvs->counts_avg);

	prvs->unitk         = 0;
	prvs->rho           = 0;
	prvs->rho_sum       = 0;
	prvs->rho_avg       = 0;
	prvs->rho_fit       = 0;
	prvs->rho_dev       = 0;
	prvs->delta_rho     = 0;
	prvs->delta_rho_avg = 0;
	prvs->counts        = 0;
	prvs->counts_sum    = 0;
	prvs->counts_avg    = 0;

	free_second_order_cssm_stats(&prvs->H_2cos);
	free_first_order_stats      (&prvs->D_1os);
	free_first_order_stats      (&prvs->V_1os);
	free_second_order_cssm_stats(&prvs->r2_2cos);
	free_second_order_cssm_stats(&prvs->lmax_2cos);
	free_second_order_stats(&prvs->jumplenbar_2os);
	free_second_order_stats(&prvs->ell0_2os);
	free_second_order_cssm_stats(&prvs->ellbar_2cos);
	free_second_order_stats(&prvs->cycle_0_spatlen_2os);
	free_second_order_cssm_stats(&prvs->meanspatlen_2cos);
	free_second_order_cssm_stats(&prvs->recipmeanspatlen_2cos);
	free_second_order_stats(&prvs->ellfbar_2os);
	free_second_order_stats(&prvs->meanspatlenf_2os);
	free_second_order_stats(&prvs->winding_number_2os);
	free_second_order_cssm_stats(&prvs->fS_2cos);
	free_second_order_cssm_stats(&prvs->fW_2cos);
	free_second_order_pair_stats(&prvs->ellxy_stats.pair_stats);
	free_second_order_stats(&prvs->jump0_stats.sostats);
	free_second_order_stats(&prvs->kink_stats.dot_sostats);
	free_second_order_stats(&prvs->kink_stats.cos_angle_sostats);
	for (k = 0; k < NUM_ANGLE_STATS_SITES; k++) {
		free_second_order_stats(&prvs->angle_stats[k].dot_sostats);
		free_second_order_stats(&prvs->angle_stats[k].cos_angle_sostats);
	}

	ihisto_free(&prvs->Wx_histo);
}
