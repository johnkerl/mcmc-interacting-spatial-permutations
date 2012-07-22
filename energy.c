// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "points.h"
#include "pmt.h"
#include "energy.h"
#include "interactions.h"
#include "util.h"
#include "worm_moves.h"
#include "checks.h"

// ----------------------------------------------------------------
double get_Delta_V_SO(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int interaction_type,
	int* pdxy, int* pdyx, int* px_y_same_cycle)
{
	if (interaction_type != RELL_INTERACTIONS) {
		// get_Delta_rell computes dxy, dyx, and x_y_same_cycle as a side
		// effect.  The others do not.
		find_dxy_dyx_xoy(ppoints, x, y, pdxy, pdyx, px_y_same_cycle);
	}

	if (interaction_type == R2_INTERACTIONS) {
		return alphas[2] * get_Delta_r2(ppoints, x, pix, y, piy);
	}
	else if (interaction_type == RELL_INTERACTIONS) {
		return get_Delta_rell(ppoints, x, pix, y, piy, alphas,
			pdxy, pdyx, px_y_same_cycle);
	}
	else if (interaction_type == VZWL_INTERACTIONS) {
		return get_Delta_Vzw_SO(ppoints, x, pix, y, piy, beta, alphas,
			interaction_type);
	}
	else if (interaction_type == VZWU_INTERACTIONS) {
		return get_Delta_Vzw_SO(ppoints, x, pix, y, piy, beta, alphas,
			interaction_type);
	}
	else {
		return 0.0;
	}
}

// ----------------------------------------------------------------
// Assumptions about the order of arguments:
// * Open:      get_Delta_V_worm(ppoints, x,    pix, w,    w,   ...);
// * Head swap: get_Delta_V_worm(ppoints, x,    pix, piiw, w,,  ...);
// * Tail swap: get_Delta_V_worm(ppoints, x,    pix, w,    piw, ...);
// * Close:     get_Delta_V_worm(ppoints, piiw, w,   w,    piw, ...);
double get_Delta_V_worm(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	int worm_move, double beta, double* alphas, int interaction_type,
	int* pdxy, int* pdyx, int* px_y_same_cycle)
{
	if (interaction_type != RELL_INTERACTIONS) {
		// get_Delta_rell computes dxy, dyx, and x_y_same_cycle as a side
		// effect.  The others do not.
		find_dxy_dyx_xoy(ppoints, x, y, pdxy, pdyx, px_y_same_cycle);
	}

	if (interaction_type == R2_INTERACTIONS) {
		return alphas[2] * get_Delta_r2(ppoints, x, pix, y, piy);
	}
	else if (interaction_type == RELL_INTERACTIONS) {
		return get_Delta_rell(ppoints, x, pix, y, piy, alphas,
			pdxy, pdyx, px_y_same_cycle);
	}
	else if ((interaction_type == VZWL_INTERACTIONS) ||
		(interaction_type == VZWU_INTERACTIONS))
	{
		switch (worm_move) {
		case WORM_OPEN:
			return get_Delta_Vzw_worm_open     (ppoints, x, pix, y, piy,
				beta, alphas, interaction_type);
			break;
		case WORM_CLOSE:
			return get_Delta_Vzw_worm_close    (ppoints, x, pix, y, piy,
				beta, alphas, interaction_type);
			break;
		case WORM_HEAD_SWAP:
			return get_Delta_Vzw_worm_head_swap(ppoints, x, pix, y, piy,
				beta, alphas, interaction_type);
			break;
		case WORM_TAIL_SWAP:
			return get_Delta_Vzw_worm_tail_swap(ppoints, x, pix, y, piy,
				beta, alphas, interaction_type);
			break;
		default:
			fprintf(stderr, "Coding error detected: file %s line %d.\n",
				__FILE__, __LINE__);
			exit(1);
			break;
		}
	}
	else {
		return 0.0;
	}
}

// ----------------------------------------------------------------
// See the dissertation for an explanation.

int get_Delta_r2(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy)
{
	int Delta_r2 = 0;
	point_t* pi2x = pix->pfwd;
	point_t* pi2y = piy->pfwd;

	if (piy == pix) // Case 0
		return 0;

	if (x == pix) { // Case 1
		if      (y == piy)  Delta_r2 =  1; // Case 1a.
		else if (y == pi2y) Delta_r2 = -1; // Case 1b.
	}

	else if (y == piy) { // Case 2
		if      (x == pix)  Delta_r2 =  1; // Case 2a; same as 1a.
		else if (x == pi2x) Delta_r2 = -1; // Case 2b.
	}

	else if (x == piy) { // Case 3
		if      (pix == y)  Delta_r2 = -1; // Case 3a.
		else if (pi2x == y) Delta_r2 =  1; // Case 3b.
	}

	else if (pix == y) { // Case 4
		if      (x == piy)  Delta_r2 =  -1; // Case 4a; same as 3a.
		else if (pi2y == x) Delta_r2 =   1;
	}

	else if (pi2x == x) { // Case 5
		if (pi2y == y)      Delta_r2 = -2; // Case 5a.
		else                Delta_r2 = -1; // Case 5b.
	}

	else if (pi2y == y) { // Case 6
		if (pi2x == x)      Delta_r2 = -2; // Case 6a; same as 5a.
		else                Delta_r2 = -1; // Case 6b.
	}

	else if (pi2x == y) { // Case 7
		if (pi2y == x)      Delta_r2 =  2; // Case 7a.
		else                Delta_r2 =  1; // Case 7b.
	}

	else if (pi2y == x) { // Case 8
		if (pi2x == y)      Delta_r2 =  2; // Case 8a.
		else                Delta_r2 =  1; // Case 8b.
	}

	return Delta_r2;
}

// ----------------------------------------------------------------
double get_Delta_rell(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy, double* alphas,
	int* pdxy, int* pdyx, int* px_y_same_cycle)
{
	double Delta_rell = 0.0;

	// In the same cycle, but neighbors.
	if (pix == piy) {
		// No change.
		*pdxy = 0;
		*pdyx = 0;
		*px_y_same_cycle = 1;
	}

	// In different cycles, but one of x and y is in a one-cycle.
	else if (x == pix) {
		Delta_rell -= alphas[1];
		Delta_rell -= alphas[y->pcycinfo->cyclen];
		Delta_rell += alphas[y->pcycinfo->cyclen+1];
		*pdxy = -1;
		*pdyx = -1;
		*px_y_same_cycle = 0;
	}

	else if (y == piy) {
		Delta_rell -= alphas[1];
		Delta_rell -= alphas[x->pcycinfo->cyclen];
		Delta_rell += alphas[x->pcycinfo->cyclen+1];
		*pdxy = -1;
		*pdyx = -1;
		*px_y_same_cycle = 0;
	}

	// In the same cycle, but neighbors.
	else if (y == pix) {
		Delta_rell -= alphas[x->pcycinfo->cyclen];
		Delta_rell += alphas[x->pcycinfo->cyclen-1];
		Delta_rell += alphas[1];
		*pdxy = 1;
		*pdyx = x->pcycinfo->cyclen - 1;
		*px_y_same_cycle = 1;
	}
	else if (x == piy) {
		Delta_rell -= alphas[y->pcycinfo->cyclen];
		Delta_rell += alphas[y->pcycinfo->cyclen-1];
		Delta_rell += alphas[1];
		*pdxy = y->pcycinfo->cyclen - 1;
		*pdyx = 1;
		*px_y_same_cycle = 1;
	}

	else {
		find_dxy_dyx_xoy(ppoints, x, y, pdxy, pdyx, px_y_same_cycle);
		if (*px_y_same_cycle) { // It splits into separate cycles.
			Delta_rell -= alphas[*pdxy + *pdyx];
			Delta_rell += alphas[*pdxy];
			Delta_rell += alphas[*pdyx];
		}
		else { // They merge into a common cycle.
			Delta_rell -= alphas[x->pcycinfo->cyclen];
			Delta_rell -= alphas[y->pcycinfo->cyclen];
			Delta_rell += alphas[x->pcycinfo->cyclen + y->pcycinfo->cyclen];
		}
	}
	return Delta_rell;
}

// ----------------------------------------------------------------
// a/|w||z| (|w| + |z|) exp(-|z+w|^2 / 8 beta)
// w = xi-xj = start diff
// z = xpii-xpi = end diff
double inline get_Vzw_lower_of_jump_pair(points_t* ppoints,
	point_t* u,  point_t* piu, point_t* v, point_t* piv,
	double beta, double* alphas)
{
	double a = alphas[2];
	coord_t w       = get_coord_diff(&v->c,   &u->c,   &ppoints->dims);
	coord_t z       = get_coord_diff(&piv->c, &piu->c, &ppoints->dims);
	coord_t zpw     = get_coord_sum (&w, &z,  &ppoints->dims);
	double normw2   = get_normsq(&w);
	double normz2   = get_normsq(&z);
	double normw    = sqrt(normw2);
	double normz    = sqrt(normz2);
	double normzpw2 = get_normsq(&zpw);
	return a/normw/normz * (normw+normz) * exp(-normzpw2/8/beta);
}

// ----------------------------------------------------------------
// a/|w||z|  *
// exp(-z.w / 4 beta) *
// (|w| exp(-|z|^2/4beta) + |z| exp(-|w|^2/4beta))
double inline get_Vzw_upper_of_jump_pair(points_t* ppoints,
	point_t* u,  point_t* piu, point_t* v, point_t* piv,
	double beta, double* alphas)
{
	double a = alphas[2];
	coord_t w     = get_coord_diff(&v->c,   &u->c,   &ppoints->dims);
	coord_t z     = get_coord_diff(&piv->c, &piu->c, &ppoints->dims);
	double normw2 = get_normsq(&w);
	double normz2 = get_normsq(&z);
	double normw  = sqrt(normw2);
	double normz  = sqrt(normz2);
	double zdotw  = get_dot(&z, &w);
	return a/normw/normz *
		exp(-zdotw/4/beta) *
		(normw * exp(-normz2/4/beta) + normz * exp(-normw2/4/beta));
}

// ----------------------------------------------------------------
double get_Vzw_of_jump_pair(points_t* ppoints,
	point_t* x,  point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int interaction_type)
{
	if (interaction_type == VZWL_INTERACTIONS)
		return get_Vzw_lower_of_jump_pair(ppoints, x, pix, y, piy,
			beta, alphas);
	else if (interaction_type == VZWU_INTERACTIONS)
		return get_Vzw_upper_of_jump_pair(ppoints, x, pix, y, piy,
			beta, alphas);
	else {
		fprintf(stderr,
			"Coding error detected in file %s line %d.\n",
			__FILE__, __LINE__);
		exit(1);
	}
}

// ----------------------------------------------------------------
// See my dissertation for explanation.
//
// SO:
// + sum { v!=x,y} [ V(x, piy, v, piv) - V(x, pix,  v, piv) ]
// + sum { v!=x,y} [ V(y, pix, v, piv) - V(y, piy,  v, piv) ]
// +                 V(x, piy, y, pix) - V(x, pix,  y, piy)

double get_Delta_Vzw_SO(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type)
{
	double Delta_V = 0.0;
	int nx = ppts->dims.nx;
	int ny = ppts->dims.ny;
	int nz = ppts->dims.nz;
	int vi, vj, vk;
	point_t* v, * piv;

	if (alphas[2] == 0.0)
		return 0.0;

	for (vi = 0; vi < nx; vi++) {
		for (vj = 0; vj < ny; vj++) {
			for (vk = 0; vk < nz; vk++) {
				v = &ppts->lattice[vi][vj][vk];
				piv = v->pfwd;
				if (v == x)
					continue;
				if (v == y)
					continue;

				Delta_V += get_Vzw_of_jump_pair(
					ppts, x, piy, v, piv, beta, alphas, type);
				Delta_V += get_Vzw_of_jump_pair(
					ppts, y, pix, v, piv, beta, alphas, type);
				Delta_V -= get_Vzw_of_jump_pair(
					ppts, x, pix, v, piv, beta, alphas, type);
				Delta_V -= get_Vzw_of_jump_pair(
					ppts, y, piy, v, piv, beta, alphas, type);
			}
		}
	}

	if (x != y) {

		Delta_V += get_Vzw_of_jump_pair(ppts,x,piy,y,pix, beta,alphas,type);
		Delta_V -= get_Vzw_of_jump_pair(ppts,x,pix,y,piy, beta,alphas,type);
	}

	return Delta_V;
}

// ----------------------------------------------------------------
// See my dissertation for explanation.

double get_Delta_Vzw_worm_open(points_t* ppts,
	point_t* x, point_t* pix, point_t* w, point_t* piw_is_w,
	double beta, double* alphas, int type)
{
	double Delta_V = 0.0;
	int nx = ppts->dims.nx;
	int ny = ppts->dims.ny;
	int nz = ppts->dims.nz;
	int vi, vj, vk;
	point_t* v, * piv;

	if (alphas[2] == 0.0)
		return 0.0;

	for (vi = 0; vi < nx; vi++) {
		for (vj = 0; vj < ny; vj++) {
			for (vk = 0; vk < nz; vk++) {
				v = &ppts->lattice[vi][vj][vk];
				piv = v->pfwd;
				if (v == x)
					continue;
				// v can't equal w since w is non-spatial (not in the lattice)
				// and we are looping over spatial v (in the lattice).

				Delta_V -= get_Vzw_of_jump_pair(
					ppts, x, pix, v, piv, beta, alphas, type);
			}
		}
	}

	return Delta_V;
}

// ----------------------------------------------------------------
double get_Delta_Vzw_worm_close(points_t* ppts,
	point_t* piiw, point_t* w, point_t* also_w, point_t* piw,
	double beta, double* alphas, int type)
{
	double Delta_V = 0.0;
	int nx = ppts->dims.nx;
	int ny = ppts->dims.ny;
	int nz = ppts->dims.nz;
	int vi, vj, vk;
	point_t* v, * piv;

	if (alphas[2] == 0.0)
		return 0.0;

	for (vi = 0; vi < nx; vi++) {
		for (vj = 0; vj < ny; vj++) {
			for (vk = 0; vk < nz; vk++) {
				v = &ppts->lattice[vi][vj][vk];
				piv = v->pfwd;
				if (v == piiw)
					continue;
				// v can't equal w since w is non-spatial (not in the lattice)
				// and we are looping over spatial v (in the lattice).

				Delta_V += get_Vzw_of_jump_pair(
					ppts, piiw, piw, v, piv, beta, alphas, type);
			}
		}
	}

	return Delta_V;
}

// ----------------------------------------------------------------
double get_Delta_Vzw_worm_head_swap(points_t* ppts,
	point_t* x, point_t* pix, point_t* piiw, point_t* w,
	double beta, double* alphas, int type)
{
	double Delta_V = 0.0;
	int nx = ppts->dims.nx;
	int ny = ppts->dims.ny;
	int nz = ppts->dims.nz;
	int vi, vj, vk;
	point_t* v, * piv;

	if (alphas[2] == 0.0)
		return 0.0;

	for (vi = 0; vi < nx; vi++) {
		for (vj = 0; vj < ny; vj++) {
			for (vk = 0; vk < nz; vk++) {
				v = &ppts->lattice[vi][vj][vk];
				piv = v->pfwd;
				if (v == x)
					continue;
				if (v == piiw)
					continue;

				Delta_V += get_Vzw_of_jump_pair(
					ppts, piiw, pix, v, piv, beta, alphas, type);
				Delta_V -= get_Vzw_of_jump_pair(
					ppts, x,    pix, v, piv, beta, alphas, type);
			}
		}
	}

	return Delta_V;
}

// ----------------------------------------------------------------
double get_Delta_Vzw_worm_tail_swap(points_t* ppts,
	point_t* x, point_t* pix, point_t* w, point_t* piw,
	double beta, double* alphas, int type)
{
	double Delta_V = 0.0;
	int nx = ppts->dims.nx;
	int ny = ppts->dims.ny;
	int nz = ppts->dims.nz;
	int vi, vj, vk;
	point_t* v, * piv;
	point_t* piiw = w->pbwd;

	if (alphas[2] == 0.0)
		return 0.0;

	for (vi = 0; vi < nx; vi++) {
		for (vj = 0; vj < ny; vj++) {
			for (vk = 0; vk < nz; vk++) {
				v = &ppts->lattice[vi][vj][vk];
				piv = v->pfwd;
				if (v == x)
					continue;
				if (v == piiw)
					continue;
				// v can't equal w since w is non-spatial (not in the lattice)
				// and we are looping over spatial v (in the lattice).

				Delta_V += get_Vzw_of_jump_pair(
					ppts, x, piw, v, piv, beta, alphas, type);
				Delta_V -= get_Vzw_of_jump_pair(
					ppts, x, pix, v, piv, beta, alphas, type);
			}
		}
	}

	return Delta_V;
}

// ----------------------------------------------------------------
double get_Vzw_of_pi(points_t* ppoints, double beta, double* alphas,
	int* counts, int interaction_type)
{
	double Vzw = 0.0;
	int i, j, k, ell;
	int N = ppoints->N;
	point_t** list;
	point_t* x, * y;

	if (alphas[2] == 0.0)
		return 0.0;

	list = (point_t**)malloc_or_die(N * sizeof(point_t*));

	ell = 0;
	for (i = 0; i < ppoints->dims.nx; i++)
		for (j = 0; j < ppoints->dims.ny; j++)
			for (k = 0; k < ppoints->dims.nz; k++)
				list[ell++] = &ppoints->lattice[i][j][k];

	for (i = 0; i < N; i++) {
		x = list[i];
		for (j = i+1; j < N; j++) {
			y = list[j];
			Vzw += get_Vzw_of_jump_pair(ppoints, x, x->pfwd, y, y->pfwd,
				beta, alphas, interaction_type);
		}
	}

	free(list);
	return Vzw;
}

// ----------------------------------------------------------------
// For initial setup.
// This should be called only for closed permutations.
double get_H_of_pi(points_t* ppoints, double beta, double* alphas,
	int interaction_type, double* pD, double* pV)
{
	int i, j, k;
	double D = 0.0;
	double V = 0.0;
	int r2;
	int* counts = int_malloc_or_die(ppoints->N + 1);

	// Distance (single-jump) terms.
	for (i = 0; i < ppoints->dims.nx; i++)
		for (j = 0; j < ppoints->dims.ny; j++)
			for (k = 0; k < ppoints->dims.nz; k++)
				D += ppoints->lattice[i][j][k].fwd_dsq;
	D *= 0.25 / beta;

	// Interaction (jump-pair) terms.
	r2 = pmt_get_num_two_cycles(ppoints);
	if (interaction_type == NO_INTERACTIONS) {
	}
	else if (interaction_type == R2_INTERACTIONS) {
		V += alphas[2] * r2;
	}
	else if (interaction_type == RELL_INTERACTIONS) {
		int ell;
		int N = ppoints->N;
		int lmax_unused_here;
		pmt_get_cycle_counts_and_lmax(ppoints, counts, &lmax_unused_here);
		for (ell = 1; ell <= N; ell++)
			V += alphas[ell] * counts[ell];
	}
	else if (interaction_type == VZWL_INTERACTIONS) {
		V += get_Vzw_of_pi(ppoints, beta, alphas, counts, interaction_type);
	}
	else if (interaction_type == VZWU_INTERACTIONS) {
		V += get_Vzw_of_pi(ppoints, beta, alphas, counts, interaction_type);
	}

	if (pD) *pD = D;
	if (pV) *pV = V;

	free(counts);
	return D+V;
}

// ----------------------------------------------------------------
// For debugging.
// This should be called only for closed permutations.
void check_H(points_t* ppoints,
	double beta, double* alphas, int interaction_type,
	double H, double D, double V)
{
#ifdef CHECK_H // in checks.h
	double cmp_H, cmp_D, cmp_V;
	cmp_H = get_H_of_pi(ppoints, beta, alphas, interaction_type,
		&cmp_D, &cmp_V);
	if (fabs(H-cmp_H) > 1e-7) {
		printf("\n");
		printf("Track H=%11.7lf D=%11.7lf V=%11.7lf\n", H, D, V);
		printf("Check H=%11.7lf D=%11.7lf V=%11.7lf\n", cmp_H, cmp_D, cmp_V);
		printf("Diff  H=%11.7lf D=%11.7lf V=%11.7lf\n",
			H-cmp_H, D-cmp_D, V-cmp_V);
		printf("\n");
		exit(1);
	}
#endif
}
