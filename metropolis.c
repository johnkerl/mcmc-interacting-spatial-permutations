// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <math.h>
#include "points.h"
#include "interactions.h"
#include "energy.h"
#include "metropolis.h"
#include "rcmrand.h"
#include "wormidx.h"
#include "worm_moves.h"
#include "util.h"

// ----------------------------------------------------------------
#define   PRINT_ZHT           0
#define   PRINT_HDVW          0
#define   PRINT_OP_CL_ACC_REJ 0
#define   PRINT_OP_CL_ACC     0

#define   NWNOJUMP            6
//#define NWNOJUMP            7

#define   NNEI                6
//#define NNEI                18
//#define NNEI                26

// ----------------------------------------------------------------
#if NWNOJUMP == 6
static int wno_jumps[NWNOJUMP][3] = {
	{ 1, 0, 0 },
	{-1, 0, 0 },
	{ 0, 1, 0 },
	{ 0,-1, 0 },
	{ 0, 0, 1 },
	{ 0, 0,-1 }
};
#endif

#if NWNOJUMP == 7
static int wno_jumps[NWNOJUMP][3] = {
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{-1, 0, 0 },
	{ 0, 1, 0 },
	{ 0,-1, 0 },
	{ 0, 0, 1 },
	{ 0, 0,-1 }
};
#endif

// ----------------------------------------------------------------
#if NNEI == 6
static int fwd_jumps[NNEI][3] = {
	{ 1, 0, 0 },
	{-1, 0, 0 },
	{ 0, 1, 0 },
	{ 0,-1, 0 },
	{ 0, 0, 1 },
	{ 0, 0,-1 }
};
#endif

#if NNEI == 18
static int fwd_jumps[NNEI][3] = {
	               { 0,-1,-1 },
	{-1,-1, 0 },   { 0,-1, 0 },   { 1,-1, 0 },
	               { 0,-1, 1 },

	{-1, 0,-1 },   { 0, 0,-1 },   { 1, 0,-1 },
	{-1, 0, 0 },                  { 1, 0, 0 },
	{-1, 0, 1 },   { 0, 0, 1 },   { 1, 0, 1 },

	               { 0, 1,-1 },
	{-1, 1, 0 },   { 0, 1, 0 },   { 1, 1, 0 },
	               { 0, 1, 1 },
};
#endif

#if NNEI == 26
static int fwd_jumps[NNEI][3] = {
	{-1,-1,-1 },   { 0,-1,-1 },   { 1,-1,-1 },
	{-1,-1, 0 },   { 0,-1, 0 },   { 1,-1, 0 },
	{-1,-1, 1 },   { 0,-1, 1 },   { 1,-1, 1 },

	{-1, 0,-1 },   { 0, 0,-1 },   { 1, 0,-1 },
	{-1, 0, 0 },                  { 1, 0, 0 },
	{-1, 0, 1 },   { 0, 0, 1 },   { 1, 0, 1 },

	{-1, 1,-1 },   { 0, 1,-1 },   { 1, 1,-1 },
	{-1, 1, 0 },   { 0, 1, 0 },   { 1, 1, 0 },
	{-1, 1, 1 },   { 0, 1, 1 },   { 1, 1, 1 }
};
#endif

// ----------------------------------------------------------------
static void try_SO_swap(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV,
	metro_type_stats_t* pmetro_type_stats);

static int try_worm_open(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan, int* wno_jump,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats);

static int try_worm_head_swap(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats);

static int try_worm_tail_swap(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats);

static int try_worm_close(points_t* ppoints, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats);

static void print_worm(point_t* w);

// ----------------------------------------------------------------
static inline point_t* random_neighbor(points_t* ppoints, point_t* a)
{
	point_t* b = a;
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int bi, bj, bk;
	int* fwd_jump;
	int jidx;

	// Choose y from a bubble around pi(x).
	// Select random jump.
	jidx = (unsigned)(IRANDOM() >> 16) % NNEI;
	fwd_jump = fwd_jumps[jidx];
	bi = a->c.selfi + fwd_jump[0];
	bj = a->c.selfj + fwd_jump[1];
	bk = a->c.selfk + fwd_jump[2];

#ifdef PERIODIC_BOUNDARY_CONDITIONS
	if      (bi <   0) bi += nx;
	else if (bi >= nx) bi -= nx;

	if      (bj <   0) bj += ny;
	else if (bj >= ny) bj -= ny;

	if      (bk <   0) bk += nz;
	else if (bk >= nz) bk -= nz;
#else
	if      (bi <   0) return a;
	else if (bi >= nx) return a;

	if      (bj <   0) return a;
	else if (bj >= ny) return a;

	if      (bk <   0) return a;
	else if (bk >= nz) return a;
#endif

	b = &ppoints->lattice[bi][bj][bk];

	return b;
}

// ----------------------------------------------------------------
void SO_sweep(points_t* ppoints, mcmc_params_t* pparams,
	double* alphas, double* pH, double* pD, double* pV,
	metro_stats_t* pmetro_stats)
{
	int  nx = ppoints->dims.nx;
	int  ny = ppoints->dims.ny;
	int  nz = ppoints->dims.nz;
	int  i, j, k;
	point_t *x, *pix, *y, *piy;
	int   xi,    xj,    xk;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				// Randomized or sequential selection of site x.
				if (pparams->do_random_sweep) {
					xi = IMODRANDOM(nx);
					xj = IMODRANDOM(ny);
					xk = IMODRANDOM(nz);
				}
				else {
					xi = i;
					xj = j;
					xk = k;
				}
				x = &ppoints->lattice[xi][xj][xk];

				// Locate pi(x).
				pix = x->pfwd;

				// Choose pi(y) near pi(x).
				piy = random_neighbor(ppoints, pix);

				// Find y.
				y = piy->pbwd;

				// Conditionally swap.
				try_SO_swap(ppoints, x, pix, y, piy,
					pparams, alphas, pH, pD, pV, &pmetro_stats->SO_stats);
			}
		}
	}
}

// ----------------------------------------------------------------
static void try_SO_swap(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, metro_type_stats_t* pmetro_type_stats)
{
	int dsq_x_piy, dsq_pix_y;
	int dxy = 0, dyx = 0, x_y_same_cycle = 0;
	double Delta_H = 0.0;
	double Delta_D = 0.0;
	double Delta_V = 0.0;
	double E;

	if (piy == pix) {
		pmetro_type_stats->num_self++;
		return;
	}

	// D terms
	// Compute squared distances |x - pi(y)|^2 and |pi(x) - y|^2
	dsq_x_piy = get_distance_squaredp(&x->c, &piy->c, &ppoints->dims);
	dsq_pix_y = get_distance_squaredp(&pix->c, &y->c, &ppoints->dims);
	Delta_D = 0.25 / pparams->beta *
		(dsq_x_piy + dsq_pix_y - x->fwd_dsq - y->fwd_dsq);

	// V terms
	// This routine must be called even if there are no interactions,
	// since it also computes dxy, dyx, and x_y_same_cycle.
	Delta_V = get_Delta_V_SO(ppoints, x, pix, y, piy,
		pparams->beta, alphas, pparams->interaction_type,
		&dxy, &dyx, &x_y_same_cycle);

	// Total energy change
	Delta_H = Delta_D + Delta_V;

	// Metropolis decision.
	E = exp(-Delta_H);
	if (URANDOM() < E) {

		x   ->pfwd = piy;
		piy ->pbwd = x;
		y   ->pfwd = pix;
		pix ->pbwd = y;

		x->fwd_dsq = dsq_x_piy;
		y->fwd_dsq = dsq_pix_y;

		x->fwd_d   = sqrt(x->fwd_dsq);
		y->fwd_d   = sqrt(y->fwd_dsq);

		update_cycinfo(ppoints, x, y, dxy, dyx, x_y_same_cycle);

		*pH  += Delta_H;
		*pD  += Delta_D;
		*pV  += Delta_V;

		pmetro_type_stats->num_change++;
	}
	else {
		pmetro_type_stats->num_keep++;
	}
}

// ----------------------------------------------------------------
void R_sweep(points_t* ppoints, mcmc_params_t* pparams)
{
	point_t* x;
	point_t* pcurr;
	point_t* ptemp;
	int dsqtemp;
	double dtemp;
	cycinfo_t* pci;

	if (!pparams->allow_cycle_reversals)
		return;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {

		// Don't reverse this cycle.
		if (URANDOM() < 0.5)
			continue;

		x = pci->psite;

		// One-cycles and two-cycles have trivial reversals.
		if (x->pfwd == x)
			continue;
		if (x->pfwd->pfwd == x)
			continue;

		// Change the fwd_dsq's.
		pcurr   = x->pbwd;
		dsqtemp = pcurr->fwd_dsq;
		dtemp   = pcurr->fwd_d;
		while (pcurr != x) {
			pcurr->fwd_dsq = pcurr->pbwd->fwd_dsq;
			pcurr->fwd_d   = pcurr->pbwd->fwd_d;
			pcurr = pcurr->pbwd;
		};
		x->fwd_dsq = dsqtemp;
		x->fwd_d   = dtemp;

		// Reverse the permutation arrows.
		pcurr = x;
		do {
			ptemp = pcurr->pfwd;
			pcurr->pfwd = pcurr->pbwd;
			pcurr->pbwd = ptemp;
			pcurr = pcurr->pbwd;
		} while (pcurr != x);
	}
}

// ----------------------------------------------------------------
// See my dissertation for details.

// ----------------------------------------------------------------
// SO
// Before:          # After:
//  x ---> pi(x)    #  x      pi(x)
//                  #  |       ^
//                  #  |       |
//                  #  v       |
//  pi(y) <--- y    #  pi(y)   y
// Sequential choice of x.
// Random choice of pi(y) near pi(x).
// ----------------------------------------------------------------
// OPEN                               | CLOSE
// Before:          # After:          | Before:         # After:
//  x  --->  pi(x)  #  x       pi(x)  | piinv(w)  pi(w) # piinv(w) --> pi(w)
//                  #   \       ^     |      \       ^  #
//   +-----+        #    \     /      |       \     /   #    +-----+
//   |     |        #     v   /       |        v   /    #    |     |
//   +-w <-+        #       w         |          w      #    +-w <-+
// Uniform random choice of x.        |  No randomness in choice of points.
// ----------------------------------------------------------------
// HEAD SWAP                          | TAIL SWAP
// Before:          # After:          | Before:        # After:
//  x --->  pi(x)   #  x     pi(x)    |  x --->  pi(x) #  x      pi(x)
//                  #  |       ^      |                #  |       ^
//                  #  |       |      |                #  |       |
//                  #  v       |      |                #  v       |
//  w <--- piinv(w) #  w    piinv(w)  | pi(w) <---  w  # pi(w)    w
// Random choice of x near piinv(w).  | Random choice of pi(x) near pi(w).

void worm_sweep(points_t* ppoints, mcmc_params_t* pparams,
	int do_print_worm, int do_print_wormspan,
	double* alphas, double* pH, double* pD, double* pV,
	metro_stats_t* pmetro_stats)
{
	double W = 0.0;
	int stepi = 0, nclrej = 0;
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	point_t* x = 0, * pix = 0, * piiw = 0, * piw = 0;
	int changed;
	coord_t zheadc = {0, 0, 0}, ztailc = {0, 0, 0};
	coord_t* zhead = &zheadc;
	coord_t* ztail = &ztailc;
	int zwormspan  = 0;
	int* wno_jump = 0;

	point_t* w = &ppoints->wormhole;

	do {
		stepi++;

		if (do_print_worm) {
			print_worm(w);
		}
#if PRINT_HDVW
		printf("%11.7lf %11.7lf %11.7lf %11.7lf# H D V W\n", *pH, *pD, *pV, W);
#endif

		// All loops are closed:  open a loop at a random lattice point.
		if (w->pfwd == w) {
			int xi = IMODRANDOM(nx);
			int xj = IMODRANDOM(ny);
			int xk = IMODRANDOM(nz);
			wno_jump = wno_jumps[(unsigned)(IRANDOM() >> 16) % NWNOJUMP];

			if (do_print_wormspan) {
				printf("%11.7lf %11.7lf # wsl\n",
					-ppoints->dims.nx/2.0, (double)ppoints->dims.nx);
			}

			x   = &ppoints->lattice[xi][xj][xk];
			pix = x->pfwd;
			changed = try_worm_open(ppoints, x, w, zhead, ztail, &zwormspan,
				wno_jump, pparams, alphas, pH, pD, pV, &W,
				&pmetro_stats->open_stats);
#if PRINT_OP_CL_ACC_REJ
			if (changed)
				printf("\n");
			printf("OPEN  %s wnoj = %3d %3d %3d\n",
				changed ? "accept" : "reject",
				wno_jump[0], wno_jump[1], wno_jump[2]);
#endif
#if PRINT_OP_CL_ACC
			if (changed) {
				printf("OPEN  accept wnoj = %3d %3d %3d\n",
					wno_jump[0], wno_jump[1], wno_jump[2]);
			}
#endif
		}

		// There is an open loop (a worm).
		else {
			// Detailed balance requires the following:
			// * Try closes 1/N of the time.
			// * Else, try head swaps and tail swaps equally often.
			// I.e. for U uniform on [0,1]:
			// * U in [0,    cut1]:  try head swap.
			// * U in [cut1, cut2]:  try tail swap.
			// * U in [cut2,    1]:  try close.

			double cut2 = 1.0 - 1.0 / ppoints->N;
			double cut1 = cut2/2;
			double U = URANDOM();

			if (do_print_wormspan) {
				double wsl = get_distancep(&w->pbwd->c, &w->pfwd->c,
					&ppoints->dims);
				printf("%11.7lf %11.7lf # wsl\n",
					wsl, (double)ppoints->dims.nx);
			}

			if (U < cut1) {
				piiw = w->pbwd;
				x    = random_neighbor(ppoints, piiw);
				pix  = x->pfwd;
				changed = try_worm_head_swap(ppoints, x, w, zhead, ztail,
					&zwormspan, pparams, alphas, pH, pD, pV, &W,
					&pmetro_stats->HS_stats);
			}

			else if (U < cut2) {
				piw = w->pfwd;
				pix = random_neighbor(ppoints, piw);
				x   = pix->pbwd;
				changed = try_worm_tail_swap(ppoints, x, w, zhead, ztail,
					&zwormspan, pparams, alphas, pH, pD, pV, &W,
					&pmetro_stats->TS_stats);
			}

			else {
				piiw = w->pbwd;
				piw  = w->pfwd;
				changed = try_worm_close(ppoints, w, zhead, ztail, &zwormspan,
					pparams, alphas, pH, pD, pV, &W,
					&pmetro_stats->close_stats);

				if (!changed)
					nclrej++;
#if PRINT_OP_CL_ACC_REJ
				if (changed) {
					int Wx, Wy, Wz;
					get_pmt_winding_numbers(ppoints, &Wx, &Wy, &Wz);
					printf("CLOSE accept W    = %3d %3d %3d\n", Wx, Wy, Wz);
				}
				else {
					printf("CLOSE reject wnoj = %3d %3d %3d\n",
						wno_jump[0], wno_jump[1], wno_jump[2]);
				}
#endif
#if PRINT_OP_CL_ACC
				if (changed) {
					int Wx, Wy, Wz;
					get_pmt_winding_numbers(ppoints, &Wx, &Wy, &Wz);
					printf("CLOSE accept W    = %3d %3d %3d\n", Wx, Wy, Wz);
				}
#endif
			}
		}
	} while (w->pfwd != w);
}

// ----------------------------------------------------------------
// OPEN
// Before:          # After:
//  x  --->  pi(x)  #  x       pi(x)
//                  #   \       ^
//   +-----+        #    \     /
//   |     |        #     v   /
//   +-w <-+        #       w
// Uniform random choice of x.

static int try_worm_open(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan, int* wno_jump,
	mcmc_params_t* pparams, double* alphas,
	double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats)
{
	point_t* pix = x->pfwd;
	coord_t newzhead;
	coord_t newztail;
	int     newzwormspan;
	int dxw = 0, dwx = 0, x_w_same_cycle = 0;
	double Delta_H = 0.0;
	double Delta_D = 0.0;
	double Delta_V = 0.0;
	double Delta_W = 0.0;
	double E;

	// D terms
	Delta_D = -0.25 / pparams->beta * x->fwd_dsq;

	// V terms
	Delta_V = get_Delta_V_worm(ppoints, x, pix, w, w,
		WORM_OPEN, pparams->beta, alphas, pparams->interaction_type,
		&dxw, &dwx, &x_w_same_cycle);

	// W terms

	//   Lift to Z^3.  Then, bias by L times +/- ihat, jhat, khat, or 0.
	newzhead = x  ->c;
	newztail = pix->c;
	newzhead.selfi += ppoints->dims.nx * wno_jump[0];
	newzhead.selfj += ppoints->dims.ny * wno_jump[1];
	newzhead.selfk += ppoints->dims.nz * wno_jump[2];

	newzwormspan = get_Zd_distance_squaredp(&newzhead, &newztail);
	//Delta_W = newzwormspan * pparams->gamma / ppoints->dims.nx;
	Delta_W = sqrt(newzwormspan) * pparams->gamma / ppoints->dims.nx;

	// Total energy change
	Delta_H = Delta_D + Delta_V + Delta_W;
#if PRINT_HDVW
	printf("OP %11.7lf %11.7lf %11.7lf %11.7lf# Delta H D V W\n",
		Delta_H, Delta_D, Delta_V, Delta_W);
#endif

	// Metropolis decision.
	E = exp(-Delta_H);
	if (URANDOM() < E) {
		x  ->pfwd = w;
		w  ->pbwd = x;
		w  ->pfwd = pix;
		pix->pbwd = w;

		x->fwd_dsq = 0;
		w->fwd_dsq = 0;
		x->fwd_d   = 0.0;
		w->fwd_d   = 0.0;

		update_cycinfo(ppoints, x, w, dxw, dwx, x_w_same_cycle);

		*zhead = newzhead;
		*ztail = newztail;
		*pzwormspan = newzwormspan;

#if PRINT_ZHT
		printf("OPEN  zhead %3d %3d %3d ztail %3d %3d %3d wnoj %3d %3d %3d zd %11.7lf\n",
			zhead->selfi, zhead->selfj, zhead->selfk,
			ztail->selfi, ztail->selfj, ztail->selfk,
			wno_jump[0], wno_jump[1], wno_jump[2], sqrt(newzwormspan));
#endif // PRINT_ZHT

		*pH  += Delta_H;
		*pD  += Delta_D;
		*pV  += Delta_V;
		*pW  += Delta_W;

		pmetro_type_stats->num_change++;
		return 1;
	}
	else {
		pmetro_type_stats->num_keep++;
		return 0;
	}
}

// ----------------------------------------------------------------
// HEAD SWAP
// Before:          # After:
//  x --->  pi(x)   #  x     pi(x)
//                  #  |       ^
//                  #  |       |
//                  #  v       |
//  w <--- piinv(w) #  w    piinv(w)
// Random choice of x near piinv(w).

static int try_worm_head_swap(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan, mcmc_params_t* pparams,
	double* alphas, double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats)
{
	int dsq_piiw_pix;
	point_t* piiw = w->pbwd;
	point_t* pix  = x->pfwd;
	coord_t new_head_diff;
	coord_t newzhead;
	int     newzwormspan;
	int d_x_piiw = 0, d_piiw_x = 0, x_piiw_same_cycle = 0;
	double Delta_H = 0.0;
	double Delta_D = 0.0;
	double Delta_V = 0.0;
	double Delta_W = 0.0;
	double E;

	if (x == piiw) {
		pmetro_type_stats->num_self++;
		return 1;
	}

	// D terms
	dsq_piiw_pix = get_distance_squaredp(&piiw->c, &pix->c, &ppoints->dims);
	Delta_D = 0.25 / pparams->beta * (dsq_piiw_pix - x->fwd_dsq);

	// V terms
	Delta_V = get_Delta_V_worm(ppoints, x, pix, piiw, w,
		WORM_HEAD_SWAP, pparams->beta, alphas, pparams->interaction_type,
		&d_x_piiw, &d_piiw_x, &x_piiw_same_cycle);

	// W terms

	//   Update zhead by ZCM x-piiw
	new_head_diff   = get_coord_diff(&x->c, &piiw->c, &ppoints->dims);
	newzhead        = *zhead;
	newzhead.selfi += new_head_diff.selfi;
	newzhead.selfj += new_head_diff.selfj;
	newzhead.selfk += new_head_diff.selfk;
	newzwormspan = get_Zd_distance_squaredp(&newzhead, ztail);

	//Delta_W = (newzwormspan - *pzwormspan) * pparams->gamma / ppoints->dims.nx;
	Delta_W = (sqrt(newzwormspan) - sqrt(*pzwormspan))
		* pparams->gamma / ppoints->dims.nx;

	// Total energy change
	Delta_H = Delta_D + Delta_V + Delta_W;
#if PRINT_HDVW
	printf("HS %11.7lf %11.7lf %11.7lf %11.7lf# Delta H D V W\n",
		Delta_H, Delta_D, Delta_V, Delta_W);
#endif

	// Metropolis decision.
	E = exp(-Delta_H);
	if (URANDOM() < E) {
		x   ->pfwd = w;
		w   ->pbwd = x;
		piiw->pfwd = pix;
		pix ->pbwd = piiw;

		x   ->fwd_dsq = 0;
		piiw->fwd_dsq = dsq_piiw_pix;
		x   ->fwd_d   = 0.0;
		piiw->fwd_d   = sqrt(dsq_piiw_pix);

		update_cycinfo(ppoints, x, piiw, d_x_piiw, d_piiw_x, x_piiw_same_cycle);

		*zhead = newzhead;
		*pzwormspan = newzwormspan;

#if PRINT_ZHT
		printf("HS    zhead %3d %3d %3d ztail %3d %3d %3d diff %3d %3d %3d zd %11.7lf\n",
			zhead->selfi, zhead->selfj, zhead->selfk,
			ztail->selfi, ztail->selfj, ztail->selfk,
			new_head_diff.selfi, new_head_diff.selfj, new_head_diff.selfk,
			sqrt(newzwormspan));
#endif // PRINT_ZHT

		*pH  += Delta_H;
		*pD  += Delta_D;
		*pV  += Delta_V;
		*pW  += Delta_W;

		pmetro_type_stats->num_change++;
		return 1;
	}
	else {
		pmetro_type_stats->num_keep++;
		return 0;
	}
}

// ----------------------------------------------------------------
// TAIL SWAP
// Before:        # After:
//  x --->  pi(x) #  x      pi(x)
//                #  |       ^
//                #  |       |
//                #  v       |
// pi(w) <---  w  # pi(w)    w
// Random choice of pi(x) near pi(w).

static int try_worm_tail_swap(points_t* ppoints, point_t* x, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan, mcmc_params_t* pparams,
	double* alphas, double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats)
{
	int dsq_x_piw;
	point_t* piw = w->pfwd;
	point_t* pix = x->pfwd;
	coord_t new_tail_diff;
	coord_t newztail;
	int     newzwormspan;
	int dxw = 0, dwx = 0, x_w_same_cycle = 0;
	double Delta_H = 0.0;
	double Delta_D = 0.0;
	double Delta_V = 0.0;
	double Delta_W = 0.0;
	double E;

	if (x == w) {
		pmetro_type_stats->num_self++;
		return 1;
	}

	// D terms
	dsq_x_piw = get_distance_squaredp(&x->c, &piw->c, &ppoints->dims);
	Delta_D = 0.25 / pparams->beta * (dsq_x_piw - x->fwd_dsq);

	// V terms
	Delta_V = get_Delta_V_worm(ppoints, x, pix, w, piw,
		WORM_TAIL_SWAP, pparams->beta, alphas, pparams->interaction_type,
		&dxw, &dwx, &x_w_same_cycle);

	// W terms

	//   Update ztail by ZCM pix-piw
	new_tail_diff   = get_coord_diff(&pix->c, &piw->c, &ppoints->dims);
	newztail        = *ztail;
	newztail.selfi += new_tail_diff.selfi;
	newztail.selfj += new_tail_diff.selfj;
	newztail.selfk += new_tail_diff.selfk;
	newzwormspan = get_Zd_distance_squaredp(zhead, &newztail);

	//Delta_W = (newzwormspan - *pzwormspan) * pparams->gamma / ppoints->dims.nx;
	Delta_W = (sqrt(newzwormspan) - sqrt(*pzwormspan))
		* pparams->gamma / ppoints->dims.nx;

	// Total energy change
	Delta_H = Delta_D + Delta_V + Delta_W;
#if PRINT_HDVW
	printf("TS %11.7lf %11.7lf %11.7lf %11.7lf# Delta H D V W\n",
		Delta_H, Delta_D, Delta_V, Delta_W);
#endif

	// Metropolis decision.
	E = exp(-Delta_H);
	if (URANDOM() < E) {
		x   ->pfwd = piw;
		piw ->pbwd = x;
		w   ->pfwd = pix;
		pix ->pbwd = w;

		x->fwd_dsq = dsq_x_piw;
		w->fwd_dsq = 0;
		x->fwd_d   = sqrt(dsq_x_piw);
		w->fwd_d   = 0.0;

		update_cycinfo(ppoints, x, w, dxw, dwx, x_w_same_cycle);

		*ztail = newztail;
		*pzwormspan = newzwormspan;

#if PRINT_ZHT
		printf("TS    zhead %3d %3d %3d ztail %3d %3d %3d diff %3d %3d %3d zd %11.7lf\n",
			zhead->selfi, zhead->selfj, zhead->selfk,
			ztail->selfi, ztail->selfj, ztail->selfk,
			new_tail_diff.selfi, new_tail_diff.selfj, new_tail_diff.selfk,
			sqrt(newzwormspan));
#endif // PRINT_ZHT

		*pH  += Delta_H;
		*pD  += Delta_D;
		*pV  += Delta_V;
		*pW  += Delta_W;

		pmetro_type_stats->num_change++;
		return 1;
	}
	else {
		pmetro_type_stats->num_keep++;
		return 0;
	}
}

// ----------------------------------------------------------------
// CLOSE
// Before:         # After:
// piinv(w)  pi(w) # piinv(w) --> pi(w)
//      \       ^  #
//       \     /   #    +-----+
//        v   /    #    |     |
//          w      #    +-w <-+
// No randomness in choice of points.

static int try_worm_close(points_t* ppoints, point_t* w,
	coord_t* zhead, coord_t* ztail, int* pzwormspan, mcmc_params_t* pparams,
	double* alphas, double* pH, double* pD, double* pV, double* pW,
	metro_type_stats_t* pmetro_type_stats)
{
	int dsq_piiw_piw;
	point_t* piiw = w->pbwd;
	point_t* piw  = w->pfwd;
	int d_piiw_w = 0, d_w_piiw = 0, piiw_w_same_cycle = 0;
	double Delta_H = 0.0;
	double Delta_D = 0.0;
	double Delta_V = 0.0;
	double Delta_W = 0.0;
	double E;

	if (piiw == w) {
		pmetro_type_stats->num_self++;
		return 1;
	}

	// D terms
	dsq_piiw_piw = get_distance_squaredp(&piiw->c, &piw->c, &ppoints->dims);
	Delta_D = 0.25 / pparams->beta * dsq_piiw_piw;

	// V terms
	Delta_V = get_Delta_V_worm(ppoints, piiw, w, w, piw,
		WORM_CLOSE, pparams->beta, alphas, pparams->interaction_type,
		&d_piiw_w, &d_w_piiw, &piiw_w_same_cycle);

	// W terms
	//Delta_W = -*pzwormspan * pparams->gamma / ppoints->dims.nx;
	Delta_W = -sqrt(*pzwormspan) * pparams->gamma / ppoints->dims.nx;

	// Total energy change
	Delta_H = Delta_D + Delta_V + Delta_W;
#if PRINT_HDVW
	printf("CL %11.7lf %11.7lf %11.7lf %11.7lf# Delta H D V W\n",
		Delta_H, Delta_D, Delta_V, Delta_W);
#endif

	// Metropolis decision.
	E = exp(-Delta_H);
	if (URANDOM() < E) {
		piiw->pfwd = piw;
		piw ->pbwd = piiw;
		w   ->pfwd = w;
		w   ->pbwd = w;

		piiw->fwd_dsq = dsq_piiw_piw;
		w   ->fwd_dsq = 0;
		piiw->fwd_d   = sqrt(dsq_piiw_piw);
		w   ->fwd_d   = 0.0;

		update_cycinfo(ppoints, piiw, w, d_piiw_w, d_w_piiw, piiw_w_same_cycle);

#if PRINT_ZHT
		printf("CLOSE zhead %3d %3d %3d ztail %3d %3d %3d dz %11.7lf\n",
			zhead->selfi, zhead->selfj, zhead->selfk,
			ztail->selfi, ztail->selfj, ztail->selfk,
			sqrt(*pzwormspan));
		printf("\n");
#endif // PRINT_ZHT

		zhead->selfi = 0;
		zhead->selfj = 0;
		zhead->selfk = 0;
		ztail->selfi = 0;
		ztail->selfj = 0;
		ztail->selfk = 0;
		*pzwormspan  = 0;

		*pH  += Delta_H;
		*pD  += Delta_D;
		*pV  += Delta_V;
		*pW  += Delta_W;

		pmetro_type_stats->num_change++;
		return 1;
	}
	else {
		pmetro_type_stats->num_keep++;
		return 0;
	}
}

// ----------------------------------------------------------------
static void print_worm(point_t* w)
{
	point_t* pcurr = w->pfwd;

	print_point("<", w, "");
	if (w->pfwd != w)
		printf(" ");
	while (pcurr != w) {
		print_point(" ", pcurr, "");
		pcurr = pcurr->pfwd;
	}
	if (w->pfwd != w)
		print_point(" ", w, "");
	printf("> (%d)\n", w->pcycinfo->cyclen);
}
