// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "points.h"
#include "pmt.h"
#include "wormidx.h"
#include "rcmrand.h"

// ----------------------------------------------------------------
// Steps:
// (1) Make a list of the N points and permute that.
// (2) Populate the forward pointers, backward pointers, and
//     squared distances.

void set_unif_rand_pmt(points_t* ppoints)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int N = ppoints->N;
	point_t** fwds = (point_t**)malloc_or_die(N * sizeof(point_t*));
	point_t* temp;
	int s, i, j, k;
	int unused_start = 0;
	int num_unused   = N;
	int u; // For "unused":  random index into the pool of unused images.

	// (1) Make a list of the N points and permute that.
	for (s = 0, i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			for (k = 0; k < nz; k++, s++)
				fwds[s] = &ppoints->lattice[i][j][k];

	for (s = 0; s < N; s++) {
		// Select a pseudorandom element from the pool of unused images.
		u = unused_start + IMODRANDOM(num_unused);

		// Swap it into place.
		temp    = fwds[u];
		fwds[u] = fwds[s];
		fwds[s] = temp;

		// Decrease the size of the pool by 1.
		// (Yes, unused_start and s always have the same value.  Using two
		// variables wastes neglible memory and makes the code easier to
		// understand.)
		unused_start++;
		num_unused--;
	}

	// (2) Populate the forward pointers, backward pointers, and
	//     squared distances.
	for (s = 0, i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++, s++) {
				temp = &ppoints->lattice[i][j][k];
				temp->pfwd = fwds[s];
				temp->pfwd->pbwd = temp;
				temp->fwd_dsq = get_distance_squared(
					temp->c.selfi, temp->c.selfj, temp->c.selfk,
					temp->pfwd->c.selfi, temp->pfwd->c.selfj, temp->pfwd->c.selfk,
					nx, ny, nz);
				temp->fwd_d = sqrt(temp->fwd_dsq);
			}
		}
	}

	free(fwds);
}

// ----------------------------------------------------------------
// Assumes the permutation is already the identity.

void set_wx_pmt(points_t* ppoints)
{
	int i, pi, j=0, k=0;
	int nx = ppoints->dims.nx;
	point_t* a, * b;

	for (i = 0; i < nx; i++) {
		pi = (i+1) % nx;
		a = &ppoints->lattice[ i][j][k];
		b = &ppoints->lattice[pi][j][k];

		a->pfwd = b;
		b->pbwd = a;
		a->fwd_dsq = 1;
		a->fwd_d   = 1.0;
	}
}

void set_wy_pmt(points_t* ppoints)
{
	int i=0, j, pj, k=1;
	int ny = ppoints->dims.ny;
	point_t* a, * b;

	if (ppoints->dims.nz < 2)
		return;

	for (j = 0; j < ny; j++) {
		pj = (j+1) % ny;
		a = &ppoints->lattice[i][ j][k];
		b = &ppoints->lattice[i][pj][k];

		a->pfwd = b;
		b->pbwd = a;
		a->fwd_dsq = 1;
		a->fwd_d   = 1.0;
	}
}

void set_wz_pmt(points_t* ppoints)
{
	int i=1, j=1, k, pk;
	int nz = ppoints->dims.nz;
	point_t* a, * b;

	if (ppoints->dims.nx < 2)
		return;
	if (ppoints->dims.ny < 2)
		return;

	for (k = 0; k < nz; k++) {
		pk = (k+1) % nz;
		a = &ppoints->lattice[i][j][ k];
		b = &ppoints->lattice[i][j][pk];

		a->pfwd = b;
		b->pbwd = a;
		a->fwd_dsq = 1;
		a->fwd_d   = 1.0;
	}
}

// ----------------------------------------------------------------
void pmt_get_cycle_counts_and_lmax(points_t* ppoints, int* counts, int *plmax)
{
	int i;
	int lmax = 0;
	cycinfo_t* pci;
	point_t* w = &ppoints->wormhole;

	for (i = 0; i <= ppoints->N; i++)
		counts[i] = 0;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		// Exclude the wormhole point.
		if (pci->psite == w)
			continue;
		if (pci->cyclen > lmax)
			lmax = pci->cyclen;
		counts[pci->cyclen]++;
	}

	if (plmax)
		*plmax = lmax;
}

// ----------------------------------------------------------------
int pmt_get_num_two_cycles(points_t* ppoints)
{
	int r2 = 0;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext)
		if (pci->cyclen == 2)
			r2++;

	return r2;
}

// ----------------------------------------------------------------
void pmt_print_cyclens(points_t* ppoints)
{
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext)
		printf(" %d", pci->cyclen);
	printf("\n");
}

// ----------------------------------------------------------------
void pmt_print_cycdec(points_t* ppoints)
{
	point_t* pstart;
	point_t* pcurr;
	point_t* pnext;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		pstart = pci->psite;
		pcurr  = pstart;

		// Don't print one-cycles.
		if (pstart->pfwd == pstart)
			continue;

		printf("(");
		printf("%d,%d,%d",
			pstart->c.selfi, pstart->c.selfj, pstart->c.selfk);

		while (1) {
			pnext = pcurr->pfwd;
			pcurr = pnext;
			if (pnext == pstart)
				break;
			printf(" %d,%d,%d",
				pnext->c.selfi, pnext->c.selfj, pnext->c.selfk);
		}
		printf(")");
	}
	printf("\n");
}

// ----------------------------------------------------------------
void pmt_print_cycle(point_t* px)
{
	point_t* pcurr = px;
	do {
		if (pcurr->c.selfi == WORMIDX) {
			printf("  (* * *) [%d]\n", pcurr->fwd_dsq);
		}
		else {
			printf("  (%d %d %d) [%d]\n",
				pcurr->c.selfi, pcurr->c.selfj, pcurr->c.selfk, pcurr->fwd_dsq);
		}
		pcurr = pcurr->pfwd;
	} while (pcurr != px);
}

// ----------------------------------------------------------------
void pmt_print_cycle_verbose(point_t* px)
{
	point_t* pcurr = px;
	do {
		print_point("", pcurr->pbwd, "");
		print_point(" >> ", pcurr,     " >> ");
		printf("[%d]", pcurr->fwd_dsq);
		print_point("", pcurr->pfwd, "\n");
		pcurr = pcurr->pfwd;
	} while (pcurr != px);
}

// ----------------------------------------------------------------
int get_uncached_cycle_length(point_t* px)
{
	int cycle_length = 0;
	point_t* pstart = px;
	point_t* pcurr  = pstart;
	point_t* pnext  = 0;

	while (1) {
		cycle_length++;
		pnext = pcurr->pfwd;
		pcurr = pnext;
		if (pnext == pstart)
			break;
	}
	return cycle_length;
}

// ----------------------------------------------------------------
void get_mean_cycle_length(points_t* ppoints,
	double* pellbar,  double* pmeanspatlen,
	double* pellfbar, double* pmeanspatlenf)
{
	double cutoff = sqrt((double)ppoints->N); // Cutoff for "long" cycles.
	point_t* w = &ppoints->wormhole;
	int cyclen;
	double spatlen;
	point_t* pstart;
	point_t* pcurr;
	point_t* pnext;
	double ellbar = 0.0;
	double meanspatlen = 0.0;
	double ellfbar = 0.0;
	double meanspatlenf = 0.0;
	double ellterm, spatlenterm;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		// Exclude the wormhole point.
		if (pci->psite == w)
			continue;

		spatlen = 0.0;

		cyclen = pci->cyclen;
		pstart = pci->psite;
		pcurr  = pstart;

		while (1) {
			spatlen += pcurr->fwd_d;

			pnext = pcurr->pfwd;
			pcurr = pnext;
			if (pnext == pstart)
				break;
		}

		// All sites in this cycle share the same cycle length
		// and cycle spatial length.  E.g. if we just found a 3-cycle,
		// we found 3 sites which participate in it.
		ellterm = cyclen * cyclen;
		spatlenterm = cyclen * spatlen;
		ellbar      += ellterm;
		meanspatlen += spatlenterm;
		if (cyclen < cutoff) {
			ellfbar += ellterm;
			meanspatlenf += spatlenterm;
		}
	}

	ellbar        /= ppoints->N;
	meanspatlen   /= ppoints->N;
	ellfbar       /= ppoints->N;
	meanspatlenf  /= ppoints->N;

	*pellbar       = ellbar;
	*pmeanspatlen  = meanspatlen;
	*pellfbar      = ellfbar;
	*pmeanspatlenf = meanspatlenf;
}

// ----------------------------------------------------------------
void pmt_print_images(points_t* ppoints)
{
	int i, j, k;
	point_t* ppoint;
	int num_non_triv = 0;
	for (i = 0; i < ppoints->dims.nx; i++) {
		for (j = 0; j < ppoints->dims.ny; j++) {
			for (k = 0; k < ppoints->dims.nz; k++) {
				ppoint = &ppoints->lattice[i][j][k];
				if (ppoint == ppoint->pfwd)
					continue;
				num_non_triv++;
				printf(" %d,%d,%d->%d,%d,%d[%d]",
					ppoint->c.selfi,
					ppoint->c.selfj,
					ppoint->c.selfk,
					ppoint->pfwd->c.selfi,
					ppoint->pfwd->c.selfj,
					ppoint->pfwd->c.selfk,
					ppoint->fwd_dsq);
			}
		}
	}
	if (num_non_triv == 0)
		printf("(none)\n");
	printf("\n");
}

// ----------------------------------------------------------------
// See the dissertation for details.
// L = 10:  L/2 = 5.  -L/2 = -5.  Range [-5 .. +5].
// L = 11:  L/2 = 5.  -L/2 = -5.  Range [-5 .. +5].
// L = 12:  L/2 = 6.  -L/2 = -6.  Range [-6 .. +6].
void get_cycle_winding_numbers(points_t* ppoints, point_t* px,
	int* pWx, int* pWy, int* pWz)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int half_nx = nx/2; // Same arithmetic as zero_centered_mod in points.h.
	int half_ny = ny/2;
	int half_nz = nz/2;
	int dx, dy, dz;
	int Wx, Wy, Wz;
	point_t* pcurr, * pnext;

	Wx = Wy = Wz = 0;

	pcurr  = px;
	pnext  = pcurr->pfwd;

	do {
		dx = pnext->c.selfi - pcurr->c.selfi;
		dy = pnext->c.selfj - pcurr->c.selfj;
		dz = pnext->c.selfk - pcurr->c.selfk;

		if      (dx >  half_nx) Wx--;
		else if (dx < -half_nx) Wx++;

		if      (dy >  half_ny) Wy--;
		else if (dy < -half_ny) Wy++;

		if      (dz >  half_nz) Wz--;
		else if (dz < -half_nz) Wz++;

		pcurr = pnext;
		pnext = pcurr->pfwd;

	} while (pcurr != px);

	*pWx = Wx;
	*pWy = Wy;
	*pWz = Wz;
}

// ----------------------------------------------------------------
int cycle_winds(points_t* ppoints, point_t* px)
{
	int cyWx, cyWy, cyWz;
	get_cycle_winding_numbers(ppoints, px, &cyWx, &cyWy, &cyWz);
	if (cyWx != 0 || cyWy != 0 || cyWz != 0)
		return 1;
	else
		return 0;
}

// ----------------------------------------------------------------
double get_fraction_of_cycles_which_wind(points_t* ppoints)
{
	int N = ppoints->N;
	int NW = 0;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext)
		if (cycle_winds(ppoints, pci->psite))
			NW += pci->cyclen;

	return (double)NW / (double)N;
}

// ----------------------------------------------------------------
// The winding number of a permutation is the sum of the winding numbers of all
// permutation jumps.  It is a 3-tuple:  number of wraps around the torus in
// the x, y, and z directions.

void get_pmt_winding_numbers(points_t* ppoints, int* pWx, int* pWy, int* pWz)
{
	int piWx = 0, piWy = 0, piWz = 0;
	int cyWx, cyWy, cyWz;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		get_cycle_winding_numbers(ppoints, pci->psite, &cyWx, &cyWy, &cyWz);
		piWx += cyWx;
		piWy += cyWy;
		piWz += cyWz;
	}
	*pWx = piWx;
	*pWy = piWy;
	*pWz = piWz;
}

// ----------------------------------------------------------------
void print_winding_cycles(points_t* ppoints)
{
	int Wx, Wy, Wz;
	cycinfo_t* pci;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		get_cycle_winding_numbers(ppoints, pci->psite, &Wx, &Wy, &Wz);
		if (Wx != 0 || Wy != 0 || Wz != 0) {
			printf("\n");
			printf("%6d %6d %6d # Wx Wy Wz cycle\n", Wx, Wy, Wz);
			pmt_print_cycle(pci->psite);
		}
	}
}

// ----------------------------------------------------------------
double get_cycle_0_spatlen(points_t* ppoints)
{
	double spatlen = 0.0;
	point_t* origin = &ppoints->lattice[0][0][0];
	point_t* pcurr = origin;
	do {
		spatlen += pcurr->fwd_d;
		pcurr = pcurr->pfwd;
	} while (pcurr != origin);

	return spatlen;
}

// ----------------------------------------------------------------
void write_pmt_to_file(points_t* ppoints, FILE* fp)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int i, j, k;
	point_t* px;
	unsigned char uimgs[3];
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				px = &ppoints->lattice[i][j][k];
				uimgs[0] = px->pfwd->c.selfi;
				uimgs[1] = px->pfwd->c.selfj;
				uimgs[2] = px->pfwd->c.selfk;
				write_to_file_or_die((void*)&uimgs[0], sizeof(uimgs), 1, fp);
			}
		}
	}

	// write num cycles
	// for each cycle:
	//   write psite's idx
	//   write cyclen
	// for each site:
	//   write cycidx

	// read num cycles
	// for each cycle:
	//   read psite's idx
	//   read cyclen
	// for each site:
	//   set cycidx

//num cycles
//for each cycle:
//	site idx
//	cyclen
//for each site:
//	cycidx

//xxx cycinfo to file as well:
//
//// Cycle-information data structure.
//typedef struct _cycinfo_t {
//	struct _cycinfo_t* pprev;  // Backward pointer for doubly linked list.
//	struct _cycinfo_t* pnext;  // Forward  pointer for doubly linked list.
//	point_t* psite;  // Points to one of the lattice sites in the cycle.
//	int      cyclen; // Number of lattice sites in the cycle.
//} cycinfo_t;
//
//Serialize:

}

// ----------------------------------------------------------------
void read_pmt_from_file(points_t* ppoints, FILE* fp)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int i, j, k;
	int u, v, w;
	point_t* px;
	unsigned char uimgs[3];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				read_from_file_or_die((void*)&uimgs[0], sizeof(uimgs), 1, fp);
				u = uimgs[0];
				v = uimgs[1];
				w = uimgs[2];
				px = &ppoints->lattice[i][j][k];
				px->pfwd = &ppoints->lattice[u][v][w];
				px->pfwd->pbwd = px;
				px->fwd_dsq = get_distance_squared(
					px->c.selfi, px->c.selfj, px->c.selfk,
					px->pfwd->c.selfi, px->pfwd->c.selfj, px->pfwd->c.selfk,
					nx, ny, nz);
				px->fwd_d = sqrt(px->fwd_dsq);
			}
		}
	}
}
