// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include "dotplot.h"

static void print_points_and_cycles_aux(points_t* ppoints,
	int do_min_max, int min_ell, int max_ell);

// ----------------------------------------------------------------
void print_points_and_cycles(points_t* ppoints)
{
	print_points_and_cycles_aux(ppoints, 0, -1, ppoints->N+1);
}

// ----------------------------------------------------------------
void print_points_and_cycles_min_max(points_t* ppoints,
	int min_ell, int max_ell)
{
	print_points_and_cycles_aux(ppoints, 1, min_ell, max_ell);
}

// ----------------------------------------------------------------
static void print_points_and_cycles_aux(points_t* ppoints,
	int do_min_max, int min_ell, int max_ell)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	point_t* a, * pia;
	int i, j, k;
	int print_line;
	double x1, y1, z1, x2, y2, z2;

	printf("#begin lines\n");
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				a = &ppoints->lattice[i][j][k];
				pia = a->pfwd;

				print_line = 1;
				if (do_min_max) {
					int ell = a->pcycinfo->cyclen;
					if ((ell < min_ell) || (ell > max_ell))
						print_line = 0;
				}
				if (!print_line)
					continue;

				x1 =   a->c.selfi;
				y1 =   a->c.selfj;
				z1 =   a->c.selfk;

				x2 = pia->c.selfi;
				y2 = pia->c.selfj;
				z2 = pia->c.selfk;

				// xxx temp toric hack
				if (fabs(x1-x2) > 0.5*nx) continue;
				if (fabs(y1-y2) > 0.5*ny) continue;
				if (fabs(z1-z2) > 0.5*nz) continue;

				printf("%11.7lf %11.7lf %11.7lf %11.7lf %11.7lf %11.7lf #dot\n",
					x1, y1, z1, x2, y2, z2);


//				{
//					double x1left, x2left, x1right, x2right;
//					double y1left, y2left, y1right, y2right;
//					double z1left, z2left, z1right, z2right;
//					int    xsplit = toric_split(x1, x2, (double)nx,
//						&x1left, &x2left, &x1right, &x2right);
//					int    ysplit = toric_split(y1, y2, (double)ny,
//						&y1left, &y2left, &y1right, &y2right);
//					int    zsplit = toric_split(z1, z2, (double)nz,
//						&z1left, &z2left, &z1right, &z2right);
//					if (xsplit) {
//					}
//					else {
//					}
//					if (ysplit) {
//					}
//					else {
//					}
//					if (zsplit) {
//					}
//					else {
//					}
//				}


			}
		}
	}
	printf("#end lines\n");
}
