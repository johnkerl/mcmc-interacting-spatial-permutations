// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef ENERGY_H
#define ENERGY_H

#include "points.h"

double get_Delta_V_SO(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int interaction_type,
	int* pdxy, int* pdyx, int* px_y_same_cycle);

// Assumptions about the order of arguments:
// * Open:      get_Delta_V_worm(ppoints, x,    pix, w,    w,   ...);
// * Head swap: get_Delta_V_worm(ppoints, x,    pix, piiw, w,,  ...);
// * Tail swap: get_Delta_V_worm(ppoints, x,    pix, w,    piw, ...);
// * Close:     get_Delta_V_worm(ppoints, piiw, w,   w,    piw, ...);
double get_Delta_V_worm(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	int worm_move, double beta, double* alphas, int interaction_type,
	int* pdxy, int* pdyx, int* px_y_same_cycle);

int get_Delta_r2(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy);

double get_Delta_rell(points_t* ppoints,
	point_t* x, point_t* pix, point_t* y, point_t* piy, double* alphas,
	int* pdxy, int* pdyx, int* px_y_same_cycle);

double get_Delta_Vzw_SO(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type);

double get_Delta_Vzw_worm_open(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type);

double get_Delta_Vzw_worm_close(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type);

double get_Delta_Vzw_worm_head_swap(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type);

double get_Delta_Vzw_worm_tail_swap(points_t* ppts,
	point_t* x, point_t* pix, point_t* y, point_t* piy,
	double beta, double* alphas, int type);

// For initial setup.
// This should be called only for closed permutations.
double get_H_of_pi(points_t* ppoints, double beta, double* alphas,
	int interaction_type, double* pD, double* pV);

// For debugging.
// This should be called only for closed permutations.
void check_H(points_t* ppoints,
	double beta, double* alphas, int interaction_type,
	double H, double D, double V);

#endif // ENERGY_H
