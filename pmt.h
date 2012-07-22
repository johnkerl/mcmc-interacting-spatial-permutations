// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
// ================================================================

#ifndef PMT_H
#define PMT_H
#include <stdio.h>
#include "points.h"

void set_unif_rand_pmt(points_t* ppoints);
void set_wx_pmt(points_t* ppoints);
void set_wy_pmt(points_t* ppoints);
void set_wz_pmt(points_t* ppoints);
void pmt_get_cycle_counts_and_lmax(points_t* ppoints, int* counts, int *plmax);
int  pmt_get_num_two_cycles(points_t* ppoints);

void pmt_print_cyclens(points_t* ppoints);
void pmt_print_cycdec(points_t* ppoints);
void pmt_print_cycle(point_t* px);
void pmt_print_cycle_verbose(point_t* px);
void pmt_print_images(points_t* ppoints);

int  get_uncached_cycle_length(point_t* px);
void get_mean_cycle_length(points_t* ppoints,
	double* pellbar, double* pmeanspatlen,
	double* pellfbar, double* pmeanspatlenf);

void get_cycle_winding_numbers(points_t* ppoints, point_t* px,
	int* pWx, int* pWy, int* pWz);
int cycle_winds(points_t* ppoints, point_t* px);
double get_fraction_of_cycles_which_wind(points_t* ppoints);
void get_pmt_winding_numbers(points_t* ppoints, int* pWx, int* pWy, int* pWz);
void print_winding_cycles(points_t* ppoints);

double get_cycle_0_spatlen(points_t* ppoints);

void write_pmt_to_file (points_t* ppoints, FILE* fp);
void read_pmt_from_file(points_t* ppoints, FILE* fp);

#endif // PMT_H
