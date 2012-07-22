// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef MCMC_RVS_H
#define MCMC_RVS_H

#include "points.h"
#include "mcmc_params.h"
#include "stats.h"

// ----------------------------------------------------------------
#define WX_HISTO_MAX  10
#define WX_HISTO_MIN -10
// If you change this number, change the code below which populates the
// array:
#define NUM_ANGLE_STATS_SITES     4

typedef struct _ellxy_stats_t {
	// Per realization:
	point_t* px;
	point_t* py;
	second_order_pair_stats_t pair_stats;
} ellxy_stats_t;

typedef struct _jump_stats_t {
	point_t* ppoint;
	// Per realization:
	// dsq is already in prvs->jump_stats[k].ppoint->fwd_dsq.
	// Means:
	second_order_stats_t sostats;
} jump_stats_t;

typedef struct _kink_stats_t {
	point_t* ppoint;

	// Per realization:
	double   dot;
	int      cos_angle_defined;
	double   cos_angle;

	// Means:
	second_order_stats_t dot_sostats;
	second_order_stats_t cos_angle_sostats;
} kink_stats_t;

typedef struct _angle_stats_t {
	point_t* px;
	point_t* py;

	// Per realization:
	double   dot;
	int      cos_angle_defined;
	double   cos_angle;

	// Means:
	second_order_stats_t dot_sostats;
	second_order_stats_t cos_angle_sostats;
} angle_stats_t;

// ----------------------------------------------------------------
typedef struct _mcmc_rvs_t {
	int      num_acc;

	// Per realization:
	double   H;
	double   D;
	double   V;

	int      r2;
	int      lmax;

	double   cycle_0_spatlen;
	double   ellbar;
	double   meanspatlen;
	double   recipmeanspatlen;
	double   ellfbar;      // hoping to find a two-sided correlation length.
	double   meanspatlenf; // "Finite" versions for below T_c ...
	int      Wx;
	int      Wy;
	int      Wz;
	double   winding_number; // W^2 := Wx^2 + Wy^2 + Wz^2
	double   fS;             // f_S = <W^2> L^2 / 3 beta N.
	double   fW;             // Fraction of sites in cycles which wind.
	int      ellx;
	int      elly;
	double   jumplenbar;
	double   maxjumplen;

	// Means:
	int      rho_sum_count;
	double * unitk;
	double * rho;
	double * rho_sum;
	double * rho_avg;
	double * rho_fit;
	double * rho_dev;
	int    * counts;
	double * counts_sum;
	double * counts_avg;
	double   rho_infty_T;
	double * delta_rho;
	double * delta_rho_avg;

	// The data types for the following are defined in stats.h.
	second_order_cssm_stats_t H_2cos;
	first_order_stats_t       D_1os;
	first_order_stats_t       V_1os;

	second_order_cssm_stats_t r2_2cos;
	second_order_cssm_stats_t lmax_2cos;

	second_order_stats_t      jumplenbar_2os;
	second_order_stats_t      ell0_2os;
	second_order_cssm_stats_t ellbar_2cos;
	second_order_stats_t      cycle_0_spatlen_2os;
	second_order_cssm_stats_t meanspatlen_2cos;
	second_order_cssm_stats_t recipmeanspatlen_2cos;
	second_order_stats_t      ellfbar_2os;
	second_order_stats_t      meanspatlenf_2os;
	second_order_stats_t      winding_number_2os;
	second_order_cssm_stats_t fS_2cos;
	second_order_cssm_stats_t fW_2cos;

	ellxy_stats_t ellxy_stats;
	jump_stats_t  jump0_stats;
	kink_stats_t  kink_stats;
	angle_stats_t angle_stats[NUM_ANGLE_STATS_SITES];

	ihisto_t Wx_histo;

} mcmc_rvs_t;

// Verbosity-control variables.
typedef struct _rv_verbosity_t {
	int sweepi;

	int pmt;
	int cycdec;

	int rho;
	int counts;

	// For printing counts.  Setting ello and elhi to zero means print cycle
	// counts for all ell.
	// int ello;
	// int elhi;

	int H;
	int tauint_Hs;

	int r2;

	int lmax;
	int ell;
	int ellbar;
	int cycle_0_spatlen;
	int meanspatlen;
	int recipmeanspatlen;
	int winding_cycle;
	int winding_number;

	int jump_length;
	int jump_target;
	int jump_vector;
	int kink_dot;
	int kink_cos_angle;
	int dot;
	int cos_angle;
} rv_verbosity_t;

void allocate_rvs(mcmc_rvs_t* prvs, mcmc_params_t* pparams,
	points_t* ppoints);
void compute_rvs(mcmc_rvs_t* prvs, points_t* ppoints, mcmc_params_t* pparams,
	double H, double D, double V);

void accumulate_rvs(mcmc_rvs_t* prvs, points_t* ppoints,
	mcmc_params_t* pparams, double H, double D, double V);
void display_realization_rvs(double H, double D, double V,
	mcmc_rvs_t* prvs, points_t* ppoints, rv_verbosity_t* pvb);
void compute_rv_means(mcmc_rvs_t* prvs, points_t* ppoints);
void display_rv_means(mcmc_rvs_t* prvs, mcmc_params_t* pparams,
	points_t* ppoints, rv_verbosity_t* pvb);
void free_rvs(mcmc_rvs_t* prvs);

#endif // MCMC_RVS_H
