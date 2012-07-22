// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef STATS_H
#define STATS_H

// ----------------------------------------------------------------
typedef struct _time_series_t {
	double* samples;
	int     nalloc;
	int     nused;
} time_series_t;

void init_time_series(time_series_t* pseries, int nalloc);
void time_series_append(time_series_t* pseries, double x);
void free_time_series            (time_series_t* pseries);

double ts_find_mean(time_series_t* pseries);
double ts_find_sampvar(time_series_t* pseries);
double ts_find_sampvar_given_mean(time_series_t* pseries, double mean);
double ts_find_stddev(time_series_t* pseries);
double ts_find_stddev_given_mean(time_series_t* pseries, double mean);
double ts_find_sample_cov (time_series_t* pxseries, time_series_t* pyseries);
double ts_find_sample_cov_given_means(time_series_t* pxseries,
	time_series_t* pyseries, double mean_x, double mean_y);
double ts_find_sample_corr(time_series_t* pxseries, time_series_t* pyseries);
double ts_find_tauint_flat_spot(time_series_t* pseries);
double ts_find_corrected_meanvar(time_series_t* pseries);

// ----------------------------------------------------------------
typedef struct _first_order_stats_t {
	// Time series
	time_series_t time_series;

	// Result
	double   avg;
} first_order_stats_t;

typedef struct _second_order_stats_t {
	// Time series
	time_series_t time_series;

	// Results
	double   avg;
	double   stddev;
} second_order_stats_t;

typedef struct _second_order_cssm_stats_t {
	// Time series
	time_series_t time_series;

	// Results
	double   avg;
	double   stddev;
	double   tauint;
	double   cssm; // Corrected sample standard deviation of the sample mean.
} second_order_cssm_stats_t;

typedef struct _second_order_pair_stats_t {
	// Time series
	time_series_t x_series;
	time_series_t y_series;

	// Results
	double   mx;
	double   my;
	double   sx;
	double   sy;
	double   cov_xy;
	double   corr_xy;
} second_order_pair_stats_t;

// ----------------------------------------------------------------
void init_first_order_stats      (first_order_stats_t*       pst,  int nalloc);
void init_second_order_stats     (second_order_stats_t*      pst,  int nalloc);
void init_second_order_cssm_stats(second_order_cssm_stats_t* pst,  int nalloc);
void init_second_order_pair_stats(second_order_pair_stats_t* pst,  int nalloc);

void accum_first_order_stats      (first_order_stats_t*       pst, double x);
void accum_second_order_stats     (second_order_stats_t*      pst, double x);
void accum_second_order_cssm_stats(second_order_cssm_stats_t* pst, double x);
void accum_second_order_pair_stats(second_order_pair_stats_t* pst,
	double x, double y);

void finish_first_order_stats      (first_order_stats_t*       pst);
void finish_second_order_stats     (second_order_stats_t*      pst);
void finish_second_order_cssm_stats(second_order_cssm_stats_t* pst);
void finish_second_order_pair_stats(second_order_pair_stats_t* pst);

void free_first_order_stats      (first_order_stats_t*       pst);
void free_second_order_stats     (second_order_stats_t*      pst);
void free_second_order_cssm_stats(second_order_cssm_stats_t* pst);
void free_second_order_pair_stats(second_order_pair_stats_t* pst);

// ----------------------------------------------------------------
// Display routines:
void print_name_and_ivalue (char* name, int value);
void print_name_and_dvalue (char* name, double value);
void print_name2_and_dvalue(char* a, char* b, double value);
void print_name3_and_dvalue(char* a, char* b, char* c, double value);

void print_first_order_stats      (first_order_stats_t*       pst, char* name);
void print_second_order_stats     (second_order_stats_t*      pst, char* name);
void print_second_order_cssm_stats(second_order_cssm_stats_t* pst, char* name);
void print_second_order_pair_stats(second_order_pair_stats_t* pst,
	char*xname, char* yname);

// ----------------------------------------------------------------
typedef struct _ihisto_t {
	int   count;
	int   lo_bin_point;
	int   hi_bin_point;
	int * bins;
	int   num_below_lo;
	int   num_above_hi;
	int   min_in_data;
	int   max_in_data;
} ihisto_t;

void ihisto_init  (ihisto_t* pihisto, int lo_bin_point, int hi_bin_point);
void ihisto_incr  (ihisto_t* pihisto, int value);
void ihisto_report(ihisto_t* pihisto, char* RV_name);
void ihisto_free  (ihisto_t* pihisto);

// ----------------------------------------------------------------
double find_tauint_flat_spot(double* xs, int N);
double tauint_to_tauexp(double tauint);
double tauint_to_eta(double tauint);

#endif // STATS_H
