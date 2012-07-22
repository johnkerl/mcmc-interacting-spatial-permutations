// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stats.h"
#include "util.h"

// ================================================================
void init_time_series(time_series_t* pseries, int nalloc)
{
	pseries->samples = double_malloc_or_die(nalloc);
	pseries->nalloc  = nalloc;
	pseries->nused   = 0;
}

// ----------------------------------------------------------------
void time_series_append(time_series_t* pseries, double x)
{
	if (pseries->nused >= pseries->nalloc) {
		fprintf(stderr, "time_series_append:  array overflow.\n");
		exit(1);
	}
	pseries->samples[pseries->nused++] = x;
}

// ================================================================
double ts_find_mean(time_series_t* pseries)
{
	double sum = 0.0;
	int i;
	for (i = 0; i < pseries->nused; i++)
		sum += pseries->samples[i];
	return sum / pseries->nused;
}

// ----------------------------------------------------------------
double ts_find_sampvar(time_series_t* pseries)
{
	double mean = ts_find_mean(pseries);
	return ts_find_sampvar_given_mean(pseries, mean);
}

// ----------------------------------------------------------------
double ts_find_sampvar_given_mean(time_series_t* pseries, double mean)
{
	double sum = 0.0;
	double diff;
	int i;
	for (i = 0; i < pseries->nused; i++) {
		diff = pseries->samples[i] - mean;
		sum += diff * diff;
	}
	return sum / (pseries->nused - 1);
}

// ----------------------------------------------------------------
double ts_find_stddev(time_series_t* pseries)
{
	return sqrt(ts_find_sampvar(pseries));
}

// ----------------------------------------------------------------
double ts_find_stddev_given_mean(time_series_t* pseries, double mean)
{
	return sqrt(ts_find_sampvar_given_mean(pseries, mean));
}

// ----------------------------------------------------------------
double ts_find_sample_cov(time_series_t* pxseries, time_series_t* pyseries)
{
	double mean_x = ts_find_mean(pxseries);
	double mean_y = ts_find_mean(pyseries);
	return ts_find_sample_cov_given_means(pxseries, pyseries,
		mean_x, mean_y);
}

// ----------------------------------------------------------------
double ts_find_sample_cov_given_means(time_series_t* pxseries,
	time_series_t* pyseries, double mean_x, double mean_y)
{
	double sum = 0.0;
	int i;

	for (i = 0; i < pxseries->nused; i++)
		sum += (pxseries->samples[i]-mean_x) * (pyseries->samples[i]-mean_y);

	return sum / (pxseries->nused - 1.0);
}

// ----------------------------------------------------------------
double ts_find_sample_corr(time_series_t* pxseries, time_series_t* pyseries)
{
	double cov = ts_find_sample_cov(pxseries, pyseries);
	double mean_x   = ts_find_mean(pxseries);
	double mean_y   = ts_find_mean(pyseries);
	double stddev_x = ts_find_stddev_given_mean(pxseries, mean_x);
	double stddev_y = ts_find_stddev_given_mean(pyseries, mean_y);
	return cov / (stddev_x * stddev_y);
}

// ----------------------------------------------------------------
double ts_find_tauint_flat_spot(time_series_t* pseries)
{
	return find_tauint_flat_spot(pseries->samples, pseries->nused);
}

// ----------------------------------------------------------------
double ts_find_corrected_meanvar(time_series_t* pseries)
{
	double sampvar = ts_find_sampvar(pseries);
	double naive_meanvar = sampvar / pseries->nused;
	double tauint = ts_find_tauint_flat_spot(pseries);
	return naive_meanvar * tauint;
}

// ----------------------------------------------------------------
void free_time_series(time_series_t* pseries)
{
	free(pseries->samples);
	pseries->nalloc = 0;
	pseries->nused  = 0;
}

// ================================================================
void init_first_order_stats(first_order_stats_t* p1os, int nalloc)
{
	init_time_series(&p1os->time_series, nalloc);

	p1os->avg = 0.0;
}

// ----------------------------------------------------------------
void init_second_order_stats(second_order_stats_t* pstats, int nalloc)
{
	init_time_series(&pstats->time_series, nalloc);

	pstats->avg    = 0.0;
	pstats->stddev = 0.0;
}

// ----------------------------------------------------------------
void init_second_order_cssm_stats(second_order_cssm_stats_t* pstats, int nalloc)
{
	init_time_series(&pstats->time_series, nalloc);

	pstats->avg    = 0.0;
	pstats->stddev = 0.0;
	pstats->tauint = 0.0;
	pstats->cssm   = 0.0;
}

// ----------------------------------------------------------------
void init_second_order_pair_stats(second_order_pair_stats_t* pstats, int nalloc)
{
	init_time_series(&pstats->x_series, nalloc);
	init_time_series(&pstats->y_series, nalloc);

	pstats->mx      = 0.0;
	pstats->my      = 0.0;
	pstats->sx      = 0.0;
	pstats->sy      = 0.0;
	pstats->cov_xy  = 0.0;
	pstats->corr_xy = 0.0;
}

// ================================================================
void accum_first_order_stats(first_order_stats_t* pstats, double x)
{
	time_series_append(&pstats->time_series, x);
}

// ----------------------------------------------------------------
void accum_second_order_stats(second_order_stats_t* pstats, double x)
{
	time_series_append(&pstats->time_series, x);
}

// ----------------------------------------------------------------
void accum_second_order_cssm_stats(second_order_cssm_stats_t* pstats, double x)
{
	time_series_append(&pstats->time_series, x);
}

// ----------------------------------------------------------------
void accum_second_order_pair_stats(second_order_pair_stats_t* pstats,
	double x, double y)
{
	time_series_append(&pstats->x_series, x);
	time_series_append(&pstats->y_series, y);
}

// ================================================================
void finish_first_order_stats(first_order_stats_t* pstats)
{
	pstats->avg = ts_find_mean(&pstats->time_series);
}

// ----------------------------------------------------------------
void finish_second_order_stats(second_order_stats_t* pstats)
{
	time_series_t* pseries = &pstats->time_series;
	double mean    = ts_find_mean(pseries);
	double sampvar = ts_find_sampvar_given_mean(pseries, mean);

	pstats->avg     = mean;
	pstats->stddev  = sqrt(sampvar);
}

// ----------------------------------------------------------------
void finish_second_order_cssm_stats(second_order_cssm_stats_t* pstats)
{
	time_series_t* pseries = &pstats->time_series;
	double mean    = ts_find_mean(pseries);
	double sampvar = ts_find_sampvar_given_mean(pseries, mean);
	double naive_meanvar = sampvar / pseries->nused;
	double tauint  = ts_find_tauint_flat_spot(pseries);
	double corrected_meanvar = naive_meanvar * tauint;

	pstats->avg    = mean;
	pstats->stddev = sqrt(sampvar);
	pstats->tauint = tauint;
	pstats->cssm   = sqrt(corrected_meanvar);
}

// ----------------------------------------------------------------
void finish_second_order_pair_stats(second_order_pair_stats_t* pstats)
{
	pstats->mx = ts_find_mean(&pstats->x_series);
	pstats->my = ts_find_mean(&pstats->y_series);

	pstats->sx = ts_find_stddev_given_mean(&pstats->x_series, pstats->mx);
	pstats->sy = ts_find_stddev_given_mean(&pstats->y_series, pstats->my);

	pstats->cov_xy = ts_find_sample_cov_given_means(&pstats->x_series,
		&pstats->y_series, pstats->mx, pstats->my);
	pstats->corr_xy = pstats->cov_xy / (pstats->sx * pstats->sy);
}

// ================================================================
void free_first_order_stats(first_order_stats_t* pstats)
{
	free_time_series(&pstats->time_series);
}

// ----------------------------------------------------------------
void free_second_order_stats(second_order_stats_t* pstats)
{
	free_time_series(&pstats->time_series);
}

// ----------------------------------------------------------------
void free_second_order_cssm_stats(second_order_cssm_stats_t* pstats)
{
	free_time_series(&pstats->time_series);
}

// ----------------------------------------------------------------
void free_second_order_pair_stats(second_order_pair_stats_t* p2op)
{
	free_time_series(&p2op->x_series);
	free_time_series(&p2op->y_series);
}

// ================================================================
void print_name_and_ivalue(char* name, int value)
{
	printf("# %-28s = %4d\n", name, value);
}

// ----------------------------------------------------------------
void print_name_and_dvalue(char* name, double value)
{
	printf("# %-28s = %12.7lf\n", name, value);
}

// ----------------------------------------------------------------
void print_name2_and_dvalue(char* a, char* b, double value)
{
	char* name = malloc_or_die(strlen(a) + strlen(b) + 2);
	sprintf(name, "%s_%s", a, b);
	print_name_and_dvalue(name, value);
	free(name);
}

// ----------------------------------------------------------------
void print_name3_and_dvalue(char* a, char* b, char* c, double value)
{
	char* name = malloc_or_die(strlen(a) + strlen(b) + strlen(c) + 3);
	sprintf(name, "%s_%s_%s", a, b, c);
	print_name_and_dvalue(name, value);
	free(name);
}

// ----------------------------------------------------------------
void print_first_order_stats(first_order_stats_t* pstats, char* name)
{
	print_name2_and_dvalue("mean", name, pstats->avg);
}

// ----------------------------------------------------------------
void print_second_order_stats(second_order_stats_t* pstats, char* name)
{
	print_name2_and_dvalue("mean",   name, pstats->avg);
	print_name2_and_dvalue("stddev", name, pstats->stddev);
}

// ----------------------------------------------------------------
void print_second_order_cssm_stats(second_order_cssm_stats_t* pstats,
	char* name)
{
	print_name2_and_dvalue("mean",   name, pstats->avg);
	print_name2_and_dvalue("stddev", name, pstats->stddev);
	print_name2_and_dvalue("tauint", name, pstats->tauint);
	print_name2_and_dvalue("eta",    name, tauint_to_eta(pstats->tauint));
	print_name2_and_dvalue("cssm",   name, pstats->cssm);
}

// ----------------------------------------------------------------
void print_second_order_pair_stats(second_order_pair_stats_t* pstats,
	char* xname, char* yname)
{
	print_name2_and_dvalue("mean",   xname, pstats->mx);
	print_name2_and_dvalue("mean",   yname, pstats->my);
	print_name2_and_dvalue("stddev", xname, pstats->sx);
	print_name2_and_dvalue("stddev", yname, pstats->sy);
	print_name3_and_dvalue("cov",    xname, yname, pstats->cov_xy);
	print_name3_and_dvalue("corr",   xname, yname, pstats->corr_xy);
}

// ================================================================
void ihisto_init(ihisto_t* pihisto, int lo_bin_point, int hi_bin_point)
{
	int num_bins = hi_bin_point - lo_bin_point + 1;
	int bidx;
	if (hi_bin_point <= lo_bin_point) {
		fprintf(stderr, "ihisto_init:  malformed boundaries %d, %d.",
			lo_bin_point, hi_bin_point);
		exit(1);
	}
	pihisto->count = 0;
	pihisto->lo_bin_point = lo_bin_point;
	pihisto->hi_bin_point = hi_bin_point;
	pihisto->bins         = int_malloc_or_die(num_bins);
	for (bidx = 0; bidx < num_bins; bidx++)
		pihisto->bins[bidx] = 0;
	pihisto->num_below_lo = 0;
	pihisto->num_above_hi = 0;
	pihisto->min_in_data  = 0;
	pihisto->max_in_data  = 0;
}

// ----------------------------------------------------------------
void ihisto_incr(ihisto_t* pihisto, int value)
{
	int binidx;

	pihisto->count++;
	if (value < pihisto->lo_bin_point) {
		pihisto->num_below_lo++;
	}
	else if (value > pihisto->hi_bin_point) {
		pihisto->num_above_hi++;
	}
	else {
		binidx = value - pihisto->lo_bin_point;
		pihisto->bins[binidx]++;
	}
	if ((pihisto->count == 1) || (value < pihisto->min_in_data))
		pihisto->min_in_data = value;
	if ((pihisto->count == 1) || (value > pihisto->max_in_data))
		pihisto->max_in_data = value;
}

// ----------------------------------------------------------------
void ihisto_report(ihisto_t* pihisto, char* RV_name)
{
	int bval;
	int bidx;
	int n = pihisto->count;
	int k;

	bval = pihisto->lo_bin_point;
	k    = pihisto->num_below_lo;
	printf("# %s < %3d:   %7d / %7d = %10.8lf\n",
		RV_name, bval, k, n, (double)k/(double)n);

	for (bval = pihisto->lo_bin_point; bval <= pihisto->hi_bin_point; bval++) {
		bidx = bval - pihisto->lo_bin_point;
		k = pihisto->bins[bidx];
		printf("# %s = %3d:   %7d / %7d = %10.8lf\n",
			RV_name, bval, k, n, (double)k/(double)n);
	}

	bval = pihisto->hi_bin_point;
	k    = pihisto->num_above_hi;
	printf("# %s > %3d:   %7d / %7d = %10.8lf\n",
		RV_name, bval, k, n, (double)k/(double)n);

	printf("# %s_min    = %7d\n", RV_name, pihisto->min_in_data);
	printf("# %s_max    = %7d\n", RV_name, pihisto->max_in_data);
}

// ----------------------------------------------------------------
void ihisto_free(ihisto_t* pihisto)
{
	pihisto->count        = 0;
	pihisto->lo_bin_point = 0;
	pihisto->hi_bin_point = 0;
	free(pihisto->bins);
	pihisto->bins         = 0;
	pihisto->num_below_lo = 0;
	pihisto->num_above_hi = 0;
	pihisto->min_in_data  = 0;
	pihisto->max_in_data  = 0;
}

// ================================================================
double find_tauint_flat_spot(double* xs, int N)
{
	double sumi = 0.0, sumi2 = 0.0, sumj = 0.0, sumj2 = 0.0, sumij = 0.0;
	double meani, meanj, stdi, stdj;
	double x, xi, xim, xj, xjp;
	int i, j, k, t, winsz, winszm1, denom;
	double autocorr_t, tauint;

	// Within-window sums sumi, sumi2, sumj, and sumj2 contain terms which can
	// be re-used across k.  The cross-sums sumij are different for each k and
	// must be recomputed.
	//
	// Here, compute the full-window sums.
	sumi = 0.0; sumi2 = 0.0;
	for (k = 0; k < N; k++) {
		x = xs[k];
		sumi  += x;
		sumi2 += x*x;
	}
	sumj = sumi; sumj2 = sumi2;

	tauint = 1.0;
	// Go up only to t=N-2.  The autocorr estimator for t=N-1 doesn't work
	// (only one sample), and if we haven't found a flat spot of tauint by
	// then, we aren't going to.
	for (t = 1; t < N-1; t++) {
		winsz = N - t;
		winszm1 = winsz - 1;
		denom = winszm1;
		if (winszm1 == 0)
			denom = 1;
		i = t; j = N-t-1;

		// Update the within-window sums.
		xi  = xs[ i ];  xj  = xs[ j ];
		xim = xs[i-1];  xjp = xs[j+1];
		sumi -= xim; sumi2 -= xim*xim;
		sumj -= xjp; sumj2 -= xjp*xjp;

		// Compute the cross-sum.
		sumij = 0.0;
		for (k = 0; k < winsz; k++)
			sumij += xs[k] * xs[t+k];

		// Compute the autocorrelation term for this t.
		meani = sumi / winsz; meanj = sumj / winsz;
		stdi  = sqrt((sumi2 - (sumi*sumi / winsz)) / denom);
		stdj  = sqrt((sumj2 - (sumj*sumj / winsz)) / denom);

		autocorr_t = 0.0;
		if ((stdi != 0.0) && (stdj != 0.0))
			autocorr_t = (sumij / winsz - (meani*meanj)) / (stdi*stdj);
		tauint += 2.0 * autocorr_t;

		//print '%d %11.7lf %11.7lf' % (t, autocorr_t, tauint)

		if (autocorr_t < 0.0)
			return tauint;
	}

	return tauint;
}

// ----------------------------------------------------------------
double tauint_to_eta(double tauint)
{
	double eta;

	if (tauint <= 1.0)
		eta = 0.0;
	else
		eta = (tauint - 1.0) / (tauint + 1.0);

	return eta;
}

// ----------------------------------------------------------------
double tauint_to_tauexp(double tauint)
{
	double eta = tauint_to_eta(tauint);
	double tauexp = -1.0 / log(eta);
	return tauexp;
}
