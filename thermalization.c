// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdlib.h>
#include "thermalization.h"
#include "util.h"

// ----------------------------------------------------------------
void therm_ctl_init(therm_ctl_t* ptherm, int smoothing_window_size,
	int therm_ntp, int force_therm)
{
	int k;

	ptherm->smoothing_window_size = smoothing_window_size;
	ptherm->therm_ntp             = therm_ntp;
	ptherm->force_therm           = force_therm;

	ptherm->sweepi                = 0;
	ptherm->H_lagger              = 0.0;
	ptherm->ssH                   = 0.0;
	ptherm->smoothed_H            = 0.0;
	ptherm->smoothed_Hm1          = 0.0;
	ptherm->smoothed_Hm2          = 0.0;
	ptherm->smoothed_H_diff1      = 0.0;
	ptherm->smoothed_H_diff2      = 0.0;
	ptherm->smoothing_scale       = 0;
	ptherm->H_window = double_malloc_or_die(ptherm->smoothing_window_size);
	for (k = 0; k < ptherm->smoothing_window_size; k++)
		ptherm->H_window[k] = 0.0;
	ptherm->ntp                   = 0;
	ptherm->thermalized           = 0;
}

// ----------------------------------------------------------------
void therm_ctl_update(therm_ctl_t* ptherm, double H)
{
	int idx = ptherm->sweepi % ptherm->smoothing_window_size;

	if (ptherm->thermalized)
		return;

	ptherm->sweepi++;

	ptherm->H_lagger = ptherm->H_window[idx];
	ptherm->H_window[idx] = H;

	// Compute smoothed system energy.
	ptherm->smoothed_H_diff1 = ptherm->smoothed_H   - ptherm->smoothed_Hm1;
	ptherm->smoothed_H_diff2 = ptherm->smoothed_Hm1 - ptherm->smoothed_Hm2;
	ptherm->smoothed_Hm2     = ptherm->smoothed_Hm1;
	ptherm->smoothed_Hm1     = ptherm->smoothed_H;

	if (ptherm->sweepi < ptherm->smoothing_window_size) {
		ptherm->smoothed_H += H;
		ptherm->smoothing_scale++;
	}
	else {
		ptherm->smoothed_H -= ptherm->H_lagger;
		ptherm->smoothed_H += H;
	}
	ptherm->ssH = ptherm->smoothed_H / ptherm->smoothing_scale;

	// Detect turning points of smoothed H.
	if (ptherm->sweepi >= ptherm->smoothing_window_size) {
		if (ptherm->smoothed_H_diff1 * ptherm->smoothed_H_diff2 < 0.0) {
			ptherm->ntp++;
			if (ptherm->ntp >= ptherm->therm_ntp)
				ptherm->thermalized = 1;
		}
	}

	// This is a useful manual override for when H has few or no turning
	// points.  (My thermalization criterion is good most of the time but
	// it is not perfect -- namely, its conservatism is unbounded, so at
	// high T, it never considers the system thermalized.)
	if (ptherm->sweepi == ptherm->force_therm)
		ptherm->thermalized = 1;
}

// ----------------------------------------------------------------
void therm_ctl_free(therm_ctl_t* ptherm)
{
	free(ptherm->H_window);
	ptherm->H_window = 0;
}
