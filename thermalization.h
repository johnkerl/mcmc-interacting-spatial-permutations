// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef THERMALIZATION_H
#define THERMALIZATION_H

typedef struct _therm_ctl_t {
	int     smoothing_window_size;
	int     therm_ntp;
	int     force_therm;

	int     sweepi;
	double  H_lagger;
	double  ssH; // Scaled, smoothed H
	double  smoothed_H;
	double  smoothed_Hm1;
	double  smoothed_Hm2;
	double  smoothed_H_diff1;
	double  smoothed_H_diff2;
	int     smoothing_scale;
	double* H_window;
	int     ntp;
	int     thermalized;
} therm_ctl_t;

void therm_ctl_init(therm_ctl_t* ptherm, int smoothing_window_size,
	int therm_ntp, int force_therm);
void therm_ctl_update(therm_ctl_t* ptherm, double H);
void therm_ctl_free(therm_ctl_t* ptherm);

#endif // THERMALIZATION_H
