// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "times.h"

// ----------------------------------------------------------------
double get_sys_time_double(void)
{
	struct timeval now;
	if (gettimeofday(&now, 0) < 0) {
		perror("gettimeofday");
		exit(1);
	}
	return (double)now.tv_sec + (double)now.tv_usec * 1e-6;
}

// ----------------------------------------------------------------
char* get_sys_time_string(void)
{
	struct timeval now;
	if (gettimeofday(&now, 0) < 0) {
		perror("gettimeofday");
		exit(1);
	}
	return ctime(&now.tv_sec);
}

