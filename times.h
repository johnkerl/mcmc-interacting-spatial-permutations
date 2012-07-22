// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef TIMES_H
#define TIMES_H

// This takes the timeval struct from gettimeofday() and puts it into a double.
double get_sys_time_double(void);

// This is returns ctime with gettimeofday's tv_sec.
// It's a printable string containing the time of day.
char*  get_sys_time_string(void);

#endif // TIMES_H
