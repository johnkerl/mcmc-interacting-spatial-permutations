// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef METRO_STATS_H
#define METRO_STATS_H

typedef struct _metro_type_stats_t {
	long long num_self;
	long long num_keep;
	long long num_change;
	long long num_total;
} metro_type_stats_t;

typedef struct _metro_stats_t {
	metro_type_stats_t open_stats;
	metro_type_stats_t close_stats;
	metro_type_stats_t HS_stats;
	metro_type_stats_t TS_stats;
	metro_type_stats_t SO_stats;
} metro_stats_t;

void clear_metro_type_stats (metro_type_stats_t* pmetro_type_stats);
void report_metro_type_stats(metro_type_stats_t* pmetro_type_stats);
void clear_metro_stats (metro_stats_t* pmetro_stats);
void report_metro_stats(metro_stats_t* pmetro_stats);

#endif // METRO_STATS_H
