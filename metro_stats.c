// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include "metro_stats.h"

// ----------------------------------------------------------------
void clear_metro_type_stats(metro_type_stats_t* pmetro_type_stats)
{
	pmetro_type_stats->num_self   = 0;
	pmetro_type_stats->num_keep   = 0;
	pmetro_type_stats->num_change = 0;
	pmetro_type_stats->num_total  = 0;
}

// ----------------------------------------------------------------
void report_metro_type_stats(metro_type_stats_t* pstats)
{
	pstats->num_total = pstats->num_self
		+ pstats->num_keep + pstats->num_change;
	if (pstats->num_total == 0) {
		printf("#   None.\n");
		return;
	}
	if (pstats->num_self != 0) {
		printf("#   Metropolis selves:  %9lld / %9lld (%7.3lf%%)\n",
			pstats->num_self, pstats->num_total,
			(double)pstats->num_self*100.0/(double)pstats->num_total);
	}
	printf("#   Metropolis keeps:   %9lld / %9lld (%7.3lf%%)\n",
		pstats->num_keep, pstats->num_total,
		(double)pstats->num_keep*100.0/(double)pstats->num_total);
	printf("#   Metropolis changes: %9lld / %9lld (%7.3lf%%)\n",
		pstats->num_change, pstats->num_total,
		(double)pstats->num_change*100.0/(double)pstats->num_total);
}

// ----------------------------------------------------------------
void clear_metro_stats(metro_stats_t* pmetro_stats)
{
	clear_metro_type_stats(&pmetro_stats->open_stats);
	clear_metro_type_stats(&pmetro_stats->close_stats);
	clear_metro_type_stats(&pmetro_stats->HS_stats);
	clear_metro_type_stats(&pmetro_stats->TS_stats);
	clear_metro_type_stats(&pmetro_stats->SO_stats);
}

// ----------------------------------------------------------------
void report_metro_stats(metro_stats_t* pmetro_stats)
{
	printf("# OPENS:\n");  report_metro_type_stats(&pmetro_stats->open_stats);
	printf("# CLOSES:\n"); report_metro_type_stats(&pmetro_stats->close_stats);
	printf("# HS:\n");     report_metro_type_stats(&pmetro_stats->HS_stats);
	printf("# TS:\n");     report_metro_type_stats(&pmetro_stats->TS_stats);
	printf("# SO:\n");     report_metro_type_stats(&pmetro_stats->SO_stats);
}
