// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include "interactions.h"
#include "util.h"

// ----------------------------------------------------------------
char* get_interaction_type_desc(int interaction_type)
{
	switch (interaction_type) {

	case NO_INTERACTIONS:   return "none";           break;
	case R2_INTERACTIONS:   return "r2";             break;
	case RELL_INTERACTIONS: return "constant r_ell"; break;
	case VZWL_INTERACTIONS: return "Vzw lower";      break;
	case VZWU_INTERACTIONS: return "Vzw upper";      break;
	default:
		return "(undefined interaction type)";
		break;
	}
}

// ----------------------------------------------------------------
// Other alpha_ell possibilities:
// * r2, r3 only
// * alpha0 * sqrt(ell)
// * ?

double* get_alphas(int interaction_type, int N, double alpha0)
{
	double* alphas = double_malloc_or_die(N+1);
	int k;

	for (k = 0; k <= N; k++)
		alphas[k] = 0.0;

	switch (interaction_type) {

	case NO_INTERACTIONS:
		break;
	case R2_INTERACTIONS:
		alphas[2] = alpha0;
		break;
	case RELL_INTERACTIONS:
		for (k = 1; k <= N; k++)
			alphas[k] = alpha0;
		break;
	case VZWL_INTERACTIONS:
		alphas[2] = alpha0;
		break;
	case VZWU_INTERACTIONS:
		alphas[2] = alpha0;
		break;
	default:
		fprintf(stderr, "get_alphas: undefined interaction type %d.\n",
			interaction_type);
		exit(1);
		break;
	}

	return alphas;
}
