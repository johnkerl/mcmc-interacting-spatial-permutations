// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#define NO_INTERACTIONS   11
#define R2_INTERACTIONS   22
#define RELL_INTERACTIONS 33
#define VZWL_INTERACTIONS 44
#define VZWU_INTERACTIONS 55

char* get_interaction_type_desc(int interaction_type);

// The caller must free the return value.
double* get_alphas(int interaction_type, int N, double alpha0);

#endif // INTERACTIONS_H
