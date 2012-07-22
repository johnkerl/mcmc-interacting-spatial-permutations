// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include "rcmrand.h"

char* get_rng_type_desc(int rng_type)
{
	char* rv = "Unrecognized RNG type";
	switch (rng_type) {
	case RCM_PSDES:   rv = "psdes";            break;
	case RCM_URANDOM: rv = "/dev/urandom";     break;
	case RCM_RAND48:  rv = "rand48";           break;
	case RCM_MERSTW:  rv = "Mersenne twister"; break;
	}
	return rv;
}
