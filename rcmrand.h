// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

// ================================================================
// RCMRAND.H
//
// This is a lightweight abstraction layer which allows me to experiment with
// various pseudorandom-number generators.
//
// The generators currently being abstracted here are:
// * pseudo-DES from Numerical Recipes, 2nd ed.
// * rand48() from the standard C library.
//
// The uses to which the generators are put in this software
// application are:
// * SRANDOM():      seed from 32-bit value.
// * URANDOM():      uniform double on [0.0, 1.0).
// * IMODRANDOM(m):  uniform int on {0, 1, 2, ..., m-1}.
// * IRANDOM(m):     32-bit unsigned integer.
// ================================================================
// John Kerl
// kerl.john.r@gmail.com
// 2008-02-02
// ================================================================

#ifndef RCMRAND_H
#define RCMRAND_H

// The pseudorandom number generator(s) being abstracted:
#include <stdlib.h> // rand48

#include <sys/types.h> // For time(0) and getpid()
#include <sys/unistd.h>
#include <time.h>

#include "psdes.h"   // pseudo-DES
#include "urandom.h" // Linux /dev/urandom with caching.
#include "mtrand.h"  // Mersenne twister

#define RCM_PSDES   1
#define RCM_URANDOM 2
#define RCM_RAND48  3
#define RCM_MERSTW  4

// Execution times for two runs with SAR algorithm and L=10:
// psdes:   0.778s 0.855s
// urandom: 6.642s 6.464s
// rand48:  0.736s 0.740s
// Merstw:  0.663s 0.714s

// Here is where one specifies the RNG:
//#define RCM_WHICH RCM_PSDES
//#define RCM_WHICH RCM_URANDOM
//#define RCM_WHICH RCM_RAND48
#define RCM_WHICH RCM_MERSTW

// ----------------------------------------------------------------
// pseudo-DES.
// It has 64-bit state but accepts a 32-bit seed.
#if RCM_WHICH == RCM_PSDES
#define RCM_RAND_DESC "psdes"
#define STRANDOM(s)
#define SRANDOM(s)    sran32((unsigned)(s))
#define URANDOM()     fran32()
#define IMODRANDOM(m) (iran32() % (m))
#define IRANDOM()     (iran32())
#endif

// ----------------------------------------------------------------
// /dev/urandom
#if RCM_WHICH == RCM_URANDOM
#define RCM_RAND_DESC "/dev/urandom"
#define STRANDOM(s)   // No such
#define SRANDOM(s)    // No such
#define URANDOM()     get_urandomu()
#define IMODRANDOM(m) ((int)((unsigned)get_urandom() % (m)))
#define IRANDOM()     ((unsigned)(get_urandom()))
#endif

// ----------------------------------------------------------------
// rand48.
// It has 48-bit state but accepts a 32-bit seed.
#if RCM_WHICH == RCM_RAND48
#define RCM_RAND_DESC "rand48"
#define STRANDOM(s)   srand48((long)(time(0) ^ getpid()))
#define SRANDOM(s)    srand48((long)(s))
#define URANDOM()     drand48()
#define IMODRANDOM(m) ((int)(lrand48() % (m)))
#define IRANDOM()     ((unsigned)lrand48())
#endif

// ----------------------------------------------------------------
// Mersenne twister.
// It has 19,968-bit state but accepts a 32-bit seed.
#if RCM_WHICH == RCM_MERSTW
#define RCM_RAND_DESC "Mersenne twister"
#define STRANDOM(s)   mtrand_init((unsigned)(time(0) ^ getpid()))
#define SRANDOM(s)    mtrand_init((unsigned)(s))
#define URANDOM()     get_mtrand_double()
#define IMODRANDOM(m) ((int)(get_mtrand_int32() % (m)))
#define IRANDOM()     ((unsigned)get_mtrand_int32())
#endif

// ----------------------------------------------------------------
#define RANDRANGE(lo,hi) (lo+(hi-lo)*URANDOM())
#define RANDPM()         (URANDOM() < 0.5 ? 1 : -1)

// ----------------------------------------------------------------
char* get_rng_type_desc(int rng_type);
#endif // RCMRAND_H
