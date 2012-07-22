// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include "vector_misc.h"

// ----------------------------------------------------------------
void vector_accumulate(double * sum, double * v, int length)
{
	int i;
	for (i = 0; i < length; i++)
		sum[i] += v[i];
}

// ----------------------------------------------------------------
void vector_iaccumulate(double * sum, int * v, int length)
{
	int i;
	for (i = 0; i < length; i++)
		sum[i] += (double)v[i];
}

// ----------------------------------------------------------------
void vector_scale(double* v, double* u, double a, int length)
{
	int i;
	for (i = 0; i < length; i++)
		v[i] = a * u[i];
}

// ----------------------------------------------------------------
void compute_vector_delta(double* v, double* u, int length)
{
	int i;
	int lm1 = length - 1;
	for (i = 0; i < lm1; i++)
		v[i] = u[i+1] - u[i];
	v[length-1] = 0.0;
}

// ----------------------------------------------------------------
void vector_row_printf(double* v, int length, char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		if (i > 0)
			printf(" ");
		printf(format, v[i]);
	}
	printf("\n");
}

// ----------------------------------------------------------------
void vector_col_printf (double* v, int length, char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		printf(format, v[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
void vector_col2_printf(double* v1, double* v2, int length, char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
void vector_col3_printf(double* v1, double* v2, double* v3, int length,
	char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
void vector_col4_printf(double* v1, double* v2, double* v3, double* v4,
	int length, char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf(" ");
		printf(format, v4[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
void vector_col5_printf(double* v1, double* v2, double* v3, double* v4, double* v5,
	int length, char* format)
{
	int i;
	for (i = 0; i < length; i++) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf(" ");
		printf(format, v4[i]);
		printf(" ");
		printf(format, v5[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
// Plot at most 2000 points.
void vector_col5_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5,
	int length, char* format)
{
	int stride = 1;
	int i;
	const int maxlen = 2000;

	if (length > maxlen)
		stride = length / maxlen;

	for (i = 0; i < length; i += stride) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf(" ");
		printf(format, v4[i]);
		printf(" ");
		printf(format, v5[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
// Plot at most 2000 points.
void vector_col6_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5, double* v6,
	int length, char* format)
{
	int stride = 1;
	int i;
	const int maxlen = 2000;

	if (length > maxlen)
		stride = length / maxlen;

	for (i = 0; i < length; i += stride) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf(" ");
		printf(format, v4[i]);
		printf(" ");
		printf(format, v5[i]);
		printf(" ");
		printf(format, v6[i]);
		printf("\n");
	}
}

// ----------------------------------------------------------------
// Plot at most 2000 points.
void vector_col7_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5, double* v6,
	double* v7, int length, char* format)
{
	int stride = 1;
	int i;
	const int maxlen = 2000;

	if (length > maxlen)
		stride = length / maxlen;

	for (i = 0; i < length; i += stride) {
		printf(format, v1[i]);
		printf(" ");
		printf(format, v2[i]);
		printf(" ");
		printf(format, v3[i]);
		printf(" ");
		printf(format, v4[i]);
		printf(" ");
		printf(format, v5[i]);
		printf(" ");
		printf(format, v6[i]);
		printf(" ");
		printf(format, v7[i]);
		printf("\n");
	}
}
