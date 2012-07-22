// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef VECTOR_MISC_H
#define VECTOR_MISC_H

// for i = 0 to length-1:
//   sum[i] += v[i]
void vector_accumulate(double* sum, double* v, int length);
void vector_iaccumulate(double* sum, int* v, int length);

// for i = 0 to length-1:
//   v[i] = a * u[i]
void vector_scale(double* v, double* u, double a, int length);

// for i = 0 to length-2:
//   v[i] = u[i+1] - u[i]
// v[length-1] = 0
void compute_vector_delta(double* v, double* u, int length);

void vector_row_printf(double* v, int length, char* format);

void vector_col_printf (double* v,  int length, char* format);
void vector_col2_printf(double* v1, double* v2, int length, char* format);
void vector_col3_printf(double* v1, double* v2, double* v3, int length,
	char* format);
void vector_col4_printf(double* v1, double* v2, double* v3, double* v4,
	int length, char* format);
void vector_col5_printf(double* v1, double* v2, double* v3, double* v4, double* v5,
	int length, char* format);
void vector_col5_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5,
	int length, char* format);
void vector_col6_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5, double* v6,
	int length, char* format);
void vector_col7_printf_decimated(
	double* v1, double* v2, double* v3, double* v4, double* v5, double* v6,
	double* v7, int length, char* format);

#endif // VECTOR_MISC_H
