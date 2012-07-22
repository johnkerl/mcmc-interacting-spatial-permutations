// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "util.h"
#include "times.h"

// ----------------------------------------------------------------
void* malloc_or_die(int num_bytes)
{
	void* rv;

	rv = malloc(num_bytes);
	if (rv == 0) {
		fprintf(stderr, "malloc(%d) failed.\n", num_bytes);
		exit(1);
	}

	return rv;
}

// ----------------------------------------------------------------
int* int_malloc_or_die(int num_ints)
{
	return (int*)malloc_or_die(num_ints * sizeof(int));
}

// ----------------------------------------------------------------
char* char_malloc_or_die(int num_chars)
{
	return (char*)malloc_or_die(num_chars * sizeof(char));
}

// ----------------------------------------------------------------
double* double_malloc_or_die(int num_doubles)
{
	return (double*)malloc_or_die(num_doubles * sizeof(double));
}

// ----------------------------------------------------------------
FILE* open_file_or_die(char* file_name, char* mode)
{
	FILE* fp = fopen(file_name, mode);
	if (fp == 0) {
		perror("fopen");
		fprintf(stderr, "Couldn't open file \"%s\" for mode \"%s\".\n",
			file_name, mode);
		exit(1);
	}
	return fp;
}

// ----------------------------------------------------------------
void write_to_file_or_die(void* ptr, size_t size, size_t nmemb, FILE* fp)
{
	if (fwrite(ptr, size, nmemb, fp) != nmemb) {
		perror("fwrite");
		exit(1);
	}
}

// ----------------------------------------------------------------
void read_from_file_or_die(void* ptr, size_t size, size_t nmemb, FILE* fp)
{
	int rv = fread(ptr, size, nmemb, fp);
	if (rv != nmemb) {
		perror("fread");
		fprintf(stderr, "ptr=%p size=%d nmemb=%d rv=%d\n",
			ptr, (int)size, (int)nmemb, (int)rv);
		exit(1);
	}
}

// ----------------------------------------------------------------
void close_file_or_die(FILE* fp)
{
	if (fclose(fp) != 0) {
		perror("fclose");
		exit(1);
	}
}

// ----------------------------------------------------------------
int ipow(int x, int n)
{
	int xn = 1;
	int i;
	for (i = 0; i < n; i++)
		xn *= x;
	return xn;
}

// ----------------------------------------------------------------
// If there is a file in the current directory called "__stop__", then quit.

void check_stop(void)
{
	struct stat statbuf;

	if (stat("./__stop__", &statbuf) == 0) {
		printf("Exiting due to __stop__ file.  Time:  %s",
			get_sys_time_string());
		exit(1);
	}
}
