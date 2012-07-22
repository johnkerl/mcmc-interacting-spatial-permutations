// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>

void *  malloc_or_die      (int num_bytes);
int*    int_malloc_or_die  (int num_ints);
char*   char_malloc_or_die (int num_chars);
double* double_malloc_or_die(int num_doubles);

FILE*   open_file_or_die(char* file_name, char* mode);
void    write_to_file_or_die (void* ptr, size_t size, size_t nmemb, FILE* fp);
void    read_from_file_or_die(void* ptr, size_t size, size_t nmemb, FILE* fp);
void    close_file_or_die(FILE* fp);

int     ipow(int x, int n);

// If there is a file in the current directory called "__stop__", then quit.
void check_stop(void);

#endif // UTIL_H
