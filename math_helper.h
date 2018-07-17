#ifndef MATH_HELPER_H
#define MATH_HELPER_H

#include <immintrin.h>

int max(int a, int b);
double _dist_ab(double a, double b);
void _swap_i(int *a, int *b);
void _swap_d(double *a, double *b);

void *ap_malloc(size_t size, size_t align);
void ap_free(void *ptr, size_t size);
void *ap_nmalloc(size_t size);
void *ap_realloc(void *ptr, size_t size, size_t oldsize);
void *ap_calloc(size_t nmem, size_t size);
void ap_nfree(void *ptr, size_t size);


#endif //MATH_HELPER_H