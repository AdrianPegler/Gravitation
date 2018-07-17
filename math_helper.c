#include "math_helper.h"

size_t bytes_alloced = 0;
size_t max_alloced = 0;

int max(int a, int b){
  return (a < b)?b:a;
}

double _dist_ab(double a, double b){
    return (a < b)?(b-a):(a-b);
}

void _swap_i(int *a, int *b){
    if(a == b){ return; }
    if(*a == *b){return;}
    *a ^= *b;
    *b ^= *a;
    *a ^= *b;
}

void _swap_d(double *a, double *b){
    if(a == b){ return; }
    double tmp;
    tmp = *a;
    *a  = *b;
    *b  = tmp; 
}

void *ap_malloc(size_t size, size_t align){
    bytes_alloced += size;
    max_alloced = bytes_alloced>max_alloced?bytes_alloced:max_alloced;
    return _mm_malloc(size, align);
}

void ap_free(void *ptr, size_t size){
    bytes_alloced -= size;
    _mm_free(ptr);
}