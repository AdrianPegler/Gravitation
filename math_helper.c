#include "math_helper.h"

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