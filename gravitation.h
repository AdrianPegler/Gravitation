#ifndef GRAVITATION_H
#define GRAVITATION_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "bodies.h"
#include "world.h"
#include "cluster.h"
#include "stopwatch.h"
#include "eval.h"

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <immintrin.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>

typedef char* String;

/*bodies *recv_buffer;

double F_max[3];

int recv_buffer_max_size;
int body_count;*/

/************************************************************
 * KERNEL
 ***********************************************************/
//void calc_potential(double xi, double yi, double zi, double xj, double yj, double zj, double *P);

#endif //GRAVITATION_H