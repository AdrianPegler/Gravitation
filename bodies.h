#ifndef BODIES_H
#define BODIES_H

/* define some astrophysical constants. */
#define LY 9.461e15
#define GAMMA 6.67384e-11 // gravitational constant
#define SUNMASS 1.98892e30

#define MINUTE 60.0;
#define HOUR 60.0 * MINUTE;
#define DAY 24 * HOUR;
#define WEEK 7 * DAY;
#define YEAR 365.25 * DAY;

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <immintrin.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#include "math_helper.h"
#include "world.h"


typedef struct _bodies bodies;

bodies *my_bs;

/* Struct that holds coordinates x, masses m and forces F for all particles. */
struct _bodies {
  int    *id;
  double *x;  // x-components of the masses.
  double *y;  // y-components of the masses.
  double *z;  // z-components of the masses.
  double *vx; // x-components of the velocities.
  double *vy; // y-components of the velocities.
  double *vz; // z-components of the velocities.
  double *Fx; // x-components of the forces.
  double *Fy; // y-components of the forces.
  double *Fz; // z-components of the forces.
  double *m;  // The masses.
  int     n;  // Number of masses.
  double th;  // length of one timestep.
};

bodies* new_bodies(int n);
void del_bodies(bodies *b);

void get_random_bodies(bodies *b);

void copy_body(bodies *b1, bodies *b2, int b1_i, int b2_i);
void copy_bodies(bodies *b1, bodies *b2) ;
void swap_bodies(bodies *b, int i, int j);

//void scatter_bodies(bodies *send, bodies *recv, int *counts, int *displ);
void alltoall_bodies(bodies *send, int *send_count, int *send_displ, bodies *recv, int *recv_count, int *recv_displ);
void gather_bodies(bodies *send, int count, bodies *recv, int recv_rank);
void gather_bodies_v(bodies *send, int offset, bodies *recv, int *counts, int *displ, int recv_rank);

void print_bodies(bodies *b);
void write_bodies(bodies *b, char *filename);

#endif //BODIES_H