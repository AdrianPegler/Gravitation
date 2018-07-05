#ifndef EVAL_H
#define EVAL_H

#include <stdbool.h>
#include <string.h>
#include <immintrin.h>
#include "cluster.h"
#include "vector.h"
#include "bodies.h"
#include "stopwatch.h"

vector **send_clusters, **recv_clusters, **send_bodies, **recv_bodies;
double *send_sub_ms, *recv_sub_ms;
int *send_sub_ms_count, *send_sub_ms_displ, *recv_sub_ms_count, *recv_sub_ms_displ;
double *send_x, *send_y, *send_z, *send_m, *recv_x, *recv_y, *recv_z, *recv_m;
int *send_bs_count, *send_bs_displ, *recv_bs_count, *recv_bs_displ;
int *send_n, *send_n_count, *send_n_displ, *recv_n, *recv_n_count, *recv_n_displ; 

void calc_potential(double xi, double yi, double zi, double xj, double yj, double zj, double *P);

void init_eval();
void finalize_eval();

void forward(Cluster *c);
void eval(Cluster *ct, Cluster *cs);
void backward(Cluster *c);

#endif //EVAL_H