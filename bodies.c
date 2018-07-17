#include "bodies.h"
//#include <time.h>

/************************************************************
 * Constructor
 ***********************************************************/

bodies* new_bodies(int n) {
  bodies *b;

  // Allocate memory for the struct.
  b = (bodies*) ap_malloc(sizeof(bodies), 64);

  // Set number of bodies.
  b->n = n;
  b->th = YEAR;

  n = n?n:1;

  b->id = (int*) ap_malloc(n * sizeof(int), 64);

  // Allocate aligned memory for the positions of the masses.
  b->x = (double*) ap_malloc(n * sizeof(double), 64);
  b->y = (double*) ap_malloc(n * sizeof(double), 64);
  b->z = (double*) ap_malloc(n * sizeof(double), 64);

  // Allocate aligned memory for the velocities of the masses.
  b->vx = (double*) ap_malloc(n * sizeof(double), 64);
  b->vy = (double*) ap_malloc(n * sizeof(double), 64);
  b->vz = (double*) ap_malloc(n * sizeof(double), 64);

  // Allocate aligned memory for the Forces of the masses.
  b->Fx = (double*) ap_malloc(n * sizeof(double), 64);
  b->Fy = (double*) ap_malloc(n * sizeof(double), 64);
  b->Fz = (double*) ap_malloc(n * sizeof(double), 64);

  // Allocate aligned memory for the masses.
  b->m = (double*) ap_malloc(n * sizeof(double), 64);

  return b;
}

/************************************************************
 * Deconstructor
 ***********************************************************/

void del_bodies(bodies *b) {
  ap_free(b->id, b->n * sizeof(int));
  ap_free(b->x,  b->n * sizeof(double));
  ap_free(b->y,  b->n * sizeof(double));
  ap_free(b->z,  b->n * sizeof(double));
  ap_free(b->vx, b->n * sizeof(double));
  ap_free(b->vy, b->n * sizeof(double));
  ap_free(b->vz, b->n * sizeof(double));
  ap_free(b->Fx, b->n * sizeof(double));
  ap_free(b->Fy, b->n * sizeof(double));
  ap_free(b->Fz, b->n * sizeof(double));
  ap_free(b->m,  b->n * sizeof(double));
  free(b);
}

/************************************************************
 * Methods
 ***********************************************************/

void get_random_bodies(bodies *b) {
  int n = b->n;
  //double r, a, h;

  int i, start;//, j, k, index, section;

  start = n * world.rank;
  for(i = 0; i < n; i++){
    b->id[i] = start + i;
  }

  srand(/*time(NULL)*/42 * (world.rank + 1));

  //OLD
    /*section = cbrt(n) + 1.0;

    srand(42 * world.rank);

    h = 1.0 / (double) section;

    index = 0;
    for (i = 0; i < section; ++i) {
      for (j = 0; j < section; ++j) {
        for (k = 0; k < section; ++k) {
          // Place masses into cube of length 1 lightyear. 

          r = ((double) rand() / (double) RAND_MAX);
          a = (double) k / (double) section;
          b->x[index] = LY * (1.0 * (a + h * r) - 1.0);

          r = ((double) rand() / (double) RAND_MAX);
          a = (double) j / (double) section;
          b->y[index] = LY * (2.0 * (a + h * r) - 1.0);

          r = ((double) rand() / (double) RAND_MAX);
          a = (double) i / (double) section;
          b->z[index++] = LY * (2.2 * (a + h * r) - 1.0);

          if (index >= n) {
            break;
          }
        }
        if (index >= n) {
          break;
        }
      }
      if (index >= n) {
        break;
      }
    }*/
  //END

  for (i = 0; i < n; ++i) {
    /* position between -1 and +1 LY */
    b->x[i] = 2.0 * LY * ((double) rand() / (double)RAND_MAX) - LY;
    b->y[i] = 2.0 * LY * ((double) rand() / (double)RAND_MAX) - LY;
    b->z[i] = 2.0 * LY * ((double) rand() / (double)RAND_MAX) - LY;

    /* Masses vary between 1 and 15 masses of sun. */
    b->m[i] = SUNMASS * (1.0 + 14.0 * (double) rand() / (double) RAND_MAX);

    /* Set velocities initially to zero. */
    b->vx[i] = 0.0;
    b->vy[i] = 0.0;
    b->vz[i] = 0.0;

    /* Set forces initially to zero. */
    b->Fx[i] = 0.0;
    b->Fy[i] = 0.0;
    b->Fz[i] = 0.0;
  }

}

void copy_body(bodies *b1, bodies *b2, int b1_i, int b2_i){
  b2->id[b2_i] = b1->id[b1_i];

  b2->x[b2_i]  = b1->x[b1_i];
  b2->y[b2_i]  = b1->y[b1_i];
  b2->z[b2_i]  = b1->z[b1_i];

  b2->vx[b2_i] = b1->vx[b1_i];
  b2->vy[b2_i] = b1->vy[b1_i];
  b2->vz[b2_i] = b1->vz[b1_i];

  b2->Fx[b2_i] = b1->Fx[b1_i];
  b2->Fy[b2_i] = b1->Fy[b1_i];
  b2->Fz[b2_i] = b1->Fz[b1_i];

  b2->m[b2_i] = b1->m[b1_i];
}

void copy_bodies(bodies *b1, bodies *b2) {
  int n = b1->n;
  int i;

  assert(b2->n == n);

  for (i = 0; i < n; ++i) {
    copy_body(b1, b2, i, i);
    b2->th = b1->th;
  }
}

void swap_bodies(bodies *b, int i, int j){
  assert(i < b->n && j < b->n);
  _swap_i(b->id + i, b->id + j);
  _swap_d(b->x  + i, b->x  + j);
  _swap_d(b->y  + i, b->y  + j);
  _swap_d(b->z  + i, b->z  + j);
  _swap_d(b->vx + i, b->vx + j);
  _swap_d(b->vy + i, b->vy + j);
  _swap_d(b->vz + i, b->vz + j);
  _swap_d(b->Fx + i, b->Fx + j);
  _swap_d(b->Fy + i, b->Fy + j);
  _swap_d(b->Fz + i, b->Fz + j);
  _swap_d(b->m  + i, b->m  + j);
}

/*void scatter_bodies(bodies *send, bodies *recv, int *counts, int *displ){
  MPI_Scatterv(send->id, counts, displ, MPI_INT   , recv->id, counts[world.rank], MPI_INT   , 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->x , counts, displ, MPI_DOUBLE, recv->x , counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->y , counts, displ, MPI_DOUBLE, recv->y , counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->z , counts, displ, MPI_DOUBLE, recv->z , counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->Fx, counts, displ, MPI_DOUBLE, recv->Fx, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->Fy, counts, displ, MPI_DOUBLE, recv->Fy, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->Fz, counts, displ, MPI_DOUBLE, recv->Fz, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->vx, counts, displ, MPI_DOUBLE, recv->vx, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->vy, counts, displ, MPI_DOUBLE, recv->vy, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->vz, counts, displ, MPI_DOUBLE, recv->vz, counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(send->m , counts, displ, MPI_DOUBLE, recv->m , counts[world.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}*/

void alltoall_bodies(bodies *send, int *send_count, int *send_displ, bodies *recv, int *recv_count, int *recv_displ){
  MPI_Alltoallv(send->id, send_count, send_displ, MPI_INT   , recv->id, recv_count, recv_displ, MPI_INT   , MPI_COMM_WORLD);
  MPI_Alltoallv(send->x , send_count, send_displ, MPI_DOUBLE, recv->x , recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->y , send_count, send_displ, MPI_DOUBLE, recv->y , recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->z , send_count, send_displ, MPI_DOUBLE, recv->z , recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->Fx, send_count, send_displ, MPI_DOUBLE, recv->Fx, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->Fy, send_count, send_displ, MPI_DOUBLE, recv->Fy, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->Fz, send_count, send_displ, MPI_DOUBLE, recv->Fz, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->vx, send_count, send_displ, MPI_DOUBLE, recv->vx, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->vy, send_count, send_displ, MPI_DOUBLE, recv->vy, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->vz, send_count, send_displ, MPI_DOUBLE, recv->vz, recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(send->m , send_count, send_displ, MPI_DOUBLE, recv->m , recv_count, recv_displ, MPI_DOUBLE, MPI_COMM_WORLD);
}

void gather_bodies(bodies *send, int count, bodies *recv, int recv_rank){
  MPI_Gather(send->id, count, MPI_INT   , recv->id, count, MPI_INT   , recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->x , count, MPI_DOUBLE, recv->x , count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->y , count, MPI_DOUBLE, recv->y , count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->z , count, MPI_DOUBLE, recv->z , count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->Fx, count, MPI_DOUBLE, recv->Fx, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->Fy, count, MPI_DOUBLE, recv->Fy, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->Fz, count, MPI_DOUBLE, recv->Fz, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->vx, count, MPI_DOUBLE, recv->vx, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->vy, count, MPI_DOUBLE, recv->vy, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->vz, count, MPI_DOUBLE, recv->vz, count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gather(send->m , count, MPI_DOUBLE, recv->m , count, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
}

void gather_bodies_v(bodies *send, int offset, bodies *recv, int *counts, int *displ, int recv_rank){
  MPI_Gatherv(send->id + offset, counts[world.rank], MPI_INT   , recv->id, counts, displ, MPI_INT   , recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->x  + offset, counts[world.rank], MPI_DOUBLE, recv->x , counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->y  + offset, counts[world.rank], MPI_DOUBLE, recv->y , counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->z  + offset, counts[world.rank], MPI_DOUBLE, recv->z , counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->Fx + offset, counts[world.rank], MPI_DOUBLE, recv->Fx, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->Fy + offset, counts[world.rank], MPI_DOUBLE, recv->Fy, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->Fz + offset, counts[world.rank], MPI_DOUBLE, recv->Fz, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->vx + offset, counts[world.rank], MPI_DOUBLE, recv->vx, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->vy + offset, counts[world.rank], MPI_DOUBLE, recv->vy, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->vz + offset, counts[world.rank], MPI_DOUBLE, recv->vz, counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
  MPI_Gatherv(send->m  + offset, counts[world.rank], MPI_DOUBLE, recv->m , counts, displ, MPI_DOUBLE, recv_rank, MPI_COMM_WORLD);
}

/* Printing function for debug purpose. */
void print_bodies(bodies *b) {
  int n = b->n;

  int i, n2;

  n2 = n > 10 ? 10 : n;

  for (i = 0; i < n2; ++i) {
    printf(
        " m: %.3e   POS: (%+.3e, %+.3e, %+.3e)  v: (%+.3e, %+.3e, %+.3e)  F:(%+.3e, %+.3e, %+.3e)\n",
        b->m[i], b->x[i], b->y[i], b->z[i], b->vx[i], b->vy[i], b->vz[i],
        b->Fx[i], b->Fy[i], b->Fz[i]);
  }
  if (n2 < n) {
    printf("...\n");
  }
}

void write_bodies(bodies *b, char *filename) {
  int n = b->n;

  FILE *file;
  int i;

  file = fopen(filename, "w");

  fprintf(file, "%d\n", n);
  for (i = 0; i < n; ++i) {
    fprintf(file, "%i:\t%.3e\t%+.3e\t%+.3e\t%+.3e\n", b->id[i], b->m[i], b->x[i], b->y[i], b->z[i]);
    fprintf(file, "\tF: [ %.3e, %.3e, %.3e ]\n", b->Fx[i], b->Fy[i], b->Fz[i]);
  }

  fclose(file);
}