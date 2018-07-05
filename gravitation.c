#include "gravitation.h"
#include <locale.h>
int pot;
int steps;

/************************************************************
 * Methods
 ***********************************************************/
  void calc_ref_points(){
      if(NULL == ref_points){
          ref_points = (double*) _mm_malloc(INTERPOLATION_POINTS * sizeof(double), 64);    
      }
      
      int q = INTERPOLATION_POINTS;
      double p;
      if(q == 1){
          ref_points[0] = 0.0;
      }else{
          for(int i = 0; i < q; i++){
              p = (double) (2*i + 1);
              p /= (double) 2*INTERPOLATION_POINTS;
              p = cos(p * M_PI);
              ref_points[i] = p;
          }
      }
  }

/*
* Initializes program variables from main arguments.
* Calls MPI_Init and stores world members.
*/
  void init(int argc, String* argv){
    INTERPOLATION_POINTS = 3;
    pot = 15;
    steps = 1;
    MPI_Init(0, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world.rank);

    switch(argc){
      case 4:
        steps = atoi(argv[3]);
      case 3:
        pot = atoi(argv[2]);
      case 2:
        INTERPOLATION_POINTS = atoi(argv[1]);
      case 1:
      case 0:
        break;
      default:
        printf("WARNING: Too many arguments!\n");
        break;
    }
    NUM_SUB_MASSES = INTERPOLATION_POINTS * INTERPOLATION_POINTS * INTERPOLATION_POINTS;

    SPLIT_DEPTH = (int) (log2(world.size) + 0.5);
    int leaf_pot = INTERPOLATION_POINTS - 1 ? ceil(log2(2 * NUM_SUB_MASSES)) : 4;
    if(leaf_pot > pot){
      leaf_pot = pot;
    }
    MAX_DEPTH = SPLIT_DEPTH + pot - leaf_pot;
    calc_ref_points();
    init_eval();
  }

  void finalize(){
    //del_bodies(recv_buffer);
    finalize_eval();
    _mm_free(ref_points);
    MPI_Finalize();
  }

/************************************************************
 * DEBUG
 ***********************************************************/
  bool all_there(bodies *b){
    int i, j;
    for(i = 0; i < b->n; i++){
      for(j = 0; j < b->n; j++){
        if(b->id[j] == i){
          break;
        }
      }
      if(j >= b->n){
        return false;
      }
    }
    return true;
  }

/************************************************************
 * ERROR
 ***********************************************************/
  double max_error[3];
  double max_force[3];

  bodies *global_bodies;

  void calc_error(){
    int *count, *displ, sum, i, j;
    count = malloc(world.size * sizeof(int));
    displ = malloc(world.size * sizeof(int));

    MPI_Allgather(&(my_bs->n), 1, MPI_INT, count, 1, MPI_INT, MPI_COMM_WORLD);
    sum = 0;
    for(i = 0; i < world.size; i++){
      displ[i] = sum;
      sum += count[i];
    }

    gather_bodies_v(my_bs, 0, global_bodies, count, displ, 0);

    world.rank?:printf("All there: %s\n", all_there(global_bodies)?"true":"false");

    /*******************
     * calc errors */
    if(!world.rank){
      double F[3], P[3];
      double *x, *y, *z, *m;
      double error;
      x = global_bodies->x;
      y = global_bodies->y;
      z = global_bodies->z;
      m = global_bodies->m;

      for(i = 0; i < sum; i++){
        F[0] = 0.0;
        F[1] = 0.0;
        F[2] = 0.0;
        for(j = 0; j < sum; j++){
          if(i == j){continue;}
          calc_potential(x[i], y[i], z[i], x[j], y[j], z[j], P);
          F[0] += GAMMA * m[i] * m[j] * P[0];
          F[1] += GAMMA * m[i] * m[j] * P[1];
          F[2] += GAMMA * m[i] * m[j] * P[2];
        }

        error = fabs(F[0] - global_bodies->Fx[i]);
        max_error[0] = error > max_error[0] ? error : max_error[0];
        max_force[0] = fabs(F[0]) > max_force[0] ? fabs(F[0]) : max_force[0];

        error = fabs(F[1] - global_bodies->Fy[i]);
        max_error[1] = error > max_error[1] ? error : max_error[1];
        max_force[1] = fabs(F[1]) > max_force[1] ? fabs(F[1]) : max_force[1];

        error = fabs(F[2] - global_bodies->Fz[i]);
        max_error[2] = error > max_error[2] ? error : max_error[2];
        max_force[2] = fabs(F[2]) > max_force[2] ? fabs(F[2]) : max_force[2];
      }

      world.rank?:printf("error:\n");
      world.rank?:printf("\t%g\n", (max_error[0] / max_force[0] + max_error[1] / max_force[1] + max_error[2] / max_force[2])/3);
    }

    free(count);
    free(displ);
  }

/************************************************************
 * Main
 ***********************************************************/

int main(int argc, String *argv){
  int n, i, j;
  //pstopwatch sw_1, sw_2;
  Cluster *c;
  double *work_times, *prep_times, *comm_times;
  double avg_work_time, min_work_time, max_work_time;
  double avg_prep_time, min_prep_time, max_prep_time;
  double avg_comm_time, min_comm_time, max_comm_time;

  init(argc, argv);

  /* define the number of interacting particles. */
  n = 2 << (pot-1);

  my_bs     = new_bodies(n);
  //sw_1 = world.rank?NULL:new_stopwatch();
  //sw_2 = world.rank?NULL:new_stopwatch();

  get_random_bodies(my_bs);

  // total_time_start = MPI_Wtime();
  c = constructClusterTree(my_bs);

  work_times = malloc(world.size * sizeof(double));
  prep_times = malloc(world.size * sizeof(double));
  comm_times = malloc(world.size * sizeof(double));

  // char *local_save = setlocale(LC_ALL, NULL);
  // char *local      = 
  setlocale(LC_ALL, "de_DE.UTF-8");

  //printf("old_local:%s \t new_local: %s\n\n", local_save, local);

  world.rank?:printf("Run with %d processes and 2^%d elements per node\r\n", world.size, pot);
  world.rank?:printf("step;total_time;work_time_min;work_time_avg;work_time_max;prep_time_min;prep_time_avg;prep_time_max;comm_time_min;comm_time_avg;comm_time_max\r\n");

  for(i = 1; i <= steps; i++){
    //world.rank?:printf("step %i / %i: ...\r", i, steps);
    total_time_start = MPI_Wtime();

    forward(c);
    eval(c, c);
    backward(c);

    total_time = MPI_Wtime() - total_time_start;
    // world.rank?:printf("step %i / %i:     \n", i , steps);
    // world.rank?:printf("\ttime: %.3f s\n", time);

    //start exact timing
      MPI_Gather(&work_time, 1, MPI_DOUBLE, work_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&prep_time, 1, MPI_DOUBLE, prep_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&comm_time, 1, MPI_DOUBLE, comm_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      min_work_time = work_times[0];
      avg_work_time = work_times[0];
      max_work_time = work_times[0];
      min_prep_time = prep_times[0];
      avg_prep_time = prep_times[0];
      max_prep_time = prep_times[0];
      min_comm_time = comm_times[0];
      avg_comm_time = comm_times[0];
      max_comm_time = comm_times[0];
      for(j = 1; j < world.size; j++){
        min_work_time = work_times[j] < min_work_time ? work_times[j] : min_work_time;
        avg_work_time += work_times[j];
        max_work_time = work_times[j] > max_work_time ? work_times[j] : max_work_time;
        min_prep_time = prep_times[j] < min_prep_time ? prep_times[j] : min_prep_time;
        avg_prep_time += prep_times[j];
        max_prep_time = prep_times[j] > max_prep_time ? prep_times[j] : max_prep_time;
        min_comm_time = comm_times[j] < min_comm_time ? comm_times[j] : min_comm_time;
        avg_comm_time += comm_times[j];
        max_comm_time = comm_times[j] > max_comm_time ? comm_times[j] : max_comm_time;
      }
      avg_work_time /= world.size;
      avg_prep_time /= world.size;
      avg_comm_time /= world.size;
      work_time = 0.0;
      prep_time = 0.0;
      comm_time = 0.0;

      // world.rank?:printf("\t work_time: %.3g \t %.3g\t\t %.3g \n", min_work_time, avg_work_time, max_work_time);
      // world.rank?:printf("\t prep_time: %.3g \t %.3g \t %.3g \n", min_prep_time, avg_prep_time, max_prep_time);
      // world.rank?:printf("\t comm_time: %.3g \t %.3g \t %.3g \n", min_comm_time, avg_comm_time, max_comm_time);
      world.rank?:printf("%d;%.3g;%.3g;%.3g;%.3g;%.3g;%.3g;%.3g;%.3g;%.3g;%.3g\r\n", i, total_time, min_work_time, avg_work_time, max_work_time,
                                                                                                    min_prep_time, avg_prep_time, max_prep_time,
                                                                                                    min_comm_time, avg_comm_time, max_comm_time);

    //end exact timing
  }

  // total_time = MPI_Wtime() - total_time_start;
  //world.rank?:printf("\ntime over all: %.3f s\n", time);

  /*****************
   * mathematical error:*/
  global_bodies = new_bodies(world.rank?0:world.size * n);
  calc_error();
  del_bodies(global_bodies);

  deleteCluster(c);
  del_bodies(my_bs);
  free(work_times);
  free(prep_times);
  free(comm_times);
  // world.rank?:del_stopwatch(sw_1);
  // world.rank?:del_stopwatch(sw_2);
  finalize();
  
}