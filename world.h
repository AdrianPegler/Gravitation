#ifndef WORLD_H
#define WORLD_H

#include <mpi.h>

typedef struct _World  World;

struct _World{
  int  size;
  int  rank;
  char name[MPI_MAX_PROCESSOR_NAME];
};

World world;

#endif //WORLD_H