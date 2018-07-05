#ifndef CLUSTER_H
#define CLUSTER_H

#include "bodies.h"
#include "world.h"
#include "math_helper.h"

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <immintrin.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#define semi_activ -1
#define inactiv -2

typedef struct _Cluster Cluster;

Cluster *activ_sub_tree;

double* ref_points;

int INTERPOLATION_POINTS;
int NUM_SUB_MASSES;
int MAX_DEPTH;
int SPLIT_DEPTH;
int INST;
int split_count;

struct _Cluster{
    bodies *bodies;             //related bodies
    int id;                     //id of the Cluster
    int activ;                  //number of the process that's responsible 
    int start;                  //index of first bodie that resides inside the Cluster
    int n;                      //number of bodies in this Cluster
    double a[3];                //"upper left" corner of the bounding box
    double b[3];                //"lower right" corner of the bounding box
    double center[3];           //center of the bounding box
    double* m;                  //substitution masses
    double* F;                  //substitution forces
    double* xs;                 //coordinates of the locations of submasses
    double diam;                //diameter of the bounding box
    int num_sons;               //number of son clusters
    struct _Cluster *son[2];    //pointers to the son
};

Cluster *constructClusterTree(bodies *b);
void deleteCluster(Cluster *c);

double lagrange(Cluster *c, int i, int dir, double x);

void write_clusters(Cluster *root, char *filename);

#endif //CLUSTER_H