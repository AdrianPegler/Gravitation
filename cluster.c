#include "cluster.h"

/*********************************************************
* DEFINITIONS
*********************************************************/
int *send_count, *recv_count, *send_displ, *recv_displ;
int new_n;
void _setup(Cluster *c, int depth);
Cluster *_new_bound_Cluster(int start, int n, bodies *bodies, double *a, double *b, int active);
void deleteCluster(Cluster *c);

/*********************************************************
* METHODS
*********************************************************/
double lagrange(Cluster *c, int i, int dir, double x){
    double l = 1.0;
    for(int k = 0; k < INTERPOLATION_POINTS; k++){
        if(k != i){
            l *= (x - c->xs[dir * INTERPOLATION_POINTS + k]);
            l /= (c->xs[dir * INTERPOLATION_POINTS + i] - c->xs[dir * INTERPOLATION_POINTS + k]);
        }
    }    
    assert(!isnan(l));
    assert(!isinf(l));
    return l;
}

int _largest_direction(Cluster *c){
    double dist[3];
    for(int i = 0; i < 3; i++){
        dist[i] = _dist_ab(c->a[i], c->b[i]);
    }

    //At least one direction has no extent
    int null = dist[0] == 0.0 || dist[1] == 0.0 || dist[2] == 0.0;
    if(null){
        return -2;
    }

    //Return the widest coordinate direction. Ordered by priority:
    //x-direction:
    int cond = dist[0] >= dist[1] && dist[0] >= dist[2];    
    if(cond){
        return 0;
    }
    //y-direction:
    cond = dist[1] >  dist[0] && dist[1] >= dist[2];
    if(cond){
        return 1;
    }
    //z-direction
    cond = dist[2] >  dist[0] && dist[2] >  dist[1];
    if(cond){
        return 2;
    }    
    //Something went wrong!
    return -1;
}

//returns the array of coordinates of bodies associated with given direction
double *_getAssocCoord(bodies *b, int dir){
    double *assoc_coord;
    switch(dir){
        case 0: assoc_coord = b->x; break;
        case 1: assoc_coord = b->y; break;
        case 2: assoc_coord = b->z; break;
        case -2:
            printf("Error: Cluster hat Ausdehnung 0.\n");
            exit(-1);
        break; 
        default:
            printf("Error: _split_Cluster returned -1.\n");
            exit(-1);
        break;
    }

    return assoc_coord;
}

int _sortIndices(Cluster *c, int dir){
    if(0 == c->n){
        return 0;
    }
    int j, front, back;
    double *assoc_coord;
    
    //get associated coordinates of bodies
    assoc_coord = _getAssocCoord(c->bodies, dir);
    
    //sort the bodies according to clusters center
    front = 0;
    back  = c->n;
    back -= 1;
    do {
        j = c->start + front;
        
        if(assoc_coord[j] < c->center[dir]){
            front++;
        } else {
            if(assoc_coord[c->start + back] <= c->center[dir]){
                swap_bodies(c->bodies, c->start + front, c->start + back);
            }
            back--;
        }
    } while(front <= back);

    return front;
}

void get_bounding_box(bodies *b, double *min, double *max){
    int i;
    double local_min[3], local_max[3];

    local_min[0] = b->x[0];
    local_min[1] = b->y[0];
    local_min[2] = b->z[0];
    local_max[0] = b->x[0];
    local_max[1] = b->y[0];
    local_max[2] = b->z[0];

    for(i = 1; i < b->n; i++){
        local_min[0] = b->x[i] < local_min[0]? b->x[i] : local_min[0];
        local_min[1] = b->y[i] < local_min[1]? b->y[i] : local_min[1];
        local_min[2] = b->z[i] < local_min[2]? b->z[i] : local_min[2];
        local_max[0] = b->x[i] > local_max[0]? b->x[i] : local_max[0];
        local_max[1] = b->y[i] > local_max[1]? b->y[i] : local_max[1];
        local_max[2] = b->z[i] > local_max[2]? b->z[i] : local_max[2];
    }

    MPI_Allreduce(local_min, min, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(local_max, max, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    /* printf("P0%i:\nmin: [%g, %g, %g]\nmax: [%g, %g, %g]\nlocal_min: [%g, %g, %g]\nlocal_max: [%g, %g, %g]\n",
           world.rank, min[0], min[1], min[2], max[0], max[1], max[2],
           local_min[0], local_min[1], local_min[2], local_max[0], local_max[1], local_max[2]); */
}

void _setup_nonLeafCluster(Cluster *c, int depth){
    int dir = _largest_direction(c);
    int border = _sortIndices(c, dir);
    double b1[3],a2[3];
    for(int d = 0; d < 3; d++){
        if(d == dir){
            b1[d] = c->center[d];
            a2[d] = c->center[d];
        } else {
            b1[d] = c->b[d];
            a2[d] = c->a[d];
        }
    } 

    c->num_sons = 2;
    c->son[0] = _new_bound_Cluster(c->start         , border       , c->bodies, c->a, b1  , c->active);
    _setup(c->son[0], depth+1);
    c->son[1] = _new_bound_Cluster(c->start + border, c->n - border, c->bodies, a2  , c->b, c->active);
    _setup(c->son[1], depth+1);
}

void _setup(Cluster *c, int depth){
    if(depth == SPLIT_DEPTH){
        c->active = ++split_count;
    }
    if(depth < MAX_DEPTH){
        _setup_nonLeafCluster(c, depth);
    }

    //till now all clusters are either semiactivee by default
    //or active for exactly one knot. Now all Clusters
    //that are neither really semiactive nor active for this
    //knot will be set to inactive.
    if(c->active == semi_active){
        if(c->son[0]->active == inactive && c->son[1]->active == inactive){
            c->active = inactive;
        } else{
            if(c->son[0]->active != world.rank && c->son[0]->active != semi_active &&
               c->son[1]->active != world.rank && c->son[1]->active != semi_active) {
                   c->active = inactive;
            }
        }
    }
}

void preSort(Cluster *c, int depth){
    if(depth == SPLIT_DEPTH){
        c->active = ++split_count;
    
        //get count and start index of data to send to process #split_count
        send_count[split_count] = c->n;
        send_displ[split_count] = split_count == 0 ? 0 : send_displ[split_count - 1] + send_count[split_count - 1];
    
        //if the current knot is the active one: gather counts and displs from the other knots:
        MPI_Gather(&send_count[split_count], 1, MPI_INT, recv_count, 1, MPI_INT, split_count, MPI_COMM_WORLD);
        if(world.rank == c->active){
            new_n = 0;
            for(int i = 0; i < world.size; i++){
                recv_displ[i] = new_n;
                new_n        += recv_count[i];
            }
        }
        
    } else {
        int dir = _largest_direction(c);
        int border = _sortIndices(c, dir);
        double b1[3],a2[3];
        for(int d = 0; d < 3; d++){
            if(d == dir){
                b1[d] = c->center[d];
                a2[d] = c->center[d];
            } else {
                b1[d] = c->b[d];
                a2[d] = c->a[d];
            }
        } 

        c->num_sons = 2;
        c->son[0] = _new_bound_Cluster(c->start         , border       , c->bodies, c->a, b1  , c->active);
        preSort(c->son[0], depth+1);
        c->son[1] = _new_bound_Cluster(c->start + border, c->n - border, c->bodies, a2  , c->b, c->active);
        preSort(c->son[1], depth+1);
    }
}

/*********************************************************
* CONSTRUCTORS
*********************************************************/

Cluster *constructClusterTree(bodies *b){
    int n = b->n;
    double min[3], max[3];
    bodies *new_bs;

    //get smallest root bounding box
    get_bounding_box(b, min, max);

    Cluster *root = _new_bound_Cluster(0, n, b, min, max, semi_active);

    if(world.size > 1){
        ////////////////////////////////////////////////////////
        // Initializations:
        split_count = -1;
        send_count = _mm_malloc(world.size * sizeof(int), 64);
        recv_count = _mm_malloc(world.size * sizeof(int), 64);
        send_displ = _mm_malloc(world.size * sizeof(int), 64);
        recv_displ = _mm_malloc(world.size * sizeof(int), 64);

        ////////////////////////////////////////////////////////
        // Sort the bodies; first locally, then globally.
        preSort(root, 0);
        new_bs = new_bodies(new_n);
        alltoall_bodies(my_bs, send_count, send_displ, new_bs, recv_count, recv_displ);
        del_bodies(my_bs);
        my_bs = new_bs;
        //reset roots bodies*
        root->bodies = my_bs;
        root->n      = my_bs->n;

        ////////////////////////////////////////////////////////
        // clean preSort
        //delete all but the root of the Cluster tree
        deleteCluster(root->son[0]);
        deleteCluster(root->son[1]);
        //reset counter
        split_count = -1;
        INST = 1;
    } else {
        split_count = -1;
        root->active = 0;
        // active_sub_tree = root;
    }

    ////////////////////////////////////////////////////////
    // construct the rest of the tree
    _setup(root, 0);
    
    ////////////////////////////////////////////////////////
    // clean up:
    _mm_free(send_count);
    _mm_free(recv_count);
    _mm_free(send_displ);
    _mm_free(recv_displ);
    return root;
}

Cluster *_new_bound_Cluster(int start, int n, bodies *bodies, double *a, double *b, int active){
    int i,j;
    Cluster *c;
    c     = (Cluster*) _mm_malloc(                          sizeof(Cluster), 64);
    if(active == world.rank || active == semi_active){   
        c->m  = (double*)  _mm_malloc(NUM_SUB_MASSES           * sizeof(double), 64);
        c->F  = (double*)  _mm_malloc(NUM_SUB_MASSES       * 3 * sizeof(double), 64);
        c->xs = (double*)  _mm_malloc(INTERPOLATION_POINTS * 3 * sizeof(double), 64);
    }else{
        c->m = NULL;
        c->F = NULL;
        c->xs = NULL;
    }
    //c->count = (int*)  _mm_malloc(world.size               * sizeof(int)   , 64);
    //c->displ = (int*)  _mm_malloc(world.size               * sizeof(int)   , 64);

    //init values
    c->bodies = bodies;
    c->id = INST++;
    c->active = active;
    c->start = start;
    c->n = n;
    c->num_sons = 0;
    for(i = 0; i < 3; i++){
        c->a[i] = a[i];
        c->b[i] = b[i];
    }
    /*for(i = 0; i < NUM_SUB_MASSES; i++){
        c->m[i] = 0.0;
        c->F[i * 3 + 0] = 0.0;
        c->F[i * 3 + 1] = 0.0;
        c->F[i * 3 + 2] = 0.0;
    }*/

    //calculate diameter and centerpoint(s)
    double diam, dist;
    diam = 0.0;
    for(i = 0; i < 3; i++){
        c->center[i] = (c->a[i] + c->b[i]) /2.0;
        dist = _dist_ab(c->a[i], c->b[i]);
        if(active == world.rank || active == semi_active){
            for(j = 0; j < INTERPOLATION_POINTS; j++){
                c->xs[i * INTERPOLATION_POINTS + j] = c->center[i] + dist / 2 * ref_points[j];
            }
        }
        dist *= dist;
        diam += dist;
    }
    c->diam = sqrt(diam);
    
    return c;
}

/*********************************************************
* DECONSTRUCTORS
*********************************************************/
void deleteCluster(Cluster *c){
    if(c->num_sons){
        deleteCluster(c->son[0]);
        deleteCluster(c->son[1]);
    }
    if(c->m != NULL){_mm_free(c->m);}
    if(c->F != NULL){_mm_free(c->F);}
    if(c->xs != NULL){_mm_free(c->xs);}
    //_mm_free(c->count);
    //_mm_free(c->displ);
    _mm_free(c);
}

/*********************************************************
* DEBUG
*********************************************************/
void write_cluster(Cluster *c, FILE *file){
    int  i;

    fprintf(file, "id: %d\n", c->id);
    fprintf(file, "a: [ %.3e, %.3e, %.3e ]\n", c->a[0], c->a[1], c->a[2]);
    fprintf(file, "b: [ %.3e, %.3e, %.3e ]\n", c->b[0], c->b[1], c->b[2]);
    fprintf(file, "m: [ %.3e", c->m[0]);
    for(i = 1; i < NUM_SUB_MASSES; i++){
        fprintf(file, "\n     %.3e", c->m[i]);
    }
    fprintf(file, " ]\n");
    fprintf(file, "F: [ [ %.3e, %.3e, %.3e ]", c->F[0], c->F[NUM_SUB_MASSES], c->F[2*NUM_SUB_MASSES]);
    for(i = 1; i < NUM_SUB_MASSES; i++){
        fprintf(file, "\n     [ %.3e, %.3e, %.3e ]", c->F[0 * NUM_SUB_MASSES + i], c->F[1 * NUM_SUB_MASSES + i], c->F[2 * NUM_SUB_MASSES + i]);
    }
    fprintf(file, " ]\n\n");

    for(i = 0; i < c->num_sons; i++){
        write_cluster(c->son[i], file);
    }
}

void write_clusters(Cluster *root, char *filename){
  FILE *file;
  file = fopen(filename, "w");

  write_cluster(root, file);

  fclose(file);
}