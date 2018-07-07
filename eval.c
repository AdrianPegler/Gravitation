#include "eval.h"

double ETA = 1.0;
int first_admissable = 1000;
int current_level = 0;


/************************************************************
 * KERNEL
 ***********************************************************/
    void calc_potential(double xi, double yi, double zi, double xj, double yj, double zj, double *P){
        double d[3], norm, norm3;
        d[0] = xj - xi;
        d[1] = yj - yi;
        d[2] = zj - zi;
        norm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        norm = 1.0 / norm;
        norm3 = norm * norm * norm;
        P[0] = d[0] * norm3;
        P[1] = d[1] * norm3;
        P[2] = d[2] * norm3;
    }

/************************************************************
 * Methods
 ***********************************************************/
    double _dist_II(double a0, double b0, double a1, double b1){
        return (b0 < a1) ? (a1 - b0) : (b1 < a0) ? (a0 - b1) : 0.0;
    }

    double _dist_CC(Cluster *c1, Cluster *c2){
        double x_dist = _dist_II(c1->a[0], c1->b[0], c2->a[0], c2->b[0]);
        double y_dist = _dist_II(c1->a[1], c1->b[1], c2->a[1], c2->b[1]);
        double z_dist = _dist_II(c1->a[2], c1->b[2], c2->a[2], c2->b[2]);
        double dist = x_dist * x_dist + y_dist * y_dist + z_dist * z_dist;
        dist = sqrt(dist);
        return dist;
    }

    bool admissable(Cluster *c1, Cluster *c2){
        return fmax(c1->diam, c2->diam) < 2.0*_dist_CC(c1,c2) * ETA;
    }

    void add(vector **vector_array, Cluster *c_recv, Cluster *c_add){
        int active = c_recv->active;

        if(active >= 0){
            vector_add_once(vector_array[active], c_add);
        }else{
            add(vector_array, c_recv->son[0], c_add);
            add(vector_array, c_recv->son[1], c_add);
        } 
    }

/************************************************************
 * Forward
 ***********************************************************/

    void _forward_leaf(Cluster *c){
        if(0 == c->n){
            return;
        }
        
        int i, j, k, l, m, w, wj;
        double x, y, z, *mass, l1, l2, l3;
        bodies *b = c->bodies;
        mass = b->m;

        for(i = 0; i < NUM_SUB_MASSES; i++){
            c->m[i] = 0.0;
        }

        //for all linked actual masses
        for(l = 0; l < c->n; l++){
            m = c->start + l;
            x = b->x[m];
            y = b->y[m];
            z = b->z[m];
            
            //for all interpolation points:
            for(i = 0; i < INTERPOLATION_POINTS; i++){
                l1 = lagrange(c, i, 0, x);
                for(j = 0; j < INTERPOLATION_POINTS; j++){
                    wj = i * INTERPOLATION_POINTS + j;
                    l2 = l1 * lagrange(c, j, 1, y);
                    for(k = 0; k < INTERPOLATION_POINTS; k++){
                        //actual interpolation point:
                        w = wj * INTERPOLATION_POINTS + k;
                        //actual lagrange weight:
                        l3 = l2 * lagrange(c, k, 2, z);

                        //calculate the substitution mass
                        c->m[w] += mass[m] * l3;
                    }
                }
            }
        }
    }

    //calculate substitution mass by adding masses of son clusters.
    void _forward_nonLeaf(Cluster *c){
        int i, j, k, w, wj, son, ison, json, kson, wson, wsonj;
        double sum, l1, l2, l3;

        forward(c->son[0]);
        forward(c->son[1]);

        //for all submasses
        for(i = 0; i < INTERPOLATION_POINTS; i++){
            for(j = 0; j < INTERPOLATION_POINTS; j++){
                wj = i * INTERPOLATION_POINTS + j;
                for(k = 0; k < INTERPOLATION_POINTS; k++){
                    w = wj * INTERPOLATION_POINTS + k;
                    sum = 0.0;
                    //for all sons
                    for(son = 0; son < 2; son++){
                        //for all sons submasses
                        for(ison = 0; ison < INTERPOLATION_POINTS; ison++){
                            l1 = lagrange(c, i, 0, c->son[son]->xs[ison]);
                            for(json = 0; json < INTERPOLATION_POINTS; json++){
                                wsonj = ison * INTERPOLATION_POINTS + json;
                                l2 = l1 * lagrange(c, j, 1, c->son[son]->xs[INTERPOLATION_POINTS + json]);
                                for(kson = 0; kson < INTERPOLATION_POINTS; kson++){
                                    wson = wsonj * INTERPOLATION_POINTS + kson;
                                    l3 = l2 * lagrange(c, k, 2, c->son[son]->xs[2 * INTERPOLATION_POINTS + kson]);

                                    //calculate the substitution mass
                                    sum += c->son[son]->m[wson] * l3;
                                }
                            }
                        }
                    }
                    assert(!isinf(sum));
                    assert(!isnan(sum));
                    c->m[w] = sum;
                }
            }
        }
    }

/************************************************************
 * Communication
 ***********************************************************/

    /** Structure is similar to the later evaluation. The method walks through the block tree and keeps track of
     * all admissable and inadmissable blocks and stores a pointer to the corredsponding target and source clusters
     * for later to exchange the data / respectivly fill in the received data. 
     * */
    void _prep_comm(Cluster *ct, Cluster *cs){
        if(ct->active != semi_active && ct->active != world.rank){
            return;
        }
        if(cs->active == world.rank){
            return;
        }

        //In case of an admissable block with inactive cs:
        //1. This knot needs the information of the knot where cs is active
        //2. Because of symmetry this knot needs to send its own information to the knot where cs is active
        if(admissable(ct, cs)){
            add(recv_clusters, cs, cs);
            add(send_clusters, cs, ct);            
        } else {
            if (ct->num_sons && cs->num_sons){//both clusters have sons left:
                if(ct->id <= cs->id){
                    _prep_comm(ct->son[0], cs->son[0]);
                    _prep_comm(ct->son[0], cs->son[1]);
                    _prep_comm(ct->son[1], cs->son[0]);
                    _prep_comm(ct->son[1], cs->son[1]);
                } else {
                    _prep_comm(ct->son[0], cs->son[0]);
                    _prep_comm(ct->son[1], cs->son[0]);
                    _prep_comm(ct->son[0], cs->son[1]);
                    _prep_comm(ct->son[1], cs->son[1]);
                }
            }
            else if (ct->num_sons){         //just target cluster has sons left:
                printf("WANRING: target cluster has sons left while source cluster has not!\n");
                _prep_comm(ct->son[0], cs);
                _prep_comm(ct->son[1], cs);
            } else if (cs->num_sons) {        //just source cluster has son left:
                printf("WANRING: source cluster has sons left while target cluster has not!\n");
                _prep_comm(ct, cs->son[0]);
                _prep_comm(ct, cs->son[1]);
            }
            else {                          //no son clusters left but still not admissable:
        // Similar to the admissable case the knots need to exchange information.
                assert(!ct->num_sons && !cs->num_sons);
                add(recv_bodies, cs, cs);
                add(send_bodies, cs, ct);
            }
        }
    }

    /** This method mainly (re)allocates memory for the send and recv buffers and copies the corresponding data to the send buffers.
     * */
    void _prep_buffers(){
        int i, j, n_sum, sum, knot_n, offset;
        vector *v;
        Cluster *c;

        //setup for copying all substitution masses to send buffer
        sum = 0;
        for(i = 0; i < world.size; i++){
            v = send_clusters[i];
            send_sub_ms_displ[i] = sum;
            send_sub_ms_count[i] = v->count * NUM_SUB_MASSES;
            sum += send_sub_ms_count[i];
        }
        send_sub_ms = realloc(send_sub_ms, sum * sizeof(double));
        sum = 0;
        for(i = 0; i < world.size; i++){
            v = recv_clusters[i];
            recv_sub_ms_displ[i] = sum;
            recv_sub_ms_count[i] = v->count * NUM_SUB_MASSES;
            sum += recv_sub_ms_count[i];
        }
        recv_sub_ms = realloc(recv_sub_ms, sum * sizeof(double));

        //actually copy substitution masses to send buffer
        offset = 0;
        for(i = 0; i < world.size; i++){
            v = send_clusters[i];
            for(j = 0; j < v->count; j++){
                c = (Cluster*)v->data[j];
                memcpy(send_sub_ms + offset, c->m, NUM_SUB_MASSES * sizeof(double));
                offset += NUM_SUB_MASSES;
            }
        }

        //setup for copying needed bodies data
        sum = 0;
        n_sum = 0;
        for(i = 0; i < world.size; i++){
            v = send_bodies[i];
            knot_n = 0;
            for(j = 0; j < v->count; j++){
                c = (Cluster*) v->data[j];
                knot_n += c->n;
            }
            send_bs_count[i] = knot_n;
            send_bs_displ[i] = sum;
            sum += knot_n;

            send_n_count[i] = v->count;
            send_n_displ[i] = n_sum;
            n_sum += v->count;
        }
        send_x = realloc(send_x, sum * sizeof(double));
        send_y = realloc(send_y, sum * sizeof(double));
        send_z = realloc(send_z, sum * sizeof(double));
        send_m = realloc(send_m, sum * sizeof(double));
        send_n = realloc(send_n, n_sum * sizeof(int));
        // Only send counts can be generated locally, so this also needs to be exchanged
        MPI_Alltoall(send_bs_count, 1, MPI_INT, recv_bs_count, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(send_n_count, 1, MPI_INT, recv_n_count, 1, MPI_INT, MPI_COMM_WORLD);
        sum = 0;
        n_sum = 0;
        for(i = 0; i < world.size; i++){
            recv_bs_displ[i] = sum;
            sum += recv_bs_count[i];

            recv_n_displ[i] = n_sum;
            n_sum += recv_n_count[i];
        }
        recv_x = realloc(recv_x, sum * sizeof(double));
        recv_y = realloc(recv_y, sum * sizeof(double));
        recv_z = realloc(recv_z, sum * sizeof(double));
        recv_m = realloc(recv_m, sum * sizeof(double));
        recv_n = realloc(recv_n, n_sum * sizeof(int));

        //actually get the needed bodies data
        n_sum = 0;
        for(i = 0; i < world.size; i++){
            v = send_bodies[i];
            sum  = 0;
            for(j = 0; j < v->count; j++){
                c = (Cluster*) v->data[j];
                offset = send_bs_displ[i] + sum;
                sum += c->n;
                memcpy(send_x + offset, my_bs->x + c->start, c->n * sizeof(double));
                memcpy(send_y + offset, my_bs->y + c->start, c->n * sizeof(double));
                memcpy(send_z + offset, my_bs->z + c->start, c->n * sizeof(double));
                memcpy(send_m + offset, my_bs->m + c->start, c->n * sizeof(double));
                send_n[n_sum] = c->n;
                n_sum++;
            }
        }
    }

    /** As the name implies here the mein communication takes place. The MPI_Alltoall() can be seen as a
     * matrix transposition so in the end each knot has all the data needed to evaluate the block tree locally.
     * */
    void _communicate(){
        //communicate all Clusterinformation, a.k. the substitution masses
        MPI_Alltoallv(send_sub_ms, send_sub_ms_count, send_sub_ms_displ, MPI_DOUBLE,
                      recv_sub_ms, recv_sub_ms_count, recv_sub_ms_displ, MPI_DOUBLE, MPI_COMM_WORLD);

        //communicate the needed information of the bodies of the inadmissable blocks
        MPI_Alltoallv(send_x, send_bs_count, send_bs_displ, MPI_DOUBLE,
                      recv_x, recv_bs_count, recv_bs_displ, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(send_y, send_bs_count, send_bs_displ, MPI_DOUBLE,
                      recv_y, recv_bs_count, recv_bs_displ, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(send_z, send_bs_count, send_bs_displ, MPI_DOUBLE,
                      recv_z, recv_bs_count, recv_bs_displ, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(send_m, send_bs_count, send_bs_displ, MPI_DOUBLE,
                      recv_m, recv_bs_count, recv_bs_displ, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(send_n, send_n_count , send_n_displ , MPI_INT   ,
                      recv_n, recv_n_count , recv_n_displ , MPI_INT   , MPI_COMM_WORLD);
    }

    /** This method mainly copies the recieved data to their destination. In case of the sent "bodies" the method
     * will only store the position and count of the data within the recv buffers, this saves space (aspecially as
     * the Cluster would need to store an whole bodies structure). The evaluation then only needs to retrieve the
     * corresponding data from the recv buffers.
     * */
    void _finalize_comm(){
        int i, j, k, offset;
        Cluster *c;
        vector *v;
        
        //write the received substitution masses to the right Cluster
        offset = 0;
        for(i = 0; i < world.size; i++){
            v = recv_clusters[i];
            for(j = 0; j < v->count; j++){
                c = (Cluster*) v->data[j];
                if(c->active > -1){ 
                    //All data come from exactly one other knot. So a memcpy is sufficiant
                    memcpy(c->m, recv_sub_ms + offset, NUM_SUB_MASSES * sizeof(double));
                } else {
                    //data may come from different locations so the data needs to be summed up
                    for(k = 0; k < NUM_SUB_MASSES; k++){
                        c->m[k] += recv_sub_ms[offset + k];
                    }
                }
                offset += NUM_SUB_MASSES;
            }
        }

        //change cluster data according to received bodies data
        offset = 0;
        k = 0;
        for(i = 0; i < world.size; i++){
            v = recv_bodies[i];
            for(j = 0; j < v->count; j++){
                c = (Cluster*) v->data[j];
                c->start = offset;
                c->n = recv_n[k];
                offset += recv_n[k];
                k++;
            }
        }
    }

    void _clear_inactive_clusters(){
        int i, j, k;
        Cluster *c;
        vector *v;
        for(i = 0; i < world.size; i++){
            v = recv_clusters[i];
            for(j = 0; j < v->count; j++){
                c = (Cluster*) v->data[j];
                if(c->active == inactive){ 
                    for(k = 0; k < NUM_SUB_MASSES; k++){
                        c->m[k] = 0.0;
                    }
                }
            }
        }
    }

/************************************************************
 * Eval
 ***********************************************************/
    void _eval_CC(Cluster *ct, Cluster *cs){
        int i1, i2, j1, j2, k1, k2, w1, w1j, w2, w2j;
        double F[3];
        // for all interpolation points of target Cluster
        for(i1 = 0; i1 < INTERPOLATION_POINTS; i1++){
            for(j1 = 0; j1 < INTERPOLATION_POINTS; j1++){
                w1j = i1 * INTERPOLATION_POINTS + j1;
                for(k1 = 0; k1 < INTERPOLATION_POINTS; k1++){
                    w1 = w1j * INTERPOLATION_POINTS + k1;
                    //for all interpolations points of source Cluster
                    for(i2 = 0; i2 < INTERPOLATION_POINTS; i2++){
                        for(j2 = 0; j2 < INTERPOLATION_POINTS; j2++){
                            w2j = i2 * INTERPOLATION_POINTS + j2;
                            for(k2 = 0; k2 < INTERPOLATION_POINTS; k2++){
                                w2 = w2j * INTERPOLATION_POINTS + k2;
                                //calculate substitution forces
                                calc_potential(ct->xs[0 * INTERPOLATION_POINTS + i1], 
                                               ct->xs[1 * INTERPOLATION_POINTS + j1], 
                                               ct->xs[2 * INTERPOLATION_POINTS + k1], 
                                               cs->xs[0 * INTERPOLATION_POINTS + i2], 
                                               cs->xs[1 * INTERPOLATION_POINTS + j2], 
                                               cs->xs[2 * INTERPOLATION_POINTS + k2], 
                                               F);
                                ct->F[0 * NUM_SUB_MASSES + w1] += F[0] * cs->m[w2];
                                ct->F[1 * NUM_SUB_MASSES + w1] += F[1] * cs->m[w2];
                                ct->F[2 * NUM_SUB_MASSES + w1] += F[2] * cs->m[w2];
                            }
                        }
                    }
                }
            }
        }
    }

    void _eval_full(Cluster *ct, Cluster *cs){
        int i1, j1, i2, j2;
        bodies *bt, *bs;
        double F[3];
        double * m;

        bt = my_bs;
        bs = cs->active == world.rank?my_bs:NULL;
        m  = cs->active == world.rank?my_bs->m:recv_m;

        for(i1 = 0; i1 < ct->n; i1++){
            for(j1 = 0; j1 < cs->n; j1++){
                i2 = ct->start + i1;
                j2 = cs->start + j1;
                if(bs){
                    if(bt->id[i2] == bs->id[j2]){continue;}
                    calc_potential(bt->x[i2], bt->y[i2], bt->z[i2],
                                   bs->x[j2], bs->y[j2], bs->z[j2],
                                   F);
                } else {
                    calc_potential(bt->x[i2] , bt->y[i2] , bt->z[i2] ,
                                   recv_x[j2], recv_y[j2], recv_z[j2],
                                   F);
                }      
                bt->Fx[i2] += F[0] * m[j2];
                bt->Fy[i2] += F[1] * m[j2];
                bt->Fz[i2] += F[2] * m[j2]; 
            }
        }
    }

    void _eval(Cluster *ct, Cluster *cs){
        if(ct->active != semi_active && ct->active != world.rank){
            return;
        }

        if(admissable(ct, cs)){
            first_admissable = first_admissable <= current_level ? first_admissable : current_level;
            current_level++;
            _eval_CC(ct, cs);
            current_level--;
        } else {
            if (ct->num_sons && cs->num_sons){//both clusters have sons left:
            current_level++;
            _eval(ct->son[0], cs->son[0]);
            _eval(ct->son[0], cs->son[1]);
            _eval(ct->son[1], cs->son[0]);
            _eval(ct->son[1], cs->son[1]);
            current_level--;
            } 
            /*
                else if (ct->num_sons){         //just target cluster has sons left:
                par_eval(ct->son[0], cs);
                par_eval(ct->son[1], cs);
                } else if (cs->num_sons) {        //just source cluster has son left:
                par_eval(ct, cs->son[0]);
                par_eval(ct, cs->son[1]);
                }
            */ 
            else {                          //no son clusters left but still not admissable:
                _eval_full(ct, cs);
            }
        }
        // world.rank?:current_level?:printf("\nfirst admissable block on level %d\n", first_admissable);
        // world.rank?:current_level?:printf("split depth: %d\n", SPLIT_DEPTH);
    }

/************************************************************
 * Backward
 ***********************************************************/
    void _backward_Leaf(Cluster *c){
        int i, j, k, l, m, w, wj;
        double x, y, z, mass, l1, l2, l3;

        //for all linked masses
        for(l = 0; l < c->n; l++){
            m = c->start + l;
            x = my_bs->x[m];
            y = my_bs->y[m];
            z = my_bs->z[m];
            mass = my_bs->m[m];

            //for all substitution masses
            for(i = 0; i < INTERPOLATION_POINTS; i++){
                l1 = lagrange(c, i, 0, x);
                for(j = 0; j < INTERPOLATION_POINTS; j++){
                    wj = i * INTERPOLATION_POINTS + j;
                    l2 = l1 * lagrange(c, j, 1, y);
                    for(k = 0; k < INTERPOLATION_POINTS; k++){
                        w = wj * INTERPOLATION_POINTS + k;
                        l3 = l2 * lagrange(c, k, 2, z);

                        // calculate portion of this substitution masses force
                        // effecting this actual mass.
                        my_bs->Fx[m] += l3 * c->F[0 * NUM_SUB_MASSES + w];
                        my_bs->Fy[m] += l3 * c->F[1 * NUM_SUB_MASSES + w];
                        my_bs->Fz[m] += l3 * c->F[2 * NUM_SUB_MASSES + w];
                    }
                }
            }

            my_bs->Fx[m] *= mass * GAMMA;
            my_bs->Fy[m] *= mass * GAMMA;
            my_bs->Fz[m] *= mass * GAMMA; 
        }
    }

    void _backward_nonLeaf(Cluster *c){
        int i, j, k, ison, json, kson, w, wj, wson, wsj, son;
        Cluster *s;
        double l1, l2, l3;

        // for all fathers interpolation points
        for(i = 0; i < INTERPOLATION_POINTS; i++){
            for(j = 0; j < INTERPOLATION_POINTS; j++){
                wj = i * INTERPOLATION_POINTS + j;
                for(k = 0; k < INTERPOLATION_POINTS; k++){
                    w = wj * INTERPOLATION_POINTS + k;
                    //and all his sons
                    for(son = 0; son < 2; son++){
                        s = c->son[son];
                        //and all the sons interpolation points
                        for(ison = 0; ison < INTERPOLATION_POINTS; ison++){
                            l1 = lagrange(c, i, 0, s->xs[ison]);
                            for(json = 0; json < INTERPOLATION_POINTS; json++){
                                wsj = ison * INTERPOLATION_POINTS + json;
                                l2 = l1 * lagrange(c, j, 1, s->xs[INTERPOLATION_POINTS + json]);
                                for(kson = 0; kson < INTERPOLATION_POINTS; kson++){
                                    wson = wsj * INTERPOLATION_POINTS + kson;
                                    l3 = l2 * lagrange(c, k, 2, s->xs[2 * INTERPOLATION_POINTS + kson]);

                                    s->F[0 * NUM_SUB_MASSES + wson] += c->F[0 * NUM_SUB_MASSES + w] * l3;
                                    s->F[1 * NUM_SUB_MASSES + wson] += c->F[1 * NUM_SUB_MASSES + w] * l3;
                                    s->F[2 * NUM_SUB_MASSES + wson] += c->F[2 * NUM_SUB_MASSES + w] * l3; 
                                }
                            }
                        }
                    }
                }
            }
        }

        if(world.rank == c->son[0]->active || -1 == c->son[0]->active)
            backward(c->son[0]);
            
        if(world.rank == c->son[1]->active || -1 == c->son[1]->active)
            backward(c->son[1]);
    }

/************************************************************
 * Interface
 ***********************************************************/
    void forward(Cluster *c){
        work_time_start = MPI_Wtime();
        if(c->active != semi_active && c->active != world.rank){
            return;
        }
        if(c->num_sons){
            _forward_nonLeaf(c);
        } else {
            _forward_leaf(c);
        }
        work_time += MPI_Wtime() - work_time_start;
    }

    /** Method called from outside. Delegates the successive steps.
     * */
    void eval(Cluster *ct, Cluster *cs){
        prep_time_start = MPI_Wtime();
            _prep_comm(ct, cs);
            _prep_buffers();
        prep_time +=  MPI_Wtime() - prep_time_start;

        comm_time_start = MPI_Wtime();
            _communicate();
        comm_time += MPI_Wtime() - comm_time_start;

        prep_time_start = MPI_Wtime();
            _finalize_comm();
        prep_time += MPI_Wtime() - prep_time_start;

        work_time_start = MPI_Wtime();
            _eval(ct,cs);
        work_time += MPI_Wtime() - work_time_start;

        prep_time_start = MPI_Wtime();
            _clear_inactive_clusters();
        prep_time +=  MPI_Wtime() - prep_time_start;
    }

    void backward(Cluster *c){
        work_time_start = MPI_Wtime();
        if(c->active != semi_active && c->active != world.rank){
            return;
        }

        if(c->num_sons){
            _backward_nonLeaf(c);
        } else {
            _backward_Leaf(c);
        }
        work_time += MPI_Wtime() - work_time_start;
    }

/************************************************************
 * Init/Finalize
 ***********************************************************/
    void init_eval(){
        int i;
        send_sub_ms_count = calloc(world.size, sizeof(int));
        send_sub_ms_displ = calloc(world.size, sizeof(int));
        recv_sub_ms_count = calloc(world.size, sizeof(int));
        recv_sub_ms_displ = calloc(world.size, sizeof(int));

        send_bs_count = calloc(world.size, sizeof(int));
        send_bs_displ = calloc(world.size, sizeof(int));
        recv_bs_count = calloc(world.size, sizeof(int));
        recv_bs_displ = calloc(world.size, sizeof(int));

        send_n_count = calloc(world.size, sizeof(int));
        send_n_displ = calloc(world.size, sizeof(int));
        recv_n_count = calloc(world.size, sizeof(int));
        recv_n_displ = calloc(world.size, sizeof(int));

        send_clusters = malloc(world.size * sizeof(vector*));
        recv_clusters = malloc(world.size * sizeof(vector*));
        send_bodies   = malloc(world.size * sizeof(vector*));
        recv_bodies   = malloc(world.size * sizeof(vector*));
        for(i = 0; i < world.size; i++){
            send_clusters[i] = malloc(sizeof(vector));
            recv_clusters[i] = malloc(sizeof(vector));
            send_bodies[i]   = malloc(sizeof(vector));
            recv_bodies[i]   = malloc(sizeof(vector));
            vector_init(send_clusters[i]);
            vector_init(recv_clusters[i]);
            vector_init(send_bodies[i]);
            vector_init(recv_bodies[i]);
        }

        send_sub_ms = NULL;
        recv_sub_ms = NULL;
        send_x      = NULL;
        send_y      = NULL;
        send_z      = NULL;
        send_m      = NULL;
        send_n      = NULL;
        recv_x      = NULL;
        recv_y      = NULL;
        recv_z      = NULL;
        recv_m      = NULL;
        recv_n      = NULL;

        //init stopwatches
        work_time  = 0.0;
        comm_time  = 0.0;
        prep_time  = 0.0;
    }

    void finalize_eval(){
        int i;
        for(i = 0; i < world.size; i++){
            vector_free(send_clusters[i]);
            vector_free(recv_clusters[i]);
            vector_free(send_bodies[i]);
            vector_free(recv_bodies[i]);
            free(send_clusters[i]);
            free(recv_clusters[i]);
            free(send_bodies[i]  );
            free(recv_bodies[i]  );
        }
        free(send_clusters);
        free(recv_clusters);
        free(send_bodies  );
        free(recv_bodies  );

        free(send_sub_ms_count);
        free(send_sub_ms_displ);
        free(recv_sub_ms_count);
        free(recv_sub_ms_displ);

        free(send_bs_count);
        free(send_bs_displ);
        free(recv_bs_count);
        free(recv_bs_displ);

        free(send_n_count);
        free(send_n_displ);
        free(recv_n_count);
        free(recv_n_displ);
        
        free(send_sub_ms);
        free(recv_sub_ms);
        free(send_x);
        free(send_y);
        free(send_z);
        free(send_m);
        free(recv_x);
        free(recv_y);
        free(recv_z);
        free(recv_m);
    }