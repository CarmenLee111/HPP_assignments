#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "graphics/graphics.h"


/* Galaxy simulation using Barnes-Hut Algorithm
 * Input
 * N            number of stars/particles to simulate
 * filname      file to read the initial configuration from
 * nsteps       number of timesteps
 * delta_t      time step, set to 1e-5
 * theta_max    Controls the level of approximation of Barnes-Hut 
 * graphics     1 or 0: graphics on/off
 * n_threads    Number of threads
 */

#define E 0.001 
#define debug 0

typedef struct vectors
{
    double                  x;
    double                  y;    
}vector;

typedef struct quadtree
{
    struct  quadtree*   nw;
    struct  quadtree*   ne;
    struct  quadtree*   sw;
    struct  quadtree*   se;
    vector              c;    // center of mass
    vector              p;    // position, center of the bounding box
    double              m;
    float               h;    // half box width

} qt;


void load_data(int N, vector* c, vector* v, double* m, double *b, char* filename);
void write_data(int N, vector* c, vector* v, double* m, double *b, char* filename);
qt *new_qt(vector c, double m);
void qt_add(qt **t, qt *tadd);
void qt_del(qt **t);
void qt_print(qt *t);
void force_BH_one_body(double G, vector* c, double *m, qt *t, vector *f, double theta);
int check_quadrant(vector c, vector p);
int timer(void);


int main(int argc, char *argv[]) {

    if (argc != 8) {
        printf("Usage: ./galsim N filname nsteps delta_t theta_max graphics n_threads\n");
        return -1;
    }

    const int N = atoi(argv[1]);
    double const G = (double) 100/N;
    char *filename = argv[2];
    const int n = atoi(argv[3]);                  // number of time steps
    const double dt = atof(argv[4]);              // time step
    const float theta = atof(argv[5]);
    int graphics = atoi(argv[6]);           // graphics not available    
    const int n_threads = atoi(argv[7]);
        int starttime, ttime=0;                       // Wall time log

    if (graphics) {printf("Graphics not available in this version\n");}


    vector *C = malloc(N * sizeof(vector));
    vector *V = malloc(N * sizeof(vector));
    double *M = malloc(N * sizeof(double));
    double *b = malloc(N * sizeof(double));    
    load_data(N, C, V, M, b, filename);


    // Stepping
    int i,j;
    for (j=0; j<n; j++) {
        qt *t = NULL;                           // reinitiate the tree root
        vector* F = calloc(N,sizeof(vector));   // reinitiate the force vector
        
        // build a new tree
        for (i=0; i<N; i++) {
            qt_add(&t, new_qt(C[i], M[i]));
        }
    #if debug
        qt_print(t);
    #endif

    starttime = timer();
#   pragma omp parallel num_threads(n_threads)
    {
#       pragma omp for 
        for (i=0; i<N; i++) {
            force_BH_one_body(G, &(C[i]), &(M[i]), t, &(F[i]), theta);
        }

#       pragma omp for 
        for (i=0; i<N; i++) {
            V[i].x +=  dt * F[i].x;
            V[i].y +=  dt * F[i].y;
            C[i].x +=  V[i].x * dt;
            C[i].y +=  V[i].y * dt;
        }
    }
    ttime += timer()-starttime;

        free(F);
        qt_del(&t);
        free(t);           // Added since submission
    }

    printf("Wall time for the parallelization is: %fs\n", ttime/1000000.0);

    char outputfile[11] = "result.gal";
    write_data(N, C, V, M, b, outputfile);

    free(C); free(V), free(M), free(b);
    return 0;
}



/* Effectively computes the acceleration for one body */
void force_BH_one_body(double G, vector *c, double *m, qt *t, vector *f, double theta) {
    // node does not exist
    if (!t) {
        ;
    }

    // external branch
    else if ((t->nw==NULL) && (t->ne==NULL) && (t->sw==NULL) && (t->se==NULL)) {
        double rx = c->x - t->c.x;
        double ry = c->y - t->c.y;
        double r = sqrt(rx*rx + ry*ry);
        double comm = -G * t->m / ((r+E)*(r+E)*(r+E));

        f->x += comm * rx;
        f->y += comm * ry;
    }

    // internal branch
    else {                                 // Calcuate distance to the mid-point of the box
        double dx = c->x - t->p.x;
        double dy = c->y - t->p.y;
        double d = sqrt(dx*dx + dy*dy);

        if ((t->h)*2 < theta * d) {        // theta criterion satisfies, treat as one body

            double rx = c->x - t->c.x;
            double ry = c->y - t->c.y;
            double r = sqrt(rx*rx + ry*ry);
            // double comm = - G * b->m * t->m / ((r+E)*(r+E)*(r+E));
            double comm = - G * t->m / ((r+E)*(r+E)*(r+E));

            f->x += comm * rx;
            f->y += comm * ry;
        } 
        else {                             // theta criterion fails, traverse the subnodes
            force_BH_one_body(G, c, m, t->sw, f, theta);
            force_BH_one_body(G, c, m, t->nw, f, theta);
            force_BH_one_body(G, c, m, t->se, f, theta);
            force_BH_one_body(G, c, m, t->ne, f, theta);
        } 
    }
}

// Constructor for a new quadtree node
qt *new_qt(vector c, double m) {
    qt *newqt = malloc(sizeof(qt));
    newqt->c  = c;
    newqt->m  = m;
    newqt->nw = NULL;
    newqt->ne = NULL;
    newqt->sw = NULL;
    newqt->se = NULL;
    newqt->p.x=0.5;   // Initialize to always start from the root quadrant
    newqt->p.y=0.5;
    newqt->h=0.5;
    return newqt;
}



void qt_add(qt **t, qt *tadd){
     // Quadrant free
    if(!(*t)) {                    
        *t = tadd;
        return;
    } 

    // Current node is internal
    else if((((*t)->nw) || ((*t)->ne) || ((*t)->sw) || ((*t)->se))) {  

        double m1 = (*t)->m;
        
        // update current node mass and CoM
        (*t)->m += tadd->m;
        (*t)->c.x = (((*t)->c.x)*m1 + (tadd->c.x)*(tadd->m))/((*t)->m);
        (*t)->c.y = (((*t)->c.y)*m1 + (tadd->c.y)*(tadd->m))/((*t)->m);

        // update new child half box width
        tadd->h = ((*t)->h)/2;

        // Insert the new node into the correct quadrant
        switch (check_quadrant(tadd->c, (*t)->p)) {
            case 1: {
                tadd->p.x = ((*t)->p.x) - tadd->h;
                tadd->p.y = ((*t)->p.y) + tadd->h;
                qt_add(&((*t)->nw), tadd);
                break;
            }
            case 2: {
                tadd->p.x = ((*t)->p.x) + tadd->h;
                tadd->p.y = ((*t)->p.y) + tadd->h;
                qt_add(&((*t)->ne), tadd);
                break;
            }
            case 3: {
                tadd->p.x = ((*t)->p.x) - tadd->h;
                tadd->p.y = ((*t)->p.y) - tadd->h;
                qt_add(&((*t)->sw), tadd);break;
            }
            case 4: {
                tadd->p.x = ((*t)->p.x) + tadd->h;
                tadd->p.y = ((*t)->p.y) - tadd->h;
                qt_add(&((*t)->se), tadd);break;
            }
            default:break;
        }
        return;
    }

    // Current node is an external node
    else { 

        // copying the current node with half of the box width
        // ready to move down
        qt *tnew =  new_qt((*t)->c, (*t)->m);
        tnew->h = ((*t)->h)/2;

        // update new child half box width
        tadd->h = ((*t)->h)/2;

        // update current node with new mass and CoM
        double m1 = (*t)->m;
        (*t)->m += tadd->m;
        (*t)->c.x = (((*t)->c.x)*m1 + (tadd->c.x)*(tadd->m))/((*t)->m);
        (*t)->c.y = (((*t)->c.y)*m1 + (tadd->c.y)*(tadd->m))/((*t)->m);

        // Check the quadrant and insert the original node there
        switch (check_quadrant(tnew->c, (*t)->p)) {
            case 1: {
                tnew->p.x = ((*t)->p.x) - tnew->h;
                tnew->p.y = ((*t)->p.y) + tnew->h;
                (*t)->nw = tnew;
                break;
            }
            case 2: {
                tnew->p.x = ((*t)->p.x) + tnew->h;
                tnew->p.y = ((*t)->p.y) + tnew->h;
                (*t)->ne = tnew;
                break;
            }
            case 3: {
                tnew->p.x = ((*t)->p.x) - tnew->h;
                tnew->p.y = ((*t)->p.y) - tnew->h;
                ((*t)->sw) = tnew;
                break;
            }
            case 4: {
                tnew->p.x = ((*t)->p.x) + tnew->h;
                tnew->p.y = ((*t)->p.y) - tnew->h;
                (*t)->se = tnew;
                break;
            }
            default:break;
        }

        // Check the quadrant and insert the new node there
        switch (check_quadrant(tadd->c, (*t)->p)) {
            case 1: {
                tadd->p.x = ((*t)->p.x) - tadd->h;
                tadd->p.y = ((*t)->p.y) + tadd->h;
                qt_add(&((*t)->nw), tadd);
                break;
            }
            case 2: {
                tadd->p.x = ((*t)->p.x) + tadd->h;
                tadd->p.y = ((*t)->p.y) + tadd->h;
                qt_add(&((*t)->ne), tadd);
                break;
            }
            case 3: {
                tadd->p.x = ((*t)->p.x) - tadd->h;
                tadd->p.y = ((*t)->p.y) - tadd->h;
                qt_add(&((*t)->sw), tadd);break;
            }
            case 4: {
                tadd->p.x = ((*t)->p.x) + tadd->h;
                tadd->p.y = ((*t)->p.y) - tadd->h;
                qt_add(&((*t)->se), tadd);break;
            }
            default:break;
        }
        return;
    }
}

int check_quadrant(vector c, vector p) {
    //first quadrant NW
    if (c.x <= p.x && c.y >  p.y) return 1;  // NW
    if (c.x >  p.x && c.y >= p.y) return 2;  // NE
    if (c.x <  p.x && c.y <= p.y) return 3;  // SW
    if (c.x >= p.x && c.y <  p.y) return 4;  // SE
    else return 0;
}

void qt_del(qt **t) {
    if ((*t) != NULL) {

        if ((*t)->nw) {
            qt_del(&((*t)->nw));
            free((*t)->nw);
        }
        if ((*t)->ne) {
            qt_del(&((*t)->ne));
            free((*t)->ne);
        }
        if ((*t)->sw) {
            qt_del(&((*t)->sw));
            free((*t)->sw);
        }
        if ((*t)->se) {
            qt_del(&((*t)->se));
            free((*t)->se);
        }

        free(*t);
        *t = NULL;
    }   
}


/* Print out the quadtree, for debugging purpose*/
void qt_print(qt *t) {
    if (t == NULL) {;}
    else {
        printf("Node mass: %f\n", t->m);
        printf("Node position: (%f, %f)\n", t->p.x, t->p.y);
        printf("Node center of mass: (%f, %f)\n", t->c.x, t->c.y);
        printf("Bounding box height: %f\n", t->h);
        qt_print(t->nw);
        qt_print(t->ne);
        qt_print(t->sw);
        qt_print(t->se);
    }
}


void load_data(int N, vector* C, vector* V, double* M, double* b, char* filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        printf("load_data error: failed to open input file '%s'.\n", filename);
        return;
    }
    int i;
    for (i=0; i<N; i++) {
        fread(&(C[i]), sizeof(double), 2, fp);
        fread(&(M[i]), sizeof(double), 1, fp);
        fread(&(V[i]), sizeof(double), 2, fp);
        fread(&(b[i]), sizeof(double), 1, fp);
    }
    fclose(fp);
}


void write_data(int N, vector* C, vector* V, double* M, double* b, char* filename) {
    FILE *fp = fopen(filename, "w");
    int i;
    for (i=0; i<N; i++) {
        fwrite(&(C[i]), sizeof(double), 2, fp);
        fwrite(&(M[i]), sizeof(double), 1, fp);
        fwrite(&(V[i]), sizeof(double), 2, fp);
        fwrite(&(b[i]), sizeof(double), 1, fp);
    }
    fclose(fp);
}


int timer(void)
{
  struct timeval tv;
  gettimeofday(&tv, (struct timezone*)0);
  return (tv.tv_sec*1000000+tv.tv_usec);
}



