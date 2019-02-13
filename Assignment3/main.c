#include <math.h>
#include "graphics/graphics.h"


/* Galaxy simulation
 * Input
 * N        number of stars/particles to simulate
 * filname  file to read the initial configuration from
 * nsteps   number of timesteps
 * delta_t  time step, set to 1e-5
 * graphics 1 or 0: graphics on/off
 */

#define E 0.001 
#define STRUCTURE_OPTION 0      // 1 in order, 0 out of order.

#if STRUCTURE_OPTION            // selecting the structure option
typedef struct celestial_bodies
{
    double                  x;
    double                  y;
    double                  mass;
    double                  vx;
    double                  vy;
    double                  brtness;
}body;
#else 
typedef struct celestial_bodies
{
    double                  x;
    double                  y;
    double                  vx;
    double                  vy;
    double                  brtness;
    double                  mass;
    double                  dummy1;
    double                  dummy2;
}body;
#endif


typedef struct vectors
{
    double                  x1;
    double                  x2;
}vector;

const float circleRadius=0.005, circleColor=0;
const int windowWidth=800;

/* from compare_gal_files, should add to header files */
void load_data(int n, body *B, char* fileName);
void body_info(body *b);
void step(double G, int N, double dt, body *B);
void load_data_outOfOrder(int n, body *B, char* fileName);


int main(int argc, char *argv[]) {

    if (argc != 6) {
        printf("Usage: ./galsim N filname nsteps delta_t graphics\n");
        return -1;
    }

    int N = atoi(argv[1]);
    double const G = (double) 100/N;
    char *filename = argv[2];
    body *B = malloc(N * sizeof(body));
    int n = atoi(argv[3]);                  // number of time steps
    double dt = atof(argv[4]);              // time step
    int graphics = atoi(argv[5]);           // graphics on or off  
    
    if (!graphics){

#if STRUCTURE_OPTION    // option to load into the memory as it is ordered
        load_data(N, B, filename);

        int i,j;
        for (j=0; j<n; j++) {

            step(G, N, dt, B);
        }

        /* Write the result to a output binray file */
        FILE *fp;
        fp = fopen("result.gal", "w");
        fwrite(B, N*sizeof(body), 1, fp);
        fclose(fp);
#else                   // option to load into the memory in a different order 
        load_data_outOfOrder(N, B, filename);

        int i,j;
        for (j=0; j<n; j++) {

            step(G, N, dt, B);
        }

        FILE *fp;
        fp = fopen("result.gal", "w");

        for (i=0; i<N; i++) {
            fwrite(&(B[i].x), sizeof(double), 1, fp);
            fwrite(&(B[i].y), sizeof(double), 1, fp);
            fwrite(&(B[i].mass), sizeof(double), 1, fp);
            fwrite(&(B[i].vx), sizeof(double), 1, fp);
            fwrite(&(B[i].vy), sizeof(double), 1, fp);
            fwrite(&(B[i].brtness), sizeof(double), 1, fp);
        }
        fclose(fp);
#endif

    } else {
        load_data(N, B, filename);

        float L = 1, W = 1;
        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
        printf("Ctrl C to quit.\n");

        /* time steps of the simulaiton with optional graphics stuff */
        int i,j;
        for (j=0; j<n; j++) {
            ClearScreen();
            for (i=0; i<N; i++) {
                DrawCircle(B[i].x, B[i].y, L, W, circleRadius, circleColor);
            }
            Refresh();
            usleep(3000);
            step(G, N, dt, B);
        }
        
        while(!CheckForQuit()){
            usleep(200000);
        }
        FlushDisplay();
        CloseDisplay();
    }   
 
    free(B);
    return 0;
}


/* Stepping function that computes/assigns the new velocity and position
 */
void step(double G, int N, double dt, body *B){
    int i,j;
    double rx=0, ry=0, r=0, fabs=0, ax=0, ay=0;
    vector *F = calloc(N, sizeof(vector));

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (i == j) {;}
            else {
                rx = B[i].x - B[j].x;
                ry = B[i].y - B[j].y;
                r = sqrt(rx*rx + ry*ry);
                fabs = - G * B[i].mass * B[j].mass / ((r+E)*(r+E)*(r+E));

                F[i].x1 += fabs * rx; 
                F[i].x2 += fabs * ry;
            }
        }
    }

    /* For-loop to update the value (Potential room for optimization????? */ 
    for (i=0; i<N; i++) {
        ax = F[i].x1 / B[i].mass;
        ay = F[i].x2 / B[i].mass;
        B[i].vx += dt * ax;
        B[i].vy += dt * ay;
        B[i].x += dt * B[i].vx;
        B[i].y += dt * B[i].vy;
    }
    free(F);
}


/* Printing information of bodies */
void body_info(body *b) {
    printf("Position: \t(%f, %f)\nMass: \t\t%f\nVelocity: \t(%f, %f)\nBrtness: \t%f\n", 
        b->x, b->y, b->mass, b->vx, b->vy, b->brtness);
}

/* Load binary files */
void load_data(int N, body *B, char* fileName) {
    FILE *fp = fopen(fileName, "rb");
    if (!fp) {
        printf("load_data error: failed to open input file '%s'.\n", fileName);
        return;
    }
    fread(B, N*sizeof(body), 1, fp);
    fclose(fp);
}


/* Load binary files out of order */
void load_data_outOfOrder(int N, body *B, char* fileName) {
    FILE *fp = fopen(fileName, "rb");
    if (!fp) {
        printf("load_data error: failed to open input file '%s'.\n", fileName);
        return;
    }
    int i;
    for (i=0; i<N; i++) {
        fread(&(B[i].x), sizeof(double), 1, fp);
        fread(&(B[i].y), sizeof(double), 1, fp);
        fread(&(B[i].mass), sizeof(double), 1, fp);
        fread(&(B[i].vx), sizeof(double), 1, fp);
        fread(&(B[i].vy), sizeof(double), 1, fp);
        fread(&(B[i].brtness), sizeof(double), 1, fp);
    }
    fclose(fp);
}



