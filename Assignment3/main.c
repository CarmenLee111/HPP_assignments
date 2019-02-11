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

#define e 0.001 


typedef struct celestial_bodies{
    double                  x;
    double                  y;
    double                  mass;
    double                  vx;
    double                  vy;
    double                  brtness;
}body;



/* from compare_gal_files, should add to header files */
void load_data(int n, body *B, char* fileName);
void body_info(body *b);
void step(double G, int N, double dt, body *B);


int main(int argc, char *argv[]) {

    if (argc != 6) {
        printf("Usage: ./galsim N filname nsteps delta_t graphics\n");
        return -1;
    }

    int N = atoi(argv[1]);
    double const G = (double) 100/N;
    //printf("Number of bodies: %d\n", N);
    char *filename = argv[2];
    body *B = malloc(N * sizeof(body));
    load_data(N, B, filename);
    int n = atoi(argv[3]);        // number of time steps
    double dt = atof(argv[4]);    // time step
    //int graphics = atoi(argv[5]); // graphics on or off  

    /* time steps of the simulaiton */
    /* graphics */
    const float circleRadius=0.025, circleColor=0;
    const int windowWidth=800;
    float L = 1, W = 1;
    InitializeGraphics(argv[0], windowWidth, windowWidth);
    SetCAxes(0,1);

    //printf("Hit q to quit.\n");
    // while(!CheckForQuit())

    int i,j;
    for (j=0; j<n; j++) {

        for (i=0; i<N; i++) {
            DrawCircle(B[i].x, B[i].y, L, W, circleRadius, circleColor);
        }
        ClearScreen();
        Refresh();
        usleep(3000);

        step(G, N, dt, B);
        // print results for each time step 
        /*
        printf("After time step %d\n", j+1);
        for (i=0; i<N; i++) {
            printf("-------------- Body %d --------------\n", i);
            body_info(&B[i]);
        }
        */
    }

    FlushDisplay();
    CloseDisplay();


    FILE *fp;
    fp = fopen("result.gal", "w");
    fwrite(B, N*sizeof(body), 1, fp);
    fclose(fp);

    // Check the output file
    body *B2 = malloc(N * sizeof(body));
    printf("******* result file: *******\n");
    load_data(N, B2, "result.gal");
    free(B2);

#if 1     // inspect data in ellipse_N_00010_after200steps
    body *B3 = malloc(N * sizeof(body));
    printf("******* compare file ellipse_N_00010_after200steps: ******\n");
    load_data(N, B3, "ref_output_data/ellipse_N_00010_after200steps.gal");
    free(B3);
#endif  


    
    free(B);

    return 0;
}


/* Stepping function that computes/assigns the new velocity and position
 */
void step(double G, int N, double dt, body *B){
    int i,j;
    
    for (i=0; i<N; i++) {
        double rx=0, ry=0, r=0, Fx=0, Fy=0, fabs=0, ax=0, ay=0;
        // printf("------Computation for Body %d------\n", i);
        for (j=0; j<N; j++) {
            if (i == j) {;}
            else {
                rx = B[i].x - B[j].x;
                ry = B[i].y - B[j].y;
                r = sqrt(rx*rx + ry*ry);
                // printf("Distance vector (%f, %f)\n", rx, ry);
                fabs = - G * B[i].mass * B[j].mass / ((r+e)*(r+e)*(r+e));
                Fx += fabs * rx; 
                Fy += fabs * ry;
                // printf("Force vector (%f, %f)\n", Fx, Fy);
            }
        }
        ax = Fx / B[i].mass;
        ay = Fy / B[i].mass;
        // printf("Acceleration vector (%f, %f)\n", ax, ay);
        B[i].vx += dt * ax;
        B[i].vy += dt * ay;
        B[i].x += dt * B[i].vx;
        B[i].y += dt * B[i].vy;
        // printf("New velocity vector (%f, %f)\n", B[i].vx, B[i].vy);
        // printf("New position vector (%f, %f)\n", B[i].x, B[i].y);

    }
}



void body_info(body *b) {
    printf("Position: \t(%f, %f)\nMass: \t\t%f\nVelocity: \t(%f, %f)\nBrtness: \t%f\n", 
        b->x, b->y, b->mass, b->vx, b->vy, b->brtness);
}


void load_data(int N, body *B, char* fileName) {
    FILE *fp = fopen(fileName, "rb");
    if (!fp) {
        printf("load_data error: failed to open input file '%s'.\n", fileName);
        return;
    }
    
    fread(B, N*sizeof(body), 1, fp);
    
    // Printing the info 
    
    int i;
    for (i=0; i<N; i++) {
        body_info(&(B[i]));
    }

    fclose(fp);
}



