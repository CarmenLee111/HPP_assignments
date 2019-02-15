#include <math.h>
#include "graphics/graphics.h"
#include <sys/time.h>


/* Galaxy simulation
 * Input
 * N        number of stars/particles to simulate
 * filname  file to read the initial configuration from
 * nsteps   number of timesteps
 * delta_t  time step, set to 1e-5
 * graphics 1 or 0: graphics on/off
 */

#define E 0.001 
#define STRUCTURE_OPTION 1     // 1 with structure, 0 with arrays

typedef struct vectors
{
    double                  v[2];
}vector;

typedef struct celestial_bodies
{
    double                  x;
    double                  y;
    double                  mass;
    double                  vx;
    double                  vy;
    double                  brtness;
}body;

#if STRUCTURE_OPTION            // selecting the structure option

void load_data(int n, body *B, char* fileName);
void step(double G, int N, double dt, body *B);



#else 

void load_data_noStruct(int N, vector *x, double *m, vector *v, double *b, char* fileName);
void step2(double G, int N, double dt, vector *x, double *m, vector *v);


#endif

static double get_wall_seconds();



int main(int argc, char *argv[]) {

    if (argc != 6) {
        printf("Usage: ./galsim N filname nsteps delta_t graphics\n");
        return -1;
    }

    int N = atoi(argv[1]);
    double const G = (double) 100/N;
    char *filename = argv[2];
    int n = atoi(argv[3]);                  // number of time steps
    double dt = atof(argv[4]);              // time step
    int graphics = atoi(argv[5]);           // graphics on or off  
    body *B = malloc(N * sizeof(body));
    vector *x = malloc(N*sizeof(vector));
    vector *v = malloc(N*sizeof(vector));
    double *m = malloc(N*sizeof(double));
    double *b = malloc(N*sizeof(double));

    
    if (!graphics){

#if STRUCTURE_OPTION    // option to load into the memory as it is ordered

        double startTime1 = get_wall_seconds();
        load_data(N, B, filename);
        double secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for loading all into B: %f\n", secondsTaken1);

        startTime1 = get_wall_seconds();
        int j;
        for (j=0; j<n; j++) {
            step(G, N, dt, B);
        }
        secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for stepping: %f\n", secondsTaken1);

        startTime1 = get_wall_seconds();
        /* Write the result to a output binray file */
        FILE *fp;
        fp = fopen("result.gal", "w");
        fwrite(B, N*sizeof(body), 1, fp);
        fclose(fp);
        secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for writing: %f\n", secondsTaken1);


#else                   // option to gather the same variable in one array. Much slower for some reason 

        double startTime1 = get_wall_seconds();
        load_data_noStruct(N, x, m, v, b, filename);
        double secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for loading all into B into arrays: %f\n", secondsTaken1);

        startTime1 = get_wall_seconds();
        int i,j;
        for (j=0; j<n; j++) {
            step2(G, N, dt, x, m, v);
        }
        secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for stepping with arrays %f\n", secondsTaken1);

        startTime1 = get_wall_seconds();
        FILE *fp;
        fp = fopen("result.gal", "w");

        for (i=0; i<N; i++) {
            fwrite(&(x[i]), sizeof(double), 2, fp);
            fwrite(&(m[i]), sizeof(double), 1, fp);
            fwrite(&(v[i]), sizeof(double), 2, fp);
            fwrite(&(b[i]), sizeof(double), 1, fp);
        }
        fclose(fp);
        secondsTaken1 = get_wall_seconds() - startTime1;
        printf("secondsTaken for writing %f\n", secondsTaken1);




#endif

    } else {
        
#if STRUCTURE_OPTION
        const float circleRadius=0.003, circleColor=0.2;
        const int windowWidth=800;
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
                DrawCircle(B[i].x, B[i].y, L, W, circleRadius*fabs(B[i].mass), circleColor);
            }
            Refresh();
            usleep(4000);
            step(G, N, dt, B);

        }
        
        while(!CheckForQuit()){
            usleep(200000);
        }
        FlushDisplay();
        CloseDisplay();

#else 
        printf("Change STRUCTURE_OPTION to enable graphics\n");

#endif



    }   
 


    free(x);
    free(v);
    free(m);
    free(b);
    free(B);
    return 0;
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

/* Stepping function that computes/assigns the new velocity and position
 */
void step(double G, int N, double dt, body *B){
    int i,j;
    double rx, ry, r, fab;
    vector *F = calloc(N, sizeof(vector));

    double Gdt = G * dt;

    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            rx = B[i].x - B[j].x;
            ry = B[i].y - B[j].y;
            r = sqrt(rx*rx + ry*ry);
            fab = - B[i].mass * B[j].mass / ((r+E)*(r+E)*(r+E));

            (F[i]).v[0] += fab * rx; 
            (F[i]).v[1] += fab * ry;
            (F[j]).v[0] -= fab * rx;
            (F[j]).v[1] -= fab * ry;
        }

        double invm = 1 / B[i].mass;
        B[i].vx += Gdt * (F[i]).v[0] * invm;
        B[i].vy += Gdt * (F[i]).v[1] * invm;
        B[i].x += dt * B[i].vx;
        B[i].y += dt * B[i].vy;

    }

    // for (i=0; i<N; i++) {
    //     double invm = 1 / B[i].mass;
    //     B[i].vx += Gdt * (F[i]).v[0] * invm;
    //     B[i].vy += Gdt * (F[i]).v[1] * invm;
    //     B[i].x += dt * B[i].vx;
    //     B[i].y += dt * B[i].vy;
    // }

    free(F);
}


/* load_data function for the array option */
void load_data_noStruct(int N, vector *x, double *m, vector *v, double *b, char* fileName) {
    FILE *fp = fopen(fileName, "rb");
    if (!fp) {
        printf("load_data error: failed to open input file '%s'.\n", fileName);
        return;
    }
    int i;
    for (i=0; i<N; i++) {
        fread(&(x[i]), sizeof(double), 2, fp);
        fread(&(m[i]), sizeof(double), 1, fp);
        fread(&(v[i]), sizeof(double), 2, fp);
        fread(&(b[i]), sizeof(double), 1, fp);
    }
    fclose(fp);
}

/* Stepping function for the array option */
void step2(double G, int N, double dt, vector *x, double *m, vector *v){
    int i,j;
    double rx=0, ry=0, r=0;
    vector *F = calloc(N, sizeof(vector));
    const double Gdt = G * dt;

    for (i=0; i<N; i++) {
      
        for (j=i+1; j<N; j++) {
            rx = x[i].v[0] - x[j].v[0];
            ry = x[i].v[1] - x[j].v[1];
            r = sqrt(rx*rx + ry*ry);
            double fab = - m[i] * m[j] /((r+E)*(r+E)*(r+E));

            F[i].v[0] += fab * rx; 
            F[i].v[1] += fab * ry;
            F[j].v[0] -= fab * rx;
            F[j].v[1] -= fab * ry;

        }
        double invm = 1 / m[i];

        v[i].v[0] += Gdt * F[i].v[0] * invm;
        v[i].v[1] += Gdt * F[i].v[1] * invm;

        x[i].v[0] += dt * v[i].v[0];
        x[i].v[1] += dt * v[i].v[1];
    }

    // for (i=0; i<N; i++) {
    //     double invm = 1 / m[i];

    //     v[i].v[0] += Gdt * F[i].v[0] * invm;
    //     v[i].v[1] += Gdt * F[i].v[1] * invm;

    //     x[i].v[0] += dt * v[i].v[0];
    //     x[i].v[1] += dt * v[i].v[1];
    // }

    free(F);
}


static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}