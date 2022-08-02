//
//  tokamak.c
//
//  A program which solves a non-linear diffusion equation describing the
//  temperature of a plasma in a tokamak nuclear fusion reactor with given
//  initial/boundary conditions.
//  Von Neumann/Fourier analysis was used to find a suitable time-stepping scheme.
//
//  Created by Chance Haycock on 2018-11-16.
//  Copyright 2018 Chance Haycock. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/**
 * Returns the maximum (in absolute value) element in
 * an array. Used for time step stability.
 */
double max_array(double *arr, long nelements) {
    long i;
    double maximum = fabs(arr[0]);
    for (i = 0; i < nelements; i++) {
        if (fabs(arr[i]) > maximum) {
            maximum = fabs(arr[i]);
        }
    }
    return maximum;
}

// x-Domain Boundary
double L;
// Number of grid points
int N;
// Length of time to run
double t_F;
// Time interval to output
double t_d;
// Initial condition parameter
double T_0;
// Boundary condition parameter
double A;
// Boundary condition parameter
double w;
// Non-linearity exponent
double P;

/**
 * Read data from text file.
 */
void read_input(double *L, int *N, double *t_F, double *t_d, double *T_0,
                double *A, double *w, double *P) {
    FILE *infile;
    if(!(infile=fopen("/Users/ajanaranyaksiya/Desktop/testing/testingpx3902018/testingpx3902018/input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    if(8!=fscanf(infile,"%lf %d %lf %lf %lf %lf %lf %lf", L, N, t_F, t_d, T_0, A, w, P)) {
        printf("Error reading parameters from file\n");
        exit(1);
    }
    fclose(infile);
}

int main() {
    // Read data in from file
    read_input(&L, &N, &t_F, &t_d, &T_0, &A, &w, &P);

    // Input Check
    if (L <=0) {
        printf("x Domain must have positive length.\n");
        exit(1);
    }
    if (N <= 1) {
        printf("Program requires more than 1 grid point.\n");
        exit(1);
    }

    // Space/Time Steps
    double dx = L / (N-1);
    double dt = 1e-5; // starting point - later depends on temp gradient

    double *T; // Current temperature
    double *T_next; // Temperature at next time step
    double *dTdx; // Temperature Gradient
    double *gammaa; // As specification
    double *dgamdx; // Spatial gradient of gamma.

    T = malloc(sizeof(double) * N);
    T_next = malloc(sizeof(double) * N);
    dTdx = malloc(sizeof(double) * N);
    gammaa = malloc(sizeof(double) * N);
    dgamdx = malloc(sizeof(double) * N);

    // Initial conditions
    double x = 0.0;
    double t = 0.0;
    double next_output_time = t_d;

    FILE *fp; // opening file to check plotting
    fp = fopen ("/Users/ajanaranyaksiya/Desktop/testing/testingpx3902018/testingpx3902018/output.txt", "w+");

    // Fill array with t=0 linear function from specification
    // T(x, 0) = T_0 * (1 - x/L)
    for (int j=0; j < N; j++) {
        x = j*dx;
        T[j] = T_0 * (1 - (x/L));  // Initial function at t=0
        fprintf(fp,"%g\t%g\t%g\n",t,x,T[j]); // <------ First line in file
    }

    // Main time loop
    while(t < t_F) {
        double dt0 = dt;
        if (t+dt0 > next_output_time) {
            dt0 = next_output_time - t;
        }
 
        // Filling derivative array
        dTdx[0] = (T[1] - T[0])/(dx); // non centered at the boundary

        // Calculating deriv at all other non-boundary points
        for (int j=1; j<N-1; j++) {
            x = j*dx;
            int jm = j-1; // space
            int jp = j+1;
            dTdx[j] = (T[jp] - T[jm])/(2.0*dx); // centred second order first derivative
        }
        dTdx[N-1] = -T_0 / L; // Important b.c -- dT/dx(L, t) = -T_0/L for all t

        // Filling Gamma array including sign function
        for(int i = 0; i < N; i++) {
            int sign = 0;
            if (dTdx[i] > 0) {
                sign = 1;
            }
            else {
                sign = -1;
            }
            gammaa[i] = -sign * pow(fabs(dTdx[i]), P);
        }

        T_next[0] = A * (1 - cos(w*t)) + T_0; // boundary condition at x = 0 for all t

        dgamdx[0] = (gammaa[1] - gammaa[0])/(dx); // non centered derivative at boundaries
        dgamdx[N-1] = (gammaa[N-1] - gammaa[N-2])/dx;

        for (int j=1 ; j<N-1; j++) {
            x = j*dx;
            int jm = j-1; // space
            int jp = j+1;
            dgamdx[j] = (gammaa[jp] - gammaa[jm])/(2*dx);
            T_next[j] = T[j]-dt0*dgamdx[j]; // Rearranged discrete equation for all points but the boundaries.
        }

        T_next[N-1] = T[N-1]-dt0*dgamdx[N-1]; // declared outside of loop as index gets missed

        double B = fabs(max_array(dTdx, N));
        double C = P*pow(B, P-1);

        // Found via linearisation and Von Neumann Stability Analysis
        dt = pow(dx, 2)/ ((8.0)*C); // 8 is arbritary integer here > 2

        // Array Swap
        double *temp;
        temp = T;
        T = T_next;
        T_next = temp;

        // update time
        t += dt0;

        if (t == next_output_time) {
            for (int j=0; j< N; j++) {
                x = j*dx;
                fprintf(fp,"%g\t%g\t%g\n",t,x,T[j]);
            }
        next_output_time += t_d;
        }
    }

    // Close Files and Free Mallocs
    fclose(fp);
    free(T);
    free(T_next);
    free(dTdx);
    free(gammaa);
    free(dgamdx);

    return 0;
}
