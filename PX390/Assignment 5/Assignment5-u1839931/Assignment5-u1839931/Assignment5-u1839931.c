//
//  main.c
//  Assignment5-u1839931
//
//  Created by Ajan Sittampalam on 15/02/2021.
//

#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <math.h>

#define T(i, j) T[(i) + P.I * (j)]
#define E(i, j) E[(i) + P.I * (j)]
#define B(i, j) B[(i) + P.I * (j)]
#define C(i, j) C[(i) + P.I * (j)]

struct band_mat{
  long ncol;        /* Number of columns in band matrix            */
  long nbrows;      /* Number of rows (bands in original matrix)   */
  long nbands_up;   /* Number of bands above diagonal           */
  long nbands_low;  /* Number of bands below diagonal           */
  double *array;    /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix   */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information         */
};
typedef struct band_mat band_mat;

struct parameter {
  double tf, td, xR, yH, gammaB, Tc, Tw, G, dx, dy, dt; //parameters in specification
  long I, J, K, M; // I = x steps, J = y steps, K = t steps, M = matrix size
};
typedef struct parameter parameter;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

void setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
}

void decv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) -= val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) {
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

/* Setting up the initial A matrix for the so that it can be used to solve for T
   I set the pointers to equal 0 before using decv, because it's a function that takes away the value from itself,
   so by setting it 0 it isn't using some random number. Also at the top left and bottom right the values are adapted
   so that it fits the boundary condition when delta T . n = 0.
*/
void init_A_coefficient(band_mat *A, parameter P) {
  for (int a = 0; a < P.M; ++a) {
    int i = a % P.I, j = a / P.I;
    setv(A, a, a, 1 / (P.dx * P.dx) + 1 / (P.dy * P.dy) + 1 / P.dt);
    setv(A, a, a - P.I + 2 * P.I * (j == 0), 0);
    decv(A, a, a - P.I + 2 * P.I * (j == 0), 1 / (2 * (P.dy * P.dy)));
    setv(A, a, a + P.I - 2 * P.I * (j == P.J - 1), 0);
    decv(A, a, a + P.I - 2 * P.I * (j == P.J - 1), 1 / (2 * (P.dy * P.dy)));
    setv(A, a, a - 1 + 2 * (i == 0), 0);
    decv(A, a, a - 1 + 2 * (i == 0), 1 / (2 * (P.dx * P.dx)));
    setv(A, a, a + 1 - 2 * (i == P.I - 1), 0);
    decv(A, a, a + 1 - 2 * (i == P.I - 1), 1 / (2 * (P.dx * P.dx)));
    }
}
void update(band_mat *A, double *E, double *T, double *B, double *C, parameter P) {
  // Set B for A*T=B
  for (int i = 0; i < P.I; ++i) {
    for (int j = 0; j < P.J; ++j) {
        int j_l = j - 1 + 2 * (j == 0);
        int j_h = j + 1 - 2 * (j == P.J - 1);
        int i_l = i - 1 + 2 * (i == 0);
        int i_r = i + 1 - 2 * (i == P.I - 1);
        B(i, j) = (T(i_r, j) / 2 + T(i_l, j) / 2 - T(i, j)) / (P.dx * P.dx) + (T(i, j_h) / 2 + T(i, j_l) / 2 - T(i, j)) / (P.dy * P.dy) + T(i, j) / P.dt; // using implicit method
        /* While adding the boundary term to B it was multiplied by -2/dx or -2/dy because its coefficient in front of dT/dx or dT/dy
           d^2T/dx^2 = -(2/dx*dx) T_(i) + (2/dx*dx) T_(i+1) -(2/dx) dT/dx_(i)
           d^2T/dy^2 = -(2/dy*dy) T_(j) + (2/dy*dy) T_(j+1) -(2/dy) dT/dy_(j)
        */
        if(i == 0){ //sets the boundary conditions at [0, yh) at x = 0
            if(j == 0){
                B(i, j) += sqrt(2) * P.G / P.dx + sqrt(2) * P.G / P.dy; // set normal at (0,0) to be -(1/sqrt(2))x -(1/sqrt(2))y, so dT/dx = -(1/sqrt(2))*G and dT/dy = -(1/sqrt(2))*G
            } else {
                B(i, j) += P.G * 2 / P.dx; // set normal to be -x so dT/dx = -G
            }
        }
        if(i == P.I - 1){ //sets the boundary at (0, yh] at x = xr
            if(j == P.J - 1){
                B(i, j) += -sqrt(2) * P.G / P.dx - sqrt(2) * P.G / P.dy; // set normal at (xr,yh) to be (1/sqrt(2))x +(1/sqrt(2))y, so dT/dx = (1/sqrt(2))*G and dT/dy = (1/sqrt(2))*G
            } else {
                B(i, j) += -P.G * 2 / P.dx; // set normal to be +x so dT/dx = G
            }
        }
        if(j == 0){ //sets the boundary conditions at (0, xr] at y = 0
            if(i == P.I - 1){
                B(i ,j) += -sqrt(2) * P.G / P.dx + sqrt(2) * P.G / P.dy; // set normal at (xr,0) to be (1/sqrt(2))x -(1/sqrt(2))y, so dT/dx = (1/sqrt(2))*G and dT/dy = -(1/sqrt(2))*G
            } else {
                B(i, j) += P.G * 2 / P.dy; // set normal to be -y so dT/dy = -G
            }
        }
        if(j == P.J - 1){ //sets the boundary conditions at [0, xr) at y = yh
            if(i == 0){
                B(i, j) += +sqrt(2) * P.G / P.dx - sqrt(2) * P.G / P.dy; // set normal at (0,yh) to be -(1/sqrt(2))x +(1/sqrt(2))y, so dT/dx = -(1/sqrt(2))*G and dT/dy = (1/sqrt(2))*G
            } else {
                B(i, j) += -P.G * 2 / P.dy; // set normal to be +y so dT/dy = G
            }
        }
    }
  }
  //Update T by solving A*T=B
  solve_Ax_eq_b(A, T, B);
  for (int i = 0; i < P.I; ++i) {
    for (int j = 0; j < P.J; ++j) {
        C(i, j) = E(i, j); //set C to be the E^(k) so we can calculate E^(k+1)
        E(i, j) = C(i, j)/(1 + (2 * P.dt) * (P.gammaB / 2.0) * (1.0 + tanh((T(i, j) - P.Tc) / P.Tw))); //use implicit method to update E, using 2*dt for delta t, and use T (intermediate)
    }
  }
  for (int i = 0; i < P.I; ++i) {
    for (int j = 0; j < P.J; ++j) {
        T(i, j) -= (E(i , j) - C(i, j))/ 2; //update T to its answer
    }
  }
}

/* Read the parameter given in the input text file */
void read_params(double *tf, double *td, long *I, long *J, double *xr, double *yh, double *gammaB, double *Tc, double *Tw, double *G) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(10!=fscanf(infile,"%lf %lf %ld %ld %lf %lf %lf %lf %lf %lf", tf, td, I, J, xr, yh, gammaB, Tc, Tw, G)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}

/* Read the initial T and E at time 0 using the initialprofile text file */
void read_arrays(double *T, double *E, long N) {
   FILE *infile;
   if(!(infile=fopen("initialprofile.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   for (int m = 0; m<N; ++m){
       fscanf(infile, "%lg %lg", &T[m], &E[m]);
   }
   fclose(infile);
}

int main() {
  parameter P;
  read_params(&P.tf, &P.td, &P.I, &P.J, &P.xR, &P.yH, &P.gammaB, &P.Tc, &P.Tw, &P.G);
  P.M = P.I * P.J;

  if (fabs(round(P.tf/P.td) - P.tf/P.td) <= (1e-08 +(1e-05)*fabs(P.tf/P.td))){ //deciding if tf is an 'exact' multiple of td
      P.K = (long) round(P.tf/ P.td);
  } else {
      P.K = (long) floor(P.tf / P.td);
  }
  P.dx = P.xR/(P.I - 1);
  P.dy = P.yH/(P.J - 1);
  P.dt = P.tf/(2 * P.K);
    
  double *T = malloc(P.M * sizeof(double));
  double *E = malloc(P.M * sizeof(double));
  double *B = malloc(P.M * sizeof(double));
  double *C = malloc(P.M * sizeof(double));

  read_arrays(T, E, P.M);

  band_mat A;
  init_band_mat(&A, P.I, P.I, P.M);
  init_A_coefficient(&A, P);
    
  FILE *output = fopen("output.txt", "w");

  for (int k=0; k<P.K; ++k){
      double t = 2 * P.dt * k;
      for (int i = 0; i < P.I; ++i) {
          for (int j = 0; j < P.J; ++j) {
              fprintf(output,"%lg %lg %lg %lg %lg\n", t, P.dx * i, P.dy * j, T(i, j), E(i, j));
          }
      }
      update(&A, E, T, B, C, P);
  }
  
  fclose(output);
  finalise_band_mat(&A);
  free(T);
  free(E);
  free(B);
  free(C);
  return 0;
}
