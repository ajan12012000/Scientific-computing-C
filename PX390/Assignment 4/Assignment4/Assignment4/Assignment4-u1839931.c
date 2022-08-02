//
//  main.c
//  Assignment4
//
//  Created by Ajan Sittampalam on 19/01/2021.
//

#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>

/* Define structure that holds band matrix information */
struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

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

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
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

void read_params(double *L, long *N, double *v, double *t, double *c) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(5!=fscanf(infile,"%lf %ld %lf %lf %lf",L, N ,v ,t ,c)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}

void read_arrays(double *D, double *S, long N) {
   FILE *infile;
   if(!(infile=fopen("coefficients.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
    char buf[N];
    int i = 0;
    while (fgets(buf, sizeof buf, infile) != NULL) {
        // process line here
        if (sscanf(buf, "%lf %lf", &D[i], &S[i]) < 2) {
            printf("only accepts lines with two elements");
        }
        i+=1;
    }
   fclose(infile);
}

int main() {
  double L, v, t, c;
  long N;
  read_params(&L, &N, &v, &t, &c);
  double *S = malloc(sizeof(double)*N);
  double *D = malloc(sizeof(double)*N);
  read_arrays(D, S, N);
  //for (int i; i<N; i++) {
  //      printf("%lf %lf \n", D[i], S[i]);
  //  }
  //printf("%g %ld %g %g %g \n", L, N, v, t, c);
  //return 0;
  // Grid spacing
  double dx = L/(N-1);
  band_mat pmat;
  band_mat qmat;
    
  long ncols = N;
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */
  long nbands_low = 1;
  long nbands_up  = 1;
  init_band_mat(&pmat, nbands_low, nbands_up, ncols);
  init_band_mat(&qmat, nbands_low, nbands_up, ncols);
  double *p = malloc(sizeof(double)*ncols);
  double *q = malloc(sizeof(double)*ncols);
  double *s = malloc(sizeof(double)*ncols);
  double *tp = malloc(sizeof(double)*ncols);

  long i;
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values.
     Note that boundaries are treated with special cases */
  setv(&pmat,0,0,-(D[0]/(dx*dx)) -(D[1]/(dx*dx)) +(v/dx)-(2*v)/(dx*dx));//using taylor series and dp/dx = vp/D we get d2p/dx2 = (p1-p0-vp0/D0)*2/(dx*dx) at x=0
  setv(&pmat,0,1,((D[1]+D[0])/(dx*dx)) -(v/dx) -t);
  s[0]= -S[0];
  for(i=1; i<ncols; i++) {
    if(i<ncols-1) {
        setv(&pmat,i,i+1, (D[i+1]/(dx*dx)) -(v/dx) -t); //also using implicit finite difference method for better stability
        setv(&pmat,i,i-1, D[i]/(dx*dx));
        setv(&pmat,i,i, -(D[i]/(dx*dx)) -(D[i+1]/(dx*dx)) +(v/dx));
    }
    s[i] = -S[i];
    if(i>=ncols-1) {
        setv(&pmat,i,i-1, 0);
        setv(&pmat,i,i, 1);
        s[i] = c;
    }
  }
  solve_Ax_eq_b(&pmat, p, s); //solve for p
  
  setv(&qmat,0,0,-(D[0]/(dx*dx)) -(D[1]/(dx*dx)) +(v/dx)-(2*v)/(dx*dx));//using taylor series and dq/dx = vq/D we get d2q/dx2 = (q1-q0-vq0/D0)*2/(dx*dx) at x=0
  setv(&qmat,0,1,((D[1]+D[0])/(dx*dx)) -(v/dx));
  tp[0]= -t*p[0];
  for(i=1; i<ncols; i++) {
    if(i<ncols-1) {
        setv(&qmat,i,i+1, (D[i+1]/(dx*dx)) -(v/dx)); //also using implicit finite difference method for better stability
        setv(&qmat,i,i-1, D[i]/(dx*dx));
        setv(&qmat,i,i, -((D[i]+D[i+1])/(dx*dx)) +(v/dx));
    }
    tp[i] = -t*p[i];
    if(i>=ncols-1) {
        setv(&qmat,i,i-1, 0);
        setv(&qmat,i,i, 1);
        tp[i] = c;
    }
  }
  
  solve_Ax_eq_b(&qmat, q, tp); // solve for q
  //printmat(&pmat);
  //printmat(&qmat);
  FILE *out_file;
  out_file = fopen("output.txt", "w");
  if (out_file == NULL)
  {
    printf("Failed to open output file!\n");
    exit(1);
  }
  double x;
  for(i=0; i<ncols; i++) {
    x = i*dx;
    fprintf(out_file, "%g %g %g \n", x, p[i] ,q[i]);
    //printf("%g %g %g \n", x, p[i], q[i]);
  }

  finalise_band_mat(&pmat);
  finalise_band_mat(&qmat);
  free(p);
  free(q);
  free(s);
  free(tp);
  return(0);
}
