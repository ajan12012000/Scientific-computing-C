/*****************************************************************
 * This programme defines functions used to generate band
 * matrix for problems in which a solution at some grid point
 * k is coupled to a solution at a number of adjasont points.
 * The file also includes a main() function with an
 * example of how to solve a stationary 2D heat equation with a
 * constant source, for two different sources.
 * 
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o heq2d heateq2d.c -lmkl -liomp5 -lm 
 * Run: ./bandu
 * ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>

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

/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P) {
  return (j>=0)&&(j<J)&&(p>=0)&&(p<P);
}


/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long k:  The Y grid point index
   long P:  The number of Y points.
*/
long indx( long j, long p, long P) {
  return j*P + p;
}

/* Return the 2D point corresponding to a particular 1D grid index */
void gridp(long indx, long P, long *j, long *p) {
  *j = indx%P;
  *p = indx - (*j)*P;
}

/* An example of how to use the band matrix routines to solve a PDE:
   The equation solved is related to the steady state solution of the heat 
   diffusion equation.   
*/
int main() {
  band_mat bmat;
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */

  long nx = 200;
  long ny = 200;
  /* Total size of problem is number of grid points on 2D plane */
  long ncols = nx*ny; 
  /* The matrix is block tridiagonal: one block on either
  side of the blocks of diagonals */ 
  long nbands_low = 2*ny-1; 
  long nbands_up  = nbands_low;

  init_band_mat(&bmat, nbands_low, nbands_up, ncols);

  /*Storage for the source term and the solution */
  double *x = malloc(sizeof(double)*ncols);
  double *b = malloc(sizeof(double)*ncols);

  double xlen = 1.0;
  double ylen = 1.0;
  /* For zero boundary conditions, the number of intervals is one
     more than the number of unknowns */
  double dx = xlen/(nx+1);
  double dy = ylen/(ny+1);
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values 
     note boundaries treated with special cases           */
  long j, p;
  long index;
  long unknown_indx;
  /* chanage these values to set surce and coupling
     coefficients at particular grid points. Set j0=p0=-1
     for Dirichlet boundary conditions.
     Set j0=nx/2, p0=ny/2 to set a source in the middle
     of the domain. Coupling coefficient for this
     point is then set to 1 on the diagonal. */
  long j0=-1, p0=-1;

  /* Main loops start here...*/
  for(j=0; j<nx; j++) {
    for(p=0; p<ny; p++) {
      unknown_indx = indx(j,p,ny);
      if(j!=j0 || p!=p0){
        if(is_valid(j-1,p,nx,ny)) {
          setv(&bmat,unknown_indx,indx(j-1,p,ny),1.0/(dx*dx));
        }
        if(is_valid(j+1,p,nx,ny)) {
          setv(&bmat,unknown_indx,indx(j+1,p,ny),1.0/(dx*dx));
        }
        if(is_valid(j,p-1,nx,ny)) {
          setv(&bmat,unknown_indx,indx(j,p-1,ny),1.0/(dy*dy));
        }
        if(is_valid(j,p+1,nx,ny)) {
          setv(&bmat,unknown_indx,indx(j,p+1,ny),1.0/(dy*dy));
        }
        setv(&bmat,unknown_indx,indx(j,p,ny),-2.0/(dx*dx) - 2.0/(dy*dy));
        /* Uniform source term in heat equation */
        b[indx(j,p,ny)] = 1.0;   
      }
      else{
        setv(&bmat,unknown_indx, unknown_indx, 1.0);
        b[indx(j,p,ny)] = 10.0;
      }
    }
  }
  
  /* Output sim parameters to a file */
  FILE *fp = fopen("heq2d_params","w");
  fprintf(fp, "%ld %ld \n", nx, ny);

  /* Solve the matrix equation */
  solve_Ax_eq_b(&bmat, x, b);
    
  for(j=0; j<nx; j++) {
    for(p=0; p<ny; p++) {   
      index = indx(j,p,ny);
      printf("%ld %ld %g %g \n",j,p,x[index],b[index]);
    }
  }
  /* Free memory */
  finalise_band_mat(&bmat);
  free(x);
  free(b);
  return 0;
}
