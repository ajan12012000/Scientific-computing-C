#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <math.h>

#define ANSI_RESET  "\e[0m"
#define ANSI_BOLD  "\x1B[1m"

#define ANSI_RED  "\e[31m"
#define ANSI_GREEN  "\e[32m"
#define ANSI_YELLOW  "\e[33m"
#define ANSI_MAGENTA  "\e[35m"
#define ANSI_L_CYAN  "\e[96m"

#define FLOAT_CMP_RTOL 1e-05
#define FLOAT_CMP_ATOL 1e-08
#define LOG

#define T_(i, j) T[(i) + P.I * (j)]
#define E_(i, j) E[(i) + P.I * (j)]
#define B_(i, j) B[(i) + P.I * (j)]

#define TIME_GRANULARITY 5

#undef I

/* Parameters */
struct params {
    /* Parameters specified in specification */
    double t_f, t_d, x_R, y_H, gamma_B, T_C, T_w,
            dx, dy, dt;
    /*
     * I: x steps
     * J: y steps
     * K: t steps
     * S: matrix size (I * J)
     */
    int I, J, K, S;
};
typedef struct params params;

/* Band matrix structure */
struct band_mat {
    int ncol;           /* Number of columns in band matrix             */
    int nbrows;         /* Number of rows (bands in original matrix)    */
    int nbands_up;      /* Number of bands above diagonal               */
    int nbands_low;     /* Number of bands below diagonal               */
    double *array;      /* Storage for the matrix in banded format      */
    /* Internal temporary storage for solving inverse problem           */
    int nbrows_inv;     /* Number of rows of inverse matrix             */
    double *array_inv;  /* Store the inverse if this is generated       */
    int *ipiv;          /* Additional inverse information               */
};
typedef struct band_mat band_mat;


/* Main update code */
void update(band_mat *A, double *E, double *T, double *B, params P);

/* Initialise a band matrix of a certain size, allocate memory, and set them parameters. */
int init_band_mat(band_mat *bmat, int nbands_lower, int nbands_upper, int n_columns);

/* Initialiase coefficient matrix A */
void init_A(band_mat *A, params P);

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat);

/* Get a pointer to a location in the band matrix, using the row and column indexes of the full matrix. */
double *getp(band_mat *bmat, int row, int column);

/* Set value to a location in the band matrix, using the row and column indexes of the full matrix. */
void setv(band_mat *bmat, int row, int column, double val);

/* Decrease value at a location in the band matrix, using the row and column indexes of the full matrix. */
void decv(band_mat *bmat, int row, int column, double val);

/* Solve the equation Ax = b for a matrix a stored in band format and x and b real arrays */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b);

/* Print band matrix */
void print_bmat(band_mat *bmat);

/* Print matrix */
void print_mat(double* bmat, params P);

/* Equality check using absolute and relative tolerance for floating-point numbers */
int is_close(double a, double b);

/* Nearest neighbour resampling */
void coarsen(double *input, double *output, int I_new, int J_new);

int main(int argc, const char* argv[]) {
    /* Parameters */
    params P;

    /* Read inputs */
    FILE *input = fopen("input.txt", "r");
    if (input == NULL) {
        fprintf(stderr, "Can't open input file\n");
        return 1;
    }
    fscanf(input, "%lg\n%lg\n%d\n%d\n%lg\n%lg\n%lg\n%lg\n%lg\n",
           &P.t_f, &P.t_d, &P.I, &P.J, &P.x_R, &P.y_H, &P.gamma_B, &P.T_C, &P.T_w);
    fclose(input);

    /* Set remaining parameters */
    P.S = P.I * P.J;

    /*
     * Round total steps if t_f is nicely divisible by t_d otherwise floor
     * to prevent unnecessary time cycles which would never get logged
     */
    if (is_close(round(P.t_f / P.t_d), P.t_f / P.t_d)) {
        P.K = (int) round(P.t_f / P.t_d);
#ifdef LOG
        printf(ANSI_YELLOW "Rounding to %d steps (%g/%g)\n" ANSI_RESET, P.K, P.t_f, P.t_d);
#endif
    } else {
        P.K = (int) floor(P.t_f / P.t_d);
#ifdef LOG
        printf(ANSI_YELLOW "Flooring to %d steps (%g/%g)\n" ANSI_RESET, P.K, P.t_f, P.t_d);
#endif
    }
    P.K *= TIME_GRANULARITY;
#ifdef LOG
    printf(ANSI_YELLOW "TIME_GRANULARITY = %d => %d steps\n" ANSI_RESET, TIME_GRANULARITY, P.K);
#endif

    /* deltas */
    P.dx = P.x_R / (P.I - 1);
    P.dy = P.y_H / (P.J - 1);
    P.dt = P.t_f / P.K;
#ifdef LOG
    printf(ANSI_YELLOW "dx: %g\ndy: %g\ndt: %g\n" ANSI_RESET, P.dx, P.dy, P.dt);
#endif

    /* Reading coefficients for T and E at t = 0 */
    FILE *coefficients = fopen("coefficients.txt", "r");
    if (coefficients == NULL) {
        fprintf(stderr, "Can't open coefficients file\n");
        return 1;
    }

    /*
     * T and E contain current T and E, T_init and E_init preserve initial
     * state for the second simulation on a coarsened grid.
     * B represents an interim calculation used when solving A * T = B.
     */
    double *T = malloc(P.S * sizeof(double)),
           *E = malloc(P.S * sizeof(double)),
           *T_init = malloc(P.S * sizeof(double)),
           *E_init = malloc(P.S * sizeof(double)),
           *B = malloc(P.S * sizeof(double));
    for (int s = 0; s < P.S; ++s) {
        fscanf(coefficients, "%lg %lg", &T[s], &E[s]);
        T_init[s] = T[s];
        E_init[s] = E[s];
    }
    fclose(coefficients);

    /* Create A used in A * T = B */
    band_mat A;
    init_band_mat(&A, P.I, P.I, P.S);
    init_A(&A, P);

#ifdef LOG
    printf(ANSI_GREEN "Coefficient matrix A:" ANSI_RESET "\n");
    print_bmat(&A);
    printf(ANSI_GREEN "Initial T:\n" ANSI_RESET);
    print_mat(T, P);
    printf(ANSI_GREEN "Initial E:\n" ANSI_RESET);
    print_mat(E, P);
#endif

    /* Output files */
    FILE *output = fopen("output.txt", "w");
    FILE *error = fopen("errorest.txt", "w");

    /* Perform simulation */
    for (int k = 0; k < P.K + 1; ++k) {
        double t = P.dt * k;

#ifdef DEBUG
        printf(ANSI_YELLOW "Time: %.2g (%d/%d)\n" ANSI_RESET, t, k, P.K);
        if (!(k % TIME_GRANULARITY)) {
            printf(ANSI_L_CYAN " Logging output\n" ANSI_RESET);
        }
        printf(ANSI_GREEN " T:" ANSI_RESET "\n");
        print_mat(T, P);
        printf(ANSI_GREEN " E:" ANSI_RESET "\n");
        print_mat(E, P);
#endif

        /* Log simulation state */
        if (!(k % TIME_GRANULARITY)) {
            for (int i = 0; i < P.I; ++i) {
                for (int j = 0; j < P.J; ++j) {
                    fprintf(output, "%lg %lg %lg %lg %lg\n", t, P.dx * i, P.dy * j, T_(i, j), E_(i, j));
                }
            }
        }

        /* Update simulation state */
        update(&A, E, T, B, P);

#ifdef DEBUG
        printf(ANSI_MAGENTA " B:" ANSI_RESET "\n");
        print_mat(B, P);
#endif
    }

#ifdef LOG
    printf(ANSI_GREEN "Final T:\n" ANSI_RESET);
    print_mat(T, P);
    printf(ANSI_GREEN "Final E:\n" ANSI_RESET);
    print_mat(E, P);
#endif

    /* Update parameters for coarse grid */
    P.I = (P.I + 1) / 2;
    P.J = (P.J + 1) / 2;
    P.S = P.I * P.J;
    P.dx = P.x_R / (P.I - 1);
    P.dy = P.y_H / (P.J - 1);

    /* Resize A for coarse grid */
    finalise_band_mat(&A);
    init_band_mat(&A, P.I, P.I, P.S);
    init_A(&A, P);

    /*
     * Resize initial grids using nearest neighbours.
     * Tc and Ec contain current T and E, Bc is an interim calculation used
     * to solve A * Tc = Bc.
     * Cc is used to compare results with non-coarsened simulation.
     */
    double *Tc = malloc(P.S * sizeof(double)),
           *Ec = malloc(P.S * sizeof(double)),
           *Cc = malloc(P.S * sizeof(double)),
           *Bc = malloc(P.S * sizeof(double));

    coarsen(T_init, Tc, P.I, P.J);
    coarsen(E_init, Ec, P.I, P.J);

#ifdef LOG
    printf(ANSI_GREEN "Coefficient matrix A (coarse):" ANSI_RESET "\n");
    print_bmat(&A);
    printf(ANSI_GREEN "Initial T (coarse):\n" ANSI_RESET);
    print_mat(Tc, P);
    printf(ANSI_GREEN "Initial E (coarse):\n" ANSI_RESET);
    print_mat(Ec, P);
#endif

    /* Perform coarse simulation */
    for (int k = 0; k < P.K + 1; ++k) {
        update(&A, Ec, Tc, Bc, P);
    }

#ifdef LOG
    printf(ANSI_GREEN "Final coarse T:\n" ANSI_RESET);
    print_mat(Tc, P);
    printf(ANSI_GREEN "Final coarse E:\n" ANSI_RESET);
    print_mat(Ec, P);
#endif

    /* Make comparisons */
    coarsen(T, Cc, P.I, P.J);
    double T_error = 0;
    for(int s = 0; s < P.S; ++s) {
        T_error += pow(Tc[s] - Cc[s], 2);
    }
    T_error /= P.S;
    fprintf(error, "%g\n", T_error);

#ifdef LOG
    printf(ANSI_GREEN "Resized T:\n" ANSI_RESET);
    print_mat(Cc, P);
    printf(ANSI_YELLOW "Final error: " ANSI_RESET "%g\n", T_error);
#endif

    /* Cleanup */
    fclose(output);
    fclose(error);
    finalise_band_mat(&A);
    free(T);
    free(T_init);
    free(Tc);
    free(E);
    free(E_init);
    free(Ec);
    free(B);
    free(Bc);
    free(Cc);

    return 0;
}

/* Main update code */
void update(band_mat *A, double *E, double *T, double *B, params P) {
    /* Update E and calculate B for A*T=B */
    for (int i = 0; i < P.I; ++i) {
        for (int j = 0; j < P.J; ++j) {
            int j_l = j - 1 + 2 * (j == 0),
                    j_h = j + 1 - 2 * (j == P.J - 1),
                    i_l = i - 1 + 2 * (i == 0),
                    i_r = i + 1 - 2 * (i == P.I - 1);
            double dEdt = -E_(i, j) * (P.gamma_B / 2.0) * (1.0 + tanh((T_(i, j) - P.T_C) / P.T_w));
            E_(i, j) += dEdt * P.dt;
            B_(i, j) = (T_(i_r, j) / 2 + T_(i_l, j) / 2 - T_(i, j)) / pow(P.dx, 2) +
                       (T_(i, j_h) / 2 + T_(i, j_l) / 2 - T_(i, j)) / pow(P.dy, 2) +
                       T_(i, j) / P.dt - dEdt;
        }
    }

    /* Update T by solving A*T=B */
    solve_Ax_eq_b(A, T, B);
}

/* Initialise a band matrix of a certain size, allocate memory, and set the parameters.  */
int init_band_mat(band_mat *bmat, int nbands_lower, int nbands_upper, int n_columns) {
    bmat->nbrows     = nbands_lower + nbands_upper + 1;
    bmat->ncol       = n_columns;
    bmat->nbands_up  = nbands_upper;
    bmat->nbands_low = nbands_lower;
    bmat->array      = (double *) malloc(sizeof(double) * bmat->nbrows * bmat->ncol);
    bmat->nbrows_inv = bmat->nbands_up * 2 + bmat->nbands_low + 1;
    bmat->array_inv  = (double *) malloc(sizeof(double) * (bmat->nbrows + bmat->nbands_low) * bmat->ncol);
    bmat->ipiv       = (int *) malloc(sizeof(int) * bmat->ncol);
    if (bmat->array == NULL || bmat->array_inv == NULL) {
        return 0;
    }
    /* Initialise array to zero */
    for (int i = 0; i < bmat->nbrows * bmat->ncol; i++) {
        bmat->array[i] = 0.0;
    }
    return 1;
}

/* Initialiase coefficient matrix A */
void init_A(band_mat *A, params P) {
    for (int s = 0; s < P.S; ++s) {
        int i = s % P.I,
                j = s / P.I;
        setv(A, s, s, 1 / pow(P.dx, 2) + 1 / pow(P.dy, 2) + 1 / P.dt);
        decv(A, s, s - P.I + 2 * P.I * (j == 0), 1 / (2 * pow(P.dy, 2)));
        decv(A, s, s + P.I - 2 * P.I * (j == P.J - 1), 1 / (2 * pow(P.dy, 2)));
        decv(A, s, s - 1 + 2 * (i == 0), 1 / (2 * pow(P.dx, 2)));
        decv(A, s, s + 1 - 2 * (i == P.I - 1), 1 / (2 * pow(P.dx, 2)));
    }
}

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using he row and column indexes of the full matrix. */
double *getp(band_mat *bmat, int row, int column) {
    int bandno = bmat->nbands_up + row - column;
    if (row < 0 || column < 0 || row >= bmat->ncol || column >= bmat->ncol) {
        fprintf(stderr, "Indexes out of matrix bounds:\n row %d\n col %d\n", row, column);
        exit(1);
    }
    return &bmat->array[bmat->nbrows * column + bandno];
}

/* Set value to a location in the band matrix, using the row and column indexes of the full matrix. */
void setv(band_mat *bmat, int row, int column, double val) {
    double *valr = getp(bmat, row, column);
    int bandno = bmat->nbands_up + row - column;
    if (bandno < 0 || bandno >= bmat->nbrows) {
        printf(ANSI_RED "Setting (%d, %d) out of band:\n band %d\n" ANSI_RESET,
               row, column, bandno);
        exit(1);
    }
    *valr = val;
}

/* Decrease value at a location in the band matrix, using the row and column indexes of the full matrix. */
void decv(band_mat *bmat, int row, int column, double val) {
    double *valr = getp(bmat, row, column);
    int bandno = bmat->nbands_up + row - column;
    if (bandno < 0 || bandno >= bmat->nbrows) {
        printf(ANSI_RED "Increasing (%d, %d) out of band:\n band %d\n" ANSI_RESET,
               row, column, bandno);
        exit(1);
    }
    *valr -= val;
}

/* Solve the equation Ax = b for a matrix a stored in band format and x and b real arrays */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
    /* Copy bmat array into the temporary store */
    for(int i = 0; i < bmat->ncol; i++) {
        for (int bandno = 0; bandno < bmat->nbrows; bandno++) {
            bmat->array_inv[bmat->nbrows_inv * i + (bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows * i + bandno];
        }
        x[i] = b[i];
    }

    int nrhs = 1;
    int ldab = bmat->nbands_low * 2 + bmat->nbands_up + 1;
    int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
    return info;
}

/* Print band matrix */
void print_bmat(band_mat *bmat) {

    /* Column header */
    printf("       ");
    for (int i = 0; i < bmat->ncol; ++i) {
        printf("%11d ", i);
    }
    printf("\n");
    printf("     ???");
    for (int i = 0; i < bmat->ncol; ++i) {
        printf("????????????????????????????????????");
    }
    printf("??????\n");

    /* Rows */
    for (int i = 0; i < bmat->ncol; ++i) {
        printf("%4d ??? ", i);
        for (int j = 0; j < bmat->ncol; ++j) {
            int bandno = bmat->nbands_up + i - j;
            if (bandno < 0 || bandno >= bmat->nbrows) {
                printf(ANSI_RED "%11.4g " ANSI_RESET, 0.0);
            } else if (i == j) {
                printf(ANSI_BOLD "%11.4g " ANSI_RESET, *getp(bmat, i, j));
            } else {
                printf("%11.4g ", *getp(bmat, i, j));
            }
        }
        printf("???\n");
    }
    printf("     ???");
    for (int i = 0; i < bmat->ncol; ++i) {
        printf("????????????????????????????????????");
    }
    printf("??????\n");
}

/* Print matrix */
void print_mat(double* mat, params P) {

    /* Column header */
    printf("       ");
    for (int i = 0; i < P.I; ++i) {
        printf("%11d ", i);
    }
    printf("\n");
    printf("     ???");
    for (int i = 0; i < P.I; ++i) {
        printf("????????????????????????????????????");
    }
    printf("??????\n");

    /* Rows */
    for (int i = 0; i < P.J; ++i) {
        printf("%4d ??? ", i);
        for (int j = 0; j < P.J; ++j) {
            printf("%11.4g ", mat[i + P.I * j]);
        }
        printf("???\n");
    }
    printf("     ???");
    for (int i = 0; i < P.J; ++i) {
        printf("????????????????????????????????????");
    }
    printf("??????\n");
}

/* Equality check using absolute and relative tolerance for floating-point numbers */
int is_close(double a, double b) {
    return fabs(a - b) <= (FLOAT_CMP_ATOL + FLOAT_CMP_RTOL * fabs(b));
}

/* Nearest neighbour resampling */
void coarsen(double *input, double *output, int I_new, int J_new) {
    int I_old = 2 * I_new - 1;
    for (int i = 0; i < I_new; ++i) {
        for (int j = 0; j < J_new; ++j) {
            output[i + I_new * j] = input[2 * i + I_old * 2 * j];
        }
    }
}
