#include <stdlib.h>
#include <stdio.h>

void init_vect(double *M, int N)
{
  int j;
  // random numbers in the range [0.5, 1.5)

  for (j=0; j<N; j++)
    M[j] = 0.5 + (double)rand()/RAND_MAX;  
}

void init_mat(double *M, int N)
{
  int j, k;
  // random numbers in the range [0.5, 1.5)

  for (k=0; k<N; k++) 
    for (j=0; j<N; j++)
      M[k*N+j] = 0.5 + (double)rand()/RAND_MAX;
}

void zero_mat(double *M, int N)
{
  int j, k;

  for (k=0; k<N; k++) 
    for (j=0; j<N; j++)
      M[k*N+j] = 0.0;
}


double checksum_vect ( double *const c, int N )
{
  int i;
  double S= 0.0;

  for (i=0; i<N; i++)
    S += c[i];
  return S;
}

double checksum_mat ( double *const c, int N )
{
  int i, j;
  double S= 0.0;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      S += c[i*N+j];
  return S;
}

void f1_mat ( double *const x, double *const y, double *restrict a, int N )
{
  int i, j;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i*N+j] = x[i] * y[j];
}

void f2_mat ( double *const x, double *const y, double *restrict a, int N )
{
  int i, j;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i+N*j] = x[i] * y[j];
}


void f1_vect ( double *x, double r, int N )
{
  int i;

  for (i=0; i<N; i++)
    x[i] = x[i] / r;
}

void do_op ( double *const a, double *const b, double *const bt, double *restrict c, int N )
{
  int i, j, k;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      for (k=0; k<N; k++) {
        c[i*N + j] += a [i*N + k] * b[j*N + k];
        c[i*N + j] += bt[i*N + k] * a[j*N + k];
      }
}

void mat_transpose (double *const M, double *Mt, int N)
{
  int j, k;

  for (k=0; k<N; k++) 
    for (j=0; j<N; j++) {
      Mt[k*N + j] = M[j*N + k];
    }
}

//////// MAIN ////////////
int main (int argc, char **argv)
{
  int N=2000;
  double *A, *B, *Bt, *C, *X, *Y, R;

  if (argc>1) {  N  = atoll(argv[1]); }
  if (N<1 || N>20000) {
     printf("input parameter: N (1-20000)\n");
     return 0;
  }

  // dynamic allocation of 2-D matrices
  A  = (double *) malloc ( N*N*sizeof(double));
  B  = (double *) malloc ( N*N*sizeof(double));
  Bt = (double *) malloc ( N*N*sizeof(double));
  C  = (double *) malloc ( N*N*sizeof(double));

  // Dynamic allocation of vectors
  X = (double *) malloc ( N*sizeof(double));
  Y = (double *) malloc ( N*sizeof(double));
  
  // initial seed for random generation
  srand(1);

  // Initialize input data with random data
  init_vect (X, N);
  init_vect (Y, N);
  R=0.0;

  // Main computation
  f1_mat        (X, Y, A, N);
  f2_mat        (X, Y, B, N);
  R +=          checksum_vect (Y, N);
  f1_vect       (X, R, N);
  R +=          checksum_vect (X, N);
  zero_mat      (C, N);
  mat_transpose (B, Bt, N);
  do_op         (A, B, Bt, C, N);
  R +=          checksum_mat(C, N);

  // Output a single value
  printf("Final Result  (N= %d ) = %e\n", N, R);

  free (A); free (B); free (Bt); free (C); free (X); free (Y);
  return 0;
}

