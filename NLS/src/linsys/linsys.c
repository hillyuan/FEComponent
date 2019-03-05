#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linsys.h"

int Crout_Profile ( int flag, double **a, double *b, int n, int *c )
{
 int i, j, k;

 if( flag != 3 )	// Reduction of matrix [a]
 {
  for( j = 1; j < n; j++ )
  {
	for( i = c[j]+1; i < j; i++ )
	 for( k = (c[i]>c[j] ? c[i] : c[j]); k < i; k++ )
	  a[j][i] -= a[i][k] * a[j][k];
	for (k=c[j]; k<j; k++)
	{
	 a[j][j] -= a[j][k] / a[k][k] * a[j][k];
	 a[j][k] /= a[k][k];
	}
  }
 }

 if( flag == 2 ) return 1;

 for( i = 1; i < n; i++ )	// Reduction of vector {b}
  for( k = c[i]; k < i; k++ )
	b[i] -= a[i][k] * b[k];

 for( i = 0; i < n; i++ ) b[i] /= a[i][i];

 for( i = n-1; i > 0; i-- )	// Back-substitution
  for( k = c[i]; k < i; k++ )
	b[k] -= a[i][k] * b[i];
 return 1;
}

static void daypx(int n, double a, double *x, double *y)
{
  int i;
  for (i=0;i<n;i++)
    y[i] = x[i]+a*y[i];
}

static void dcopy(int n, double *x, double *y)
{
  memcpy(y, x, sizeof(double)*n);
}

static double ddot(int n, double *x, double *y)
{
  int i;
  double sum = 0;
  for (i=0;i<n;i++) sum += x[i]*y[i];
  return sum;
}

static double dnrm2(int n, double *x)
{
  int i;
  double sum = 0;
  for (i=0;i<n;i++) sum += x[i]*x[i];
  return sqrt(sum);
}

static void daxpy(int n, double a, double *x, double *y)
{
  int i;
  for (i=0;i<n;i++)
    y[i] += a*x[i];
}

static void diagonalAx(double *diag, int n, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++) y[i] = diag[i]*x[i];
}

typedef struct
{
  int *c;
  double **a;
} profile_datatype;

static void profileAx(profile_datatype *s, int n, double *x, double *y)
{
  int i,j;
  for (i=0; i<n; i++)
  {
    y[i] = 0.;
    for (j = s->c[i]; j < i; j++)
    {
      y[i] += s->a[i][j]*x[j];  // inferior elements
      y[j] += s->a[i][j]*x[i];  // superior elements
    }
    y[i] += s->a[i][i]*x[i];    // diagonal element
  }
}

typedef struct
{
  int *row_ptr;
  int *col_ind;
  double *val;
} csr_datatype;

static void csrAx(csr_datatype *s, int n, double *x, double *y)
{
  int i,j,p;
  for(i=0;i<n;i++) {
    y[i]=0;
    for(p=s->row_ptr[i];p<s->row_ptr[i+1];p++) {
      j = s->col_ind[p];
      y[i] += s->val[p]*x[j];
    }
  }
}

static int pcg(void (*Ax)(), void *AxData, void (*Minvx)(), void *MinvxData,
        int n, double *b, double *x, int *iter, int maxiter, double tol)
{
  int k;
  double *r,*p,*q;
  double rr1, rr2, normb, alpha, beta, omega;

  r = (double*) calloc(n, sizeof(double));
  p = (double*) calloc(n, sizeof(double));
  q = (double*) calloc(n, sizeof(double));

  // r = b-A*x;
  (*Ax)(AxData, n, x, r); // A->matvec(x, r); // r = A*x;
  daypx(n, -1.0, b, r); // r = b - r;
  (*Minvx)(MinvxData, n, r, q); // M->solve(r, q); // Mq = r
  dcopy(n, q, p); // p = q;

  rr1 = ddot(n, r, q); // rr1 = r'*r;
  normb = dnrm2(n, q); // normb = norm(b)

  //printf("norm(b_tilde)=%lf, \tnorm(r_tilde)=%lf\n",normb,sqrt(rr1));
  if (sqrt(rr1)<normb*tol) {
    (*iter) = 0;
    free(r);
    free(p);
    free(q);
    return 0;
  }

  for (k=1; k<=maxiter; k++) { //for i=1:maxit

    (*Ax)(AxData, n, p, q); // A->matvec(p, q); // q = A*p;
    omega = ddot(n, p, q); // omega = p'*q;
    alpha = rr1/omega; // alpha = rr1/omega;

    daxpy(n, alpha, p, x); // x = x + alpha*p;
    daxpy(n, -alpha, q, r); // r = r - alpha*q;
    (*Minvx)(MinvxData, n, r, q); // M->solve(r, q); // Mq = r
    rr2 = ddot(n, r, q); // rr2 = r'*q;
    beta = rr2 / rr1; // beta = rr2/rr1;

    daypx(n, beta, q, p); // p = q + beta*p;

    rr1 = rr2; // rr1 = rr2;
    //printf("iter=%d,\t norm(r_tilde)=%lf\n", k, sqrt(rr1));
    if (sqrt(rr1)<normb*tol) { // if (sqrt(rr1)<normb*tol) 
      (*iter) = k;
      free(r);
      free(p);
      free(q);
      return 0;
    }
  }
  
  (*iter) = maxiter; // iter = maxit;
  free(r);
  free(p);
  free(q);
  return 1;
}

int PCG_Profile(double **a, int *c,
        int n, double *b, double *x, int *iter, int maxiter, double tol)
{
  int i,info;
  double *diag = (double*) calloc (n, sizeof(double));
  profile_datatype prof;
  prof.a = a; prof.c = c;

  for (i=0; i<n; i++) diag[i] = 1/a[i][i];

  info = pcg(&profileAx, &prof, &diagonalAx, diag,
         n, b, x, iter, maxiter, tol);

  free(diag);

  return info;
}

int PCG_CSR(int *row_ptr, int *col_ind, double *val,
        int n, double *b, double *x, int *iter, int maxiter, double tol)
{
  int i,j,info;
  double *diag = (double*) calloc (n, sizeof(double));
  csr_datatype csr;
  csr.row_ptr = row_ptr; csr.col_ind = col_ind; csr.val = val;

  for (i=0; i<n; i++) 
    for (j=row_ptr[i]; j<row_ptr[i+1];j++)
      if (i==col_ind[j]) diag[i] = 1/val[j];

  info = pcg(&csrAx, &csr, &diagonalAx, diag,
         n, b, x, iter, maxiter, tol);

  free(diag);

  return info;
}

