#include "sparsetools.h"

void csr_eliminate_zeros(int n_row, int *row_ptr, int *col_ind, double *val)
{
  int i,j,jj;
  int nnz = 0;
  int row_end = 0;
  double x;
  for (i=0; i<n_row; i++) {
    jj = row_end;
    row_end = row_ptr[i+1];
    while (jj<row_end) {
      j = col_ind[jj];
      x = val[jj];
      if(x != 0) {
        col_ind[nnz] = j;
        val[nnz] = x;
        nnz++;
      }
      jj++;
    }
    row_ptr[i+1] = nnz;
  }
}

int csr_diagonal_optimized(int n_row, int *row_ptr, int *col_ind, double *val)
{
  int i, j, k, nnz, p, *ja, c;
  double *va, v;
  for (i=0;i<n_row;i++) {
    nnz = row_ptr[i+1] -row_ptr[i];
    ja = &(col_ind[row_ptr[i]]);
    va = &(val[row_ptr[i]]);
    p=0;
    while ((p<nnz)&&(ja[p]!=i)) p++;
    if (p==nnz) {
      return i;
    }
    v     = va[0];
    c     = ja[0];
    va[0] = va[p];
    ja[0] = ja[p];
    va[p] = v;
    ja[p] = c;
  }
  for (i=0; i<n_row; i++) {
    for (j=row_ptr[i]+1;j<row_ptr[i+1];j++) {
      for (k=j+1;k<row_ptr[i+1];k++) {
        if (col_ind[k]<col_ind[j]) {
          v          = val[k];
          c          = col_ind[k];
          val[k]     = val[j];
          col_ind[k] = col_ind[j];
          val[j]     = v;
          col_ind[j] = c;
        }
      }
    }
  }
  return 1;
}

int csr_to_profile_pass1(int n_row, int *row_ptr, int *col_ind, int *c)
{
  int i,j,jstart;
  for (i=0; i<n_row; i++) {
    // check for empty rows
    if (row_ptr[i]==row_ptr[i+1]) c[i] = i;

    // find the lowest column
    jstart = col_ind[row_ptr[i]];
    for (j=row_ptr[i]+1; j<row_ptr[i+1]; j++)
      if (col_ind[j]<jstart) jstart = col_ind[j];
    if (jstart>i) return 0; // not symmetric
    c[i] = jstart;
  }
  return 1;
}

int profile_stored_elements(int *c, int n)
{
  int i, nnz = 0;
  for (i=0; i<n; i++) {
    nnz += i - c[i] + 1;
  }
  return nnz;
}

int csr_to_profile_pass2(int n_row, int *row_ptr, int *col_ind, double *val,
                   double **a, int *c)
{
  int i,j;
  for (i=0; i<n_row; i++) {
    for (j=c[i]; j<=i; j++) {
      a[i][j] = 0.;
    }
  }

  for (i=0; i<n_row; i++) {
    for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      if (col_ind[j]<=i) {
        a[i][col_ind[j]] = val[j];
      }
    }
  }
  return 1;
}

void profile_to_csr(double **a, int *c, int n,
    int *row_ptr, int *col_ind, double *val)
{
  int i,j,nnz=0;
  row_ptr[0]=0;
  for (i=0;i<n;i++) {
    // diagonal elements
    val[nnz] = a[i][i];
    col_ind[nnz] = i;
    nnz++;
    // lower elements
    for(j=c[i];j<i;j++) {
      val[nnz] = a[i][j];
      col_ind[nnz] = j;
      nnz++;
    }
    // upper elements
    for(j=i+1;j<n-c[i];j++) {
      val[nnz] = a[j][i];
      col_ind[nnz] = j;
      nnz++;
    }
    row_ptr[i+1]=nnz;
  }
}


