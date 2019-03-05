#ifndef _SPARSETOOLS_H_
#define _SPARSETOOLS_H_

#ifdef __cplusplus
extern "C" {
#endif

void csr_eliminate_zeros(int n_row, int *row_ptr, int *col_ind, double *val);

int csr_diagonal_optimized(int n_row, int *row_ptr,
    int *col_ind, double *val);

int csr_to_profile_pass1(int n_row, int *row_ptr, int *col_ind, int *c);

int csr_to_profile_pass2(int n_row, int *row_ptr, int *col_ind, double *val,
                   double **a, int *c);

int profile_stored_elements(int *c, int n);

void profile_to_csr(double **a, int *c, int n,
    int *row_ptr, int *col_ind, double *val);

#ifdef __cplusplus
}
#endif

#endif
