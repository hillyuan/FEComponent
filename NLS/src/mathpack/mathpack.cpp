// --------------------------------------------------------------------------
//
// MathPack.cpp - This file contains the array functions.
//
// --------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ============================ MathMatCrout ================================

void MathMatCrout ( int flag, double **a, double *b, int n, int *c )
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

 if( flag == 2 ) return;

 for( i = 1; i < n; i++ )	// Reduction of vector {b}
  for( k = c[i]; k < i; k++ )
	b[i] -= a[i][k] * b[k];

 for( i = 0; i < n; i++ ) b[i] /= a[i][i];

 for( i = n-1; i > 0; i-- )	// Back-substitution
  for( k = c[i]; k < i; k++ )
	b[k] -= a[i][k] * b[i];
}

// ============================= MathVecDot =================================

double MathVecDot ( int n, double *w, double *v )
{
 double dot = 0.0;

 for( int i = 0; i < n; i++ ) dot += w[i] * v[i];

 return( dot );
}

// =========================== MathVecScale =================================

void MathVecScale ( int n, double e, double *v )
{
 for( int i = 0; i < n; i++ ) v[i] *= e;
}

// ========================== MathVecMaxNorm ================================

double MathVecMaxNorm ( int n, double *v )
{
 double MaxNorm = 0.0;

 for( int i = 0; i < n; i++ )
  if( fabs( v[i] ) > MaxNorm )
	MaxNorm = fabs( v[i] );

 return( MaxNorm );
}

// ============================ MathVecNorm =================================

double MathVecNorm ( int n, double *v )
{
 return( sqrt( MathVecDot( n, v, v ) ) );
}

// =========================== MathVecAssign ================================

void MathVecAssign ( int n, double *a, double *b )
{
 for( int i = 0; i < n; i++ ) a[i] = b[i];
}

// =========================== MathVecAssign1 ===============================

void MathVecAssign1 ( int n, double *a, double s, double *b )
{
 for( int i = 0; i < n; i++ ) a[i] = s * b[i];
}

// =========================== MathVecAssign2 ===============================

void MathVecAssign2 ( int n, double *a, double s, double *b, double *c )
{
 for( int i = 0; i < n; i++ ) a[i] = s * b[i] + c[i];
}

// ============================ MathVecZero =================================

void MathVecZero ( int n, double *a )
{
 for( int i = 0; i < n; i++ ) a[i] = 0.0;
}

// ============================= MathVecAdd =================================

void MathVecAdd ( int n, double *a, double *b )
{
 for( int i = 0; i < n; i++ ) a[i] += b[i];
}

// ============================= MathVecAdd1 ================================

void MathVecAdd1 ( int n, double *a, double s, double *b )
{
 for( int i = 0; i < n; i++ ) a[i] += s * b[i];
}

// ============================= MathVecAdd2 ================================

void MathVecAdd2 ( int n, double *a, double s, double *b, double *c )
{
 for( int i = 0; i < n; i++ ) a[i] += s * b[i] + c[i];
}

// ============================= MathVecSub =================================

void MathVecSub ( int n, double *a, double *b )
{
 for( int i = 0; i < n; i++ ) a[i] -= b[i];
}

