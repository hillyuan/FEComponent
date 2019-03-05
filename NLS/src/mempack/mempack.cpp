// ----------------------------------------------------------------------
//
// MemPack.cpp - This file contains the allocation functions.
//
// ----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

// ============================= MemAlloc ===============================

void *MemAlloc ( int n, int s )
{
 void *v = NULL;

 v = calloc ( n, s );

 if( !v )
 {
  printf( "\n\t ### Allocation error (MemAlloc) !!! ###\n\n" );
  exit( -1 );
 }

 return( v );
}

// ============================= MemFree ================================

void MemFree ( void *v )
{
 free( v );
}

// ========================= MemProfileAlloc ============================

double **MemProfileAlloc ( int n, int *profile )
{
 double *v;
 double **address = NULL;

 address = (double **) calloc ( n, sizeof(double *) );

 if( !address )
 {
  printf( "\n\t ### Allocation error (MemProfileAlloc) !!! ###\n\n" );
  exit( -1 );
 }

 for( int i = 0; i < n; i++ )
 {
  v = NULL;
  v = (double *) calloc ( i - profile[i] + 1, sizeof(double) );
  if( !v )
  {
	printf( "\n\t ### Allocation error (MemProfileAlloc) !!! ###\n\n" );
	exit( -1 );
  }
  address[i] = v - profile[i];
 }

 return( address );
}

// ========================== MemProfileFree ============================

void MemProfileFree ( double **m, int n, int *profile )
{
 for( int i = 0; i < n; i++ ) free( m[i]+profile[i] );
 free( m );
}



