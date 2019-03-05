// ----------------------------------------------------------------------
//
// Model.cpp - Metodos gerais da classe modelo.
//
// ----------------------------------------------------------------------

#include <stdio.h>
#include "mempack/mempack.h"
#include "model/model.h"

cModel :: ~cModel ( void )
{
 for( int i = 0; i < _iNumGra; i++ ) fclose( _afOut[i] );
 MemFree( _afOut );

 _iNumEq = 0;
}

int cModel :: NumEpsEq ( void )
{
 return( _iNumEpsEq );
}

int cModel :: NumEq ( void )
{
 return( _iNumEq );
}

int cModel :: Profile ( int *piProfile )
{
 if (piProfile==0) return( 0 );
 // Generate profile vector (full matrix)
 for( int i = 0; i < _iNumEq; i++ ) piProfile[i] = 0;
 return( 1 );
}

int cModel :: SparsityPattern ( int **piSparsity )
{
 if (piSparsity==0) return( 0 );
 // Generate sparticity matrix (full matrix)
 for( int i = 0; i < _iNumEq; i++ ) {
   for( int j = 0; j < _iNumEq; j++ )
     piSparsity[i][j] = 1;
 }
 return( 1 );
}

void cModel :: Convergence ( double dFactor, double *pdSol )
{
 for( int i = 0; i < _iNumGra; i++ )
  fprintf( _afOut[i], "%f   %f\n", pdSol[i], dFactor );
}

void cModel :: InitFile ( void )
{
 _afOut = (FILE **) MemAlloc( _iNumEq, sizeof(FILE *) );
 for( int i = 0; i < _iNumEq; i++ )
 {
  char fn[50];
  sprintf( fn, "nls%d.out", i );
  _afOut[i] = fopen( fn, "w" );
  fprintf( _afOut[i], "%f   %f\n", 0.0, 0.0 );
 }
}