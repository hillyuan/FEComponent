// ----------------------------------------------------------------------
//
// ModFunc.cpp - Funcao apresentada no artigo do Hong Chen.
//
// ----------------------------------------------------------------------

#include <math.h>
#include "modfunc_1d.h"

cModelFunction :: cModelFunction ( char *filename ) : cModel ( )
{
 _iNumEq = _iNumGra = 1;
}

cModelFunction :: ~cModelFunction ( void )
{
}

void cModelFunction :: Init ( void )
{
 //InitFile( );
 _afOut = (FILE **) MemAlloc( _iNumEq, sizeof(FILE *) );
 for( int i = 0; i < _iNumEq; i++ )
 {
  char fn[50];
  sprintf( fn, "function_%d.out", i );
  _afOut[i] = fopen( fn, "w" );
  fprintf( _afOut[i], "%f   %f\n", 0.0, 0.0 );
 }
}

void cModelFunction :: InternalVector ( double *u, double *f )
{
 double x = u[0] - 1.0;
 f[0] = (x < 0 ? 3 : -3) * pow( fabs(x), 1.0/3.0 ) + 4 * x + 1;
}

void cModelFunction :: Reference ( double *f )
{
 f[0] = 1;
}

void cModelFunction :: TangentMatrix ( double *u, cLinSys *kt )
{
 double x = u[0] - 1.0; //why??
 double k = x == 0 ? -1e32 : -1 / pow( fabs(x), 2/3.0 ) + 4; //why??
 kt->Zero();
 kt->AddA(0,0,k);
}