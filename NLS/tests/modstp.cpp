// ----------------------------------------------------------------------
//
// ModFunc_nobre.cpp
//
// ----------------------------------------------------------------------

#include <math.h>
#include "modstp.h"

cModelFunction_Nobre :: cModelFunction_Nobre ( char *filename ) : cModel ( )
{
 _iNumEq = _iNumGra = 1;
}

cModelFunction_Nobre :: ~cModelFunction_Nobre ( void )
{
}

void cModelFunction_Nobre :: Init ( void )
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

void cModelFunction_Nobre :: InternalVector ( double *u, double *f )
{
 double x = u[0];
 f[0] = (sqrt(2.0)/4)*(2*x - 3*pow(x,2) + pow(x,3));
}

void cModelFunction_Nobre :: Reference ( double *f )
{
 f[0] = 1;
}

void cModelFunction_Nobre :: TangentMatrix ( double *u, cLinSys *kt )
{
 double x = u[0];
 double k = (sqrt(2.0)/4)*(2 - 6*x + 3*pow(x,2));
 kt->Zero();
 kt->AddA(0,0,k);
}