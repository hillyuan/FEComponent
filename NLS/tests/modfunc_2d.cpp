// ----------------------------------------------------------------------
//
// ModFunc_MT.cpp - Function from CEE 598 Midterm
//
// ----------------------------------------------------------------------

#include <math.h>
#include "modfunc_2d.h"

cModelFunction_MT :: cModelFunction_MT ( char *filename ) : cModel ( )
{
	_iNumEq = _iNumGra = 2;
}

cModelFunction_MT :: ~cModelFunction_MT ( void )
{
}

void cModelFunction_MT :: Init ( void )
{
	InitFile( );
}

void cModelFunction_MT :: InternalVector ( double *u, double *f )
{
 	f[0] = 10.0 * u[0] + 0.4 * pow(u[1], 3) - 5.0 * pow(u[1], 2);
	f[1] = 0.4 * pow(u[0], 3) - 3.0 * pow(u[0],2) + 10 * u[1];
}

void cModelFunction_MT :: Reference ( double *f )
{
	f[0] = 40.0;
	f[1] = 15.0;
}

void cModelFunction_MT :: TangentMatrix ( double *u, cLinSys *kt )
{
   double k11 = 10.0;
   double k12 = 1.2 * u[1] * u[1] - 10.0 * u[1];
   double k21 = 1.2 * u[0] * u[0] -  6.0 * u[0];
   double k22 = 10.0;

   kt->FillA(0, 0, k11);
   kt->FillA(0, 1, k12);
   kt->FillA(1, 0, k21);
   kt->FillA(1, 1, k22);

}

