// --------------------------------------------------------------------------
//
// LinSys.cpp - Linear System class.
//
// --------------------------------------------------------------------------

#include "linsys/linearsystem.h"
#include "mempack/mempack.h"
#include "mathpack/mathpack.h"

// --------------------------------------------------------------------------
// Public functions:

// ============================= cLinSys ===================================

cLinSys :: cLinSys ( ) : _iNumEq(0)
{
}

// ============================= ~cLinSys ==================================

cLinSys :: ~cLinSys ( void )
{
}

// ============================= Solve ==================================

void cLinSys :: Solve ( double *b )
{
  // x = 0 (start solution)
  double *x = (double *) MemAlloc( _iNumEq, sizeof(double) );
  MathVecZero( _iNumEq, x );
  Solve ( b, x );
  MathVecAssign( _iNumEq, b, x ); // b = x
  MemFree( x );
}


