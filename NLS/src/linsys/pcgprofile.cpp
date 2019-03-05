// --------------------------------------------------------------------------
//
// PCGProfile.cpp - PCG Linear System class.
//
// --------------------------------------------------------------------------

#include "linsys/pcgprofile.h"
#include "mempack/mempack.h"
#include "mathpack/mathpack.h"
#include "model/model.h"
#include <linsys.h>

// --------------------------------------------------------------------------
// Public functions:

// ============================= cLinSys ===================================

cPCGProfile :: cPCGProfile (int iMaxIter, double dTol ) : _iMaxIter(iMaxIter), _dTol(dTol)
{
}

// ============================= ~cLinSys ==================================

cPCGProfile :: ~cPCGProfile ( void )
{
}

void cPCGProfile :: Solve ( double *b )
{
  // x = 0 (start solution)
  double *x = (double *) MemAlloc( _iNumEq, sizeof(double) );
  MathVecZero( _iNumEq, x );
  Solve ( b, x );
  MathVecAssign( _iNumEq, b, x ); // b = x
  MemFree( x );
}

void cPCGProfile :: Solve ( double *b, double *x )
{
  int iter;
  PCG_Profile(_paA, _piProfile,
      _iNumEq, b, x, &iter, _iMaxIter, _dTol);
  printf("PCG_Profile - iterations = %d \n", iter);
}


