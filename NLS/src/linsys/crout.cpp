// --------------------------------------------------------------------------
//
// Crout.cpp - Crout Linear System class.
//
// --------------------------------------------------------------------------

#include "linsys/crout.h"
#include "mempack/mempack.h"
#include "mathpack/mathpack.h"
#include "model/model.h"
#include <math.h>
#include <linsys.h>

// --------------------------------------------------------------------------
// Public functions:

// ============================= cLinSys ===================================

cCroutProfile :: cCroutProfile ( ) : _piProfile(0),  _paA(0), _isFactorized(0)
{
}

// ============================= ~cLinSys ==================================

cCroutProfile :: ~cCroutProfile ( void )
{
 if (_paA!=0) MemProfileFree( _paA, _iNumEq, _piProfile );
 if (_piProfile!=0) MemFree( _piProfile );
}

void cCroutProfile :: Init ( cModel *pcModel )
{
 // Release previous allocated values
 if (_paA!=0) MemProfileFree( _paA, _iNumEq, _piProfile );
 if (_piProfile!=0) MemFree( _piProfile );

 _iNumEq      = pcModel->NumEq( );
 _piProfile    = (int *) MemAlloc( _iNumEq, sizeof(int) );
 pcModel->Profile( _piProfile );
 _paA = MemProfileAlloc( _iNumEq, _piProfile );
}

void cCroutProfile :: Zero ( )
{
 int i,j;
 for(i=0;i<_iNumEq;i++) {
    for(j=_piProfile[i];j<=i;j++) {
       _paA[i][j] = 0.;
    }
 }
 _isFactorized = 0;
}

void cCroutProfile :: AddA ( int i, int j, double k )
{
 if (i <= j) {
     // in upper triangle or on main diagonal
     if (i >= _piProfile[j]) {
         _paA[j][i] += k;
     } else {
         // element not defined;
     }
 } else {
     // in lower triangle
     if (j >= _piProfile[i]) {
         _paA[i][j] += k;
     } else {
         // element not defined;
     }
 }
 _isFactorized = 0;
}

void cCroutProfile :: FillA ( int i, int j, double k )
{
 _paA[i][j] = k;

 _isFactorized = 0;
}

void cCroutProfile :: Solve ( double *_pdDx )
{
 if( _isFactorized ) {
   Crout_Profile( 3, _paA, _pdDx, _iNumEq, _piProfile );
 } else {
   if (_pdDx==0)
     Crout_Profile( 2, _paA, _pdDx, _iNumEq, _piProfile );
   else
     Crout_Profile( 1, _paA, _pdDx, _iNumEq, _piProfile );
   _isFactorized = 1;
 }
 //printf("TOP_Crout_Profile\n");
}


void cCroutProfile :: Solve ( double *b, double *_pdDx )
{
  //copy b to _pdDx, so now _pdDx = b
  MathVecAssign( _iNumEq, _pdDx, b );
  Solve( _pdDx );
}

/* ============================ Solve2x2NonSym ============================ */
/*     Solves a system of linear equations: [a].{x} = {b}
**     where [a] is a square and nonsymmetric matrix.
**
**     Solve2x2NonSym( a, b )
**
**      a  -  Coefficient matrix  ( nonsymmetric )                (I)
**      b  -  Independent terms / system solution               (I/O)
**
**      Developed by: Professor Paulino?  
**      Modified  by:   
**
**		Added to code to handle CEE598 midterm problem - SL 050709
*/
/* ========================= Solve2x2NonSym ========================= */
void cCroutProfile :: Solve2x2NonSym( double *b, double *_pdDx )
{
   double x, y, det;
   
   det = _paA[0][0] * _paA[1][1] - _paA[0][1] * _paA[1][0];
   if ( fabs ( det ) < 1.0e-15 ) {
      b[0] = b[1] = 0.0;
   }
   else {
      x = ( b[0] * _paA[1][1] - b[1] * _paA[0][1] ) / det;
      y = ( b[1] * _paA[0][0] - b[0] * _paA[1][0] ) / det;
      _pdDx[0] = x;
      _pdDx[1] = y;
   }
}  /* End of Solve2x2NonSym */

