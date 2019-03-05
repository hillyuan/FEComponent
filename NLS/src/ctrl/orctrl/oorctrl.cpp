// --------------------------------------------------------------------------
//
// OORCtrl.cpp - This file contains routines for solving the  incremental 
//               and iteractive nonlinear problem by using the Orthogonal
//               Residual method.
//
// --------------------------------------------------------------------------

#include <stdio.h>
#include "oorctrl.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "model/model.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ======================= cOldOrthResidualControl ==========================

cOldOrthResidualControl :: cOldOrthResidualControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                           cControl                ( pcModel, psControl, pcLinSys )
{
}

// ================================ Solver ==================================

void cOldOrthResidualControl :: Solver ( void )
{
 int i;
 double *ri;
 double *fi;
 double *ft;
 double *dx1;
 double *dxi;
 double *xt;
 double *dpn;
 double eps, nref, totfactor;
 double ducurr, dumax=0, nri;

 dx1 = (double *) MemAlloc( _iNumEq, sizeof(double) );
 dxi = (double *) MemAlloc( _iNumEq, sizeof(double) );
 dpn = (double *) MemAlloc( _iNumEq, sizeof(double) );
 ri  = (double *) MemAlloc( _iNumEq, sizeof(double) );
 fi  = (double *) MemAlloc( _iNumEq, sizeof(double) );
 ft  = (double *) MemAlloc( _iNumEq, sizeof(double) );
 xt  = (double *) MemAlloc( _iNumEq, sizeof(double) );

 MathVecAssign1( _iNumEq, dpn, _sControl.CtrlFactor, _pdReference );

 nref = MathVecNorm( _iNumEq, _pdReference );

 MathVecZero( _iNumEq, xt );

 _iCurrStep = 1;

 do
 {
  //MathVecAssign( _iNumEq, dx1, dpn );

  _pcModel->TangentMatrix( xt, _pcLinSys );

  _pcLinSys->Solve( dpn, dx1 );

  if( _iCurrStep == 1 ) 
   dumax = _sControl.CtrlIniFactor * MathVecMaxNorm( _iNumEq, dx1 ); 

  ducurr = MathVecMaxNorm( _iNumEq, dx1 );

  if( ducurr > dumax ) MathVecScale( _iNumEq, dumax/ ducurr, dx1 );

  if( MathVecDot( _iNumEq, dxi, dx1 ) < 0 )
  {
   MathVecScale( _iNumEq, -1.0, dx1 );
   MathVecScale( _iNumEq, -1.0, dpn );
  }

  MathVecAssign( _iNumEq, dxi, dx1 );

  _iCurrIte = 1;

  do
  {
   _iConvFlag = 0;

   for( i = 0; i < _iNumEq; i++ ) ri[i] = xt[i] + dxi[i];

   if( _sControl.UpdateType == STANDARD )
   {
    _pcModel->TangentMatrix( ri, _pcLinSys );
    _pcLinSys->Solve(0);
   }

   _pcModel->InternalVector( ri, fi );

   for( i = 0; i < _iNumEq; i++ ) ri[i] = ft[i] - fi[i];

   eps = - MathVecDot( _iNumEq, ri, dxi ) / MathVecDot( _iNumEq, dpn, dxi );

   MathVecAdd1( _iNumEq, ri, eps, dpn );

   nri = MathVecNorm( _iNumEq, ri );

   _pcLinSys->Solve(ri);

   ducurr = MathVecMaxNorm( _iNumEq, ri );

   if( ducurr > dumax ) MathVecScale( _iNumEq, dumax / ducurr, ri );

   MathVecAdd( _iNumEq, dxi, ri );

   if( nri < (_sControl.Tol * nref) )
   {
    _iConvFlag = 1;
    break;
   }

  } while( ++_iCurrIte <= _sControl.NumMaxIte );

  if( !_iConvFlag )
  {
   printf( "\n\n\t ### Convergence not achieved (Step %d) !!! ###\n\n", _iCurrStep );
   break;
  }

  MathVecAdd ( _iNumEq, xt, dxi );
  MathVecAdd1( _iNumEq, ft, eps, dpn );

  totfactor = MathVecDot( _iNumEq, ft, _pdReference ) / (nref*nref);

  // Convergence feedback

  _pcModel->Convergence( totfactor, xt );

 } while( ++_iCurrStep <= _sControl.NumMaxStep );

 // Free memory

 MemFree( dpn );
 MemFree( dx1 );
 MemFree( dxi );
 MemFree( fi  );
 MemFree( ft  );
 MemFree( ri  );
 MemFree( xt  );

}
