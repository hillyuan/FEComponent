// --------------------------------------------------------------------------
//
// SRCtrl.cpp - This file contains routines for  solving  the  incremental 
//              and iteractive nonlinear problem by using the strain ratio 
//              control method.
//
// --------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "ctrl/unisch/srctrl.h"
#include "model/model.h"
#include "linsys/linearsystem.h"

#define SRCBIG  1e300

// --------------------------------------------------------------------------
// Public functions:

// ======================== cStrainRatioControl =============================

cStrainRatioControl :: cStrainRatioControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                       cUnifiedSchemes     ( pcModel, psControl, pcLinSys )
{
 _iNumEpsEq = _pcModel->NumEpsEq( );

 _pdDe1 = (double *) MemAlloc( _iNumEpsEq, sizeof(double) );
 _pdDe2 = (double *) MemAlloc( _iNumEpsEq, sizeof(double) );

 _pdSumDx = (double *) MemAlloc( _iNumEq, sizeof(double) );

 _pdEpsCurr = (double *) MemAlloc( _iNumEpsEq, sizeof(double) );
 _pdEpsEq   = (double *) MemAlloc( _iNumEpsEq, sizeof(double) );

 _piEpsActv = (int *) MemAlloc( _iNumEpsEq, sizeof(int) );

 _dEpsZero = 1.0e-4;
}

// ======================== ~cStrainRatioControl ============================

cStrainRatioControl :: ~cStrainRatioControl ( void )
{
 MemFree( _pdDe1 );
 MemFree( _pdDe2 );

 MemFree( _pdSumDx );

 MemFree( _pdEpsCurr );
 MemFree( _pdEpsEq   );

 MemFree( _piEpsActv );
}

// ============================== Lambda ====================================

void cStrainRatioControl :: Lambda ( void )
{
 if( _iCurrStep == 1 )
 {
  if( _iCurrIte == 1 )
  {
   _dLambda = _sControl.CtrlIniFactor;
   MathVecAssign1( _iNumEq, _pdSumDx, _dLambda, _pdDx1 );
  }
  else
  {
   _dLambda = 0.0;
   MathVecAdd( _iNumEq, _pdSumDx, _pdDx2 );
  }
  return;
 }

 if( _iCurrIte == 1 )
 {
  int Dir = MathVecDot( _iNumEq, _pdSumDx, _pdDx1 ) < 0 ? -1 : 1;

  MathVecAssign( _iNumEpsEq, _pdEpsCurr, _pdEpsEq );

  _pcModel->StrainVector( _pdX, _pdEpsEq );

  _pcModel->DeltaStrainVector( _pdX, _pdDx1, _pdDe1 );

  _dLambda = Dir * SRCBIG;

  for( int i = 0; i < _iNumEpsEq; i++ )
   if( _pdDe1[i] != 0.0 && fabs( _pdEpsEq[i] ) > _dEpsZero )
   {
    double Curr = (_sControl.CtrlFactor - 1) * _pdEpsEq[i] / _pdDe1[i];
    if( Dir * Curr > 0 )
    {
     _piEpsActv[i] = 1;
     if( fabs( Curr ) < fabs( _dLambda ) ) _dLambda = Curr;
    }
    else _piEpsActv[i] = 0;
   }
   else _piEpsActv[i] = 0;

  if( _dLambda == (Dir * SRCBIG) )
  {
   _dLambda  = 0.0;
   _iCurrIte = _sControl.NumMaxIte;
   printf( "\n\n\t ### Consistent parameter not achieved !!! ###\n" );
  }

  MathVecAssign1( _iNumEq, _pdSumDx, _dLambda, _pdDx1 );
 }
 else
 {
  _pcModel->DeltaStrainVector( _pdX, _pdDx1, _pdDe1 );
  _pcModel->DeltaStrainVector( _pdX, _pdDx2, _pdDe2 );

  _pcModel->StrainVector( _pdX, _pdEpsCurr );

  int  IdLo = -1;
  int  IdHi = -1;
  double Lo = -SRCBIG;
  double Hi = +SRCBIG;

  for( int i = 0; i < _iNumEpsEq; i++ )
   if( _piEpsActv[i] && _pdDe1[i] != 0.0 && fabs( _pdEpsEq[i] ) > _dEpsZero )
   {
    double Curr = (_sControl.CtrlFactor * _pdEpsEq[i] - _pdEpsCurr[i] - _pdDe2[i]) / _pdDe1[i];
    double Prod = _pdDe1[i] * _pdEpsEq[i];
    if( Curr > Lo && Prod < 0 )
    {
     IdLo = i;
     Lo   = Curr;
    }
    else
     if( Curr < Hi && Prod > 0 )
     {
      IdHi = i;
      Hi   = Curr;
     }
   }

  if( Lo > Hi )
  {
   _dLambda  = 0.0;
   _iCurrIte = _sControl.NumMaxIte;
   printf( "\n\n\t ### Consistent parameter not achieved !!! ###" );
   printf( "\n\t EpsEqLo(%d) = %f\t Lo = %f\n\t EpsEqHi(%d) = %f\t Hi = %f\n\n",
           IdLo, _pdEpsEq[IdLo], Lo, IdHi, _pdEpsEq[IdHi], Hi );
  }
  else
   if( IdLo >= 0 && IdHi >= 0 )
    _dLambda = fabs( Hi * _pdDe1[IdLo] + _pdDe2[IdLo] ) >
               fabs( Lo * _pdDe1[IdHi] + _pdDe2[IdHi] ) ? Lo : Hi;
   else
    _dLambda = IdLo >= 0 ? Lo : Hi;

  MathVecAdd2( _iNumEq, _pdSumDx, _dLambda, _pdDx1, _pdDx2 );
 }
}

#undef SRCBIG

