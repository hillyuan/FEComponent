// --------------------------------------------------------------------------
//
// StrnCtrl.cpp - This file contains routines for solving the incremental and
//                iteractive nonlinear problem by using  the  strain  control 
//                method.
//
// --------------------------------------------------------------------------

#include <math.h>
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "ctrl/unisch/strnctrl.h"
#include "linsys/linearsystem.h"
#include "model/model.h"

// --------------------------------------------------------------------------
// Public functions:

// ============================ cStrainControl ==============================

cStrainControl :: cStrainControl  ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                  cUnifiedSchemes ( pcModel, psControl, pcLinSys )
{
 _iCtrlEq     = psControl->CtrlEq;
 _dCtrlFactor = psControl->CtrlFactor;

 _pdDe1 = (double *) MemAlloc( _pcModel->NumEpsEq(), sizeof(double) );
 _pdDe2 = (double *) MemAlloc( _pcModel->NumEpsEq(), sizeof(double) );

 if( _sControl.CtrlType == VARIABLE )
 {
  _pdL      = (double *) MemAlloc( _pcModel->NumEpsEq(), sizeof(double) );
  _pdEpsEq  = (double *) MemAlloc( _pcModel->NumEpsEq(), sizeof(double) );
  _pdEpsPrv = (double *) MemAlloc( _pcModel->NumEpsEq(), sizeof(double) );

  MathVecZero( _pcModel->NumEpsEq(), _pdL      );
  MathVecZero( _pcModel->NumEpsEq(), _pdEpsPrv );
 }
}

// ========================== ~cStrainControl ===============================

cStrainControl :: ~cStrainControl ( void )
{
 MemFree( _pdDe1 );
 MemFree( _pdDe2 );

 if( _sControl.CtrlType == VARIABLE )
 {
  MemFree( _pdL      );
  MemFree( _pdEpsEq  );
  MemFree( _pdEpsPrv );
 }
}

// ============================== Lambda ====================================

void cStrainControl :: Lambda ( void )
{
 int i;

 _pcModel->DeltaStrainVector( _pdX, _pdDx1, _pdDe1 );

 if( _iCurrIte == 1 )
 {
  if( _iCurrStep != 1 && _sControl.CtrlType == VARIABLE )
  {
   _pcModel->StrainVector( _pdX, _pdEpsEq );

   int flag = 0;

   for( i = 0; i < _pcModel->NumEpsEq(); i++ )
   {
    _pdL[i] += fabs( _pdEpsEq[i] - _pdEpsPrv[i] );

    if( _pdEpsPrv[i] != 0.0 && _pdEpsEq[i] / _pdEpsPrv[i] > 1.0 ) flag = 1;

    _pdEpsPrv[i] = _pdEpsEq[i];
   }

   if( _iCurrStep > 2 && !flag )
   {
    printf( "\n\t ### Strain ratio < 1 !!! ###\n\n" );
   }

   int EqMax = 0;
   for( i = 1; i < _iNumEq; i++ )
    if( fabs( _pdDe1[i] ) > fabs( _pdDe1[EqMax] ) ) EqMax = i;

   if( EqMax != _iCtrlEq )
   {
    int Dir = (_pdDe1[EqMax] * _pdDe1[_iCtrlEq] > 0 ? 1 : -1 );

    _dCtrlFactor *= Dir * _pdL[EqMax] / _pdL[_iCtrlEq];
    _iCtrlEq      = EqMax;
   }
  }

  _dLambda = _dCtrlFactor / _pdDe1[_iCtrlEq];
 }
 else
 {
  _pcModel->DeltaStrainVector( _pdX, _pdDx2, _pdDe2 );
  _dLambda = - _pdDe2[_iCtrlEq] / _pdDe1[_iCtrlEq];
 }
}
