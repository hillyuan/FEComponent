// --------------------------------------------------------------------------
//
// NRCtrl.cpp - This file contains routines for solving the incremental and
//              iteractive nonlinear problem by  using  the  Newton-Raphson
//              method.
//
// --------------------------------------------------------------------------

#include "ctrl/unisch/nrctrl.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ========================== cNewtonRaphson ================================

cNewtonRaphson :: cNewtonRaphson  ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                  cUnifiedSchemes ( pcModel, psControl, pcLinSys )
{
}

// ============================== Lambda ====================================

void cNewtonRaphson :: Lambda ( void )
{
 _dLambda = _iCurrIte == 1 ? _sControl.CtrlFactor : 0.0;
}
