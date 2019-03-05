// --------------------------------------------------------------------------
//
// DispCtrl.cpp - This file contains routines for solving the incremental and
//                iteractive  nonlinear  problem  by  using  the displacement
//                control method.
//
// --------------------------------------------------------------------------

#include <stdio.h>

#include <math.h>
#include "ctrl/unisch/dispctrl.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "linsys/linearsystem.h"


// --------------------------------------------------------------------------
// Public functions:

// ======================== cDisplacementControl ============================

cDisplacementControl :: cDisplacementControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                        cUnifiedSchemes      ( pcModel, psControl, pcLinSys )
{
 _iCtrlEq     = psControl->CtrlEq;
 _dCtrlFactor = psControl->CtrlFactor;

 if( _sControl.CtrlType == VARIABLE )
 {
  _pdL    = (double *) MemAlloc( _iNumEq, sizeof(double) );
  _pdXPrv = (double *) MemAlloc( _iNumEq, sizeof(double) );

  MathVecZero( _iNumEq, _pdL    );
  MathVecZero( _iNumEq, _pdXPrv );
 }
}

// ======================= ~cDisplacementControl ============================

cDisplacementControl :: ~cDisplacementControl ( void )
{
 if( _sControl.CtrlType == VARIABLE  )
 {
  MemFree( _pdL    );
  MemFree( _pdXPrv );
 }
}

// ============================== Lambda ====================================

void cDisplacementControl :: Lambda ( void )
{
	int i;

	if( _iCurrIte == 1 ) {
		if( _iCurrStep != 1 && _sControl.CtrlType == VARIABLE ) {
			for( i = 0; i < _iNumEq; i++ ) {
				_pdL[i]   += fabs( _pdX[i] - _pdXPrv[i] );
				_pdXPrv[i] = _pdX[i];
			}

			int EqMax = 0;

			//implementation of Fujii et al 1992 for best control parameter
			//find degree of freedom with the largest displacement
			for( i = 1; i < _iNumEq; i++ ) {
				if( fabs( _pdDx1[i] ) > fabs( _pdDx1[EqMax] ) ) EqMax = i;
			}
			//printf("Step: %d  ", _iCurrStep);
			//printf("DOF: 0, value %f  ", _pdDx1[0]);
			//printf("DOF: 1, value %f\n", _pdDx1[1]);
			
			//if the control DOF is not the one with the max displacement, switch the control dof
			//also keep the sign of the max dof
			if( EqMax != _iCtrlEq ) {
				int Dir = (_pdDx1[EqMax] * _pdDx1[_iCtrlEq] > 0 ? 1 : -1 ); 
				printf("Ctep: %d  ", _iCurrStep);
				printf("Control DOF: %d, value %f  ", _iCtrlEq, _pdDx1[_iCtrlEq]);
				printf("Max DOF: %d, value %f\n", EqMax, _pdDx1[EqMax]);
				_dCtrlFactor *= Dir * _pdL[EqMax] / _pdL[_iCtrlEq];
				_iCtrlEq      = EqMax;
			}
		}
		_dLambda = _dCtrlFactor / _pdDx1[_iCtrlEq];
	} else {
		_dLambda = - _pdDx2[_iCtrlEq] / _pdDx1[_iCtrlEq];
	}
}

/*
Agrees with Yang's Book:

_dCtrlFactor is DeltaU_qj
_pdDx1[_iCtrlEq] is DeltaUhat_qj
_pdDx2[_iCtrlEq] is DeltaUbar_qj

Not sure about current iteration part yet
*/
