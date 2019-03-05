// --------------------------------------------------------------------------
//
// WorkCtrl.cpp - This file contains routines for  solving  the  incremental
//                and iteractive nonlinear problem by using the work control 
//                method.
//
// --------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "workctrl.h"
#include "mathpack/mathpack.h"

// --------------------------------------------------------------------------
// Public functions:

// =========================== cWorkControl =================================

cWorkControl :: cWorkControl    ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                cUnifiedSchemes ( pcModel, psControl, pcLinSys )
{
}

// ============================== Lambda ====================================

void cWorkControl :: Lambda ( void )
{
	double Csp;
	double sign;

	if( _iCurrIte == 1 ) {
		Csp = _sControl.CtrlFactor / MathVecDot( _iNumEq, _pdReference, _pdDx1 );

		_dLambda = sqrt( fabs( Csp ) );

		if( Csp < 0.0 ){
			_dLambda *= -1;
		}
		
		if (0) { 
			printf("step: %d, ite: %d, ", _iCurrStep, _iCurrIte);
			printf("lambda: %f\n", _dLambda);
		}
	
	} else {
		_dLambda = - MathVecDot( _iNumEq, _pdReference, _pdDx2 ) /
               MathVecDot( _iNumEq, _pdReference, _pdDx1 );

		if (0) {
			printf("step: %d, ite: %d, ", _iCurrStep, _iCurrIte);
			printf("lambda: %f, ", _dLambda);
			printf("num: %f, ", MathVecDot( _iNumEq, _pdReference, _pdDx2 ));
			printf("denom: %f\n", MathVecDot( _iNumEq, _pdReference, _pdDx1 ));
		}
		
	}
}

//implementation matches N+1 theory, no add ons or modifcations