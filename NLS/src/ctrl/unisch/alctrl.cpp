// --------------------------------------------------------------------------
//
// ALCtrl.cpp - This file contains routines for  solving the incremental
//              and iteractive nonlinear problem by using the arc-length
//              control method.
//
// --------------------------------------------------------------------------

#include <math.h>
#include "ctrl/unisch/alctrl.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ========================== cArcLengthControl =============================

cArcLengthControl :: cArcLengthControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                     cUnifiedSchemes   ( pcModel, psControl, pcLinSys )
{
 _iDir = 1;

 _eCtrlType  = psControl->CtrlType;

 _pdDx1_i1   = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdDx1_i_11 = (double *) MemAlloc( _iNumEq, sizeof(double) );

 if( _eCtrlType == VARIABLE ) _pdSumDx = (double *) MemAlloc( _iNumEq, sizeof(double) );
}

// ========================== ~cArcLengthControl ============================

cArcLengthControl :: ~cArcLengthControl ( void )
{
 MemFree( _pdDx1_i1   );
 MemFree( _pdDx1_i_11 );

 if( _eCtrlType == VARIABLE ) MemFree( _pdSumDx );
}

// ============================== Lambda ====================================

void cArcLengthControl :: Lambda ( void )
{
	if( _iCurrIte == 1 ) {
		MathVecAssign( _iNumEq, _pdDx1_i_11, _pdDx1_i1 );

		_dLambda = sqrt( pow( _sControl.CtrlFactor, 2.0 ) /
			(1.0 + MathVecDot( _iNumEq, _pdDx1, _pdDx1 )) );
		
		MathVecAssign( _iNumEq, _pdDx1_i1, _pdDx1 );

		if( _iCurrStep != 1 ) 
			if( MathVecDot( _iNumEq, _pdDx1_i_11, _pdDx1_i1 ) < 0 ) 
				_iDir *= -1;
		
		_dLambda *= _iDir;
		
		if( _eCtrlType == VARIABLE ) {

			_dSumLambda = _dLambda;
			
			MathVecAssign1( _iNumEq, _pdSumDx, _dLambda, _pdDx1 );
		}
	} else {
		if( _eCtrlType == CONSTANT ) {
			_dLambda = - MathVecDot( _iNumEq, _pdDx1_i1, _pdDx2 ) /
			(MathVecDot( _iNumEq, _pdDx1_i1, _pdDx1 ) + 1);
		} else {
			_dLambda = - MathVecDot( _iNumEq, _pdSumDx, _pdDx2 ) /
				(MathVecDot( _iNumEq, _pdSumDx, _pdDx1 ) + _dSumLambda );
			
			_dSumLambda += _dLambda;
			
			MathVecAdd2( _iNumEq, _pdSumDx, _dLambda, _pdDx1, _pdDx2 );
		}
	}
}
/*
Almost agrees with Yang's Book...

_sControl.CtrlFactor is DeltaS

Neither matches the book exactly, the book calculates lambda_1 and 
DeltaU_1 and uses them, without updating, at every jth iteration

VARIABLE - updates lambda and DeltaU at every jth iteration
so lambda_1 is actually lambda_j-1 and DeltaU_1 is DeltaU_j-1
_pdDx1 is DeltaUhat
_pdDx2 is DeltaUbar
_pdDxSum is DeltaU

CONSTANT - sets lambda_1 = 1 and uses DeltaU_1 without updating
_pdDx1 is DeltaUhat
_pdDx2 is DeltaUbar
_pdDx1_i1 is DeltaU

***see hand written notes***

*/
