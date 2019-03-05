// --------------------------------------------------------------------------
//
// ORCtrl.cpp - This file contains routines for  solving the incremental
//              and iteractive nonlinear problem by using the Orthogonal
//              Residual method.
//
// --------------------------------------------------------------------------

#include <stdio.h>

#include <math.h>
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "ctrl/unisch/orctrl.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ========================== cOrthResidualControl ==========================

cOrthResidualControl :: cOrthResidualControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                        cUnifiedSchemes      ( pcModel, psControl, pcLinSys )
{
 //what is this used for?
 _pdDx    = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdSumDx = (double *) MemAlloc( _iNumEq, sizeof(double) );
}

// ======================== ~cOrthResidualControl ===========================

cOrthResidualControl :: ~cOrthResidualControl ( void )
{
 MemFree( _pdDx    );
 MemFree( _pdSumDx );
}

// ============================== Lambda ====================================

void cOrthResidualControl :: Lambda ( void )
{
	if( _iCurrIte == 1 )
	{
		_dLambda = _sControl.CtrlFactor;

		if( MathVecDot( _iNumEq, _pdSumDx, _pdDx1 ) < 0 ) {
			_dLambda *= -1;
		}

		MathVecAssign1( _iNumEq, _pdSumDx, _dLambda, _pdDx1 );

		if( _iCurrStep == 1 )
			_dMaxDx = _sControl.CtrlIniFactor * MathVecNorm( _iNumEq, _pdSumDx );
		else
		{
			_dCurrDx = MathVecNorm( _iNumEq, _pdSumDx );

			if( _dCurrDx > _dMaxDx )
			{
				_dLambda *= _dMaxDx / _dCurrDx;

				MathVecScale( _iNumEq, _dMaxDx / _dCurrDx, _pdSumDx );
			}
		}
	} else {

		_dLambda = MathVecDot( _iNumEq, _pdUi, _pdSumDx ) /
             MathVecDot( _iNumEq, _pdReference, _pdSumDx ) - _dTotFactor;


		/* Change the sign of the load parameter between displacement limit points */
		//if ( _iCurrStep >= 13676 && _iCurrStep < 18345 ) _dLambda *= -1.0;

		MathVecAssign2( _iNumEq, _pdDx, _dLambda, _pdDx1, _pdDx2 );

		_dCurrDx = MathVecNorm( _iNumEq, _pdDx );

		if( _dCurrDx > _dMaxDx )
		{
			_dLambda *= _dMaxDx / _dCurrDx;

			MathVecScale( _iNumEq, _dMaxDx / _dCurrDx, _pdDx2 );
		}

		MathVecAdd2( _iNumEq, _pdSumDx, _dLambda, _pdDx1, _pdDx2 );
	}
}
