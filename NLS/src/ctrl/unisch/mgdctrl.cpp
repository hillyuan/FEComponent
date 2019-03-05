// --------------------------------------------------------------------------
//
// MGDCtrl.cpp - This file contains routines for solving  the  incremental
//              and iteractive nonlinear problem by using the modified
//				generalized displacement control method.
//
//				Prof. Eduardo Nobre
// --------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include "ctrl/unisch/mgdctrl.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ======================= cModGenDisplacementControl ==========================

cModGenDisplacementControl :: cModGenDisplacementControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys ) :
                           cUnifiedSchemes         ( pcModel, psControl, pcLinSys )
{
 _pdDx1_i1   = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdDx1_i_11 = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdSumDx    = (double *) MemAlloc( _iNumEq, sizeof(double) );
}

// ======================= ~cModGenDisplacementControl =========================

cModGenDisplacementControl :: ~cModGenDisplacementControl ( void )
{
 MemFree( _pdDx1_i1   );
 MemFree( _pdDx1_i_11 );
 MemFree( _pdSumDx    );
}

// ============================== Lambda ====================================

void cModGenDisplacementControl :: Lambda ( void )
{
	if( _iCurrIte == 1 ) {
		if( _iCurrStep == 1 ) {
			printf("here\n");
			MathVecZero( _iNumEq, _pdSumDx );

			MathVecAssign( _iNumEq, _pdDx1_i1,   _pdDx1 );
			MathVecAssign( _iNumEq, _pdDx1_i_11, _pdDx1 );

			_dN2_Dx1_11 = MathVecDot( _iNumEq, _pdDx1, _pdDx1 );

			_dLambda = _sControl.CtrlFactor;
		} else {
			MathVecAssign( _iNumEq, _pdDx1_i_11, _pdDx1_i1 );
			MathVecAssign( _iNumEq, _pdDx1_i1, _pdDx1 );

			double Gsp = _dN2_Dx1_11 / MathVecDot( _iNumEq, _pdDx1, _pdDx1 );
			//printf("step = %d ;  sign of gsp = %d\n", _iCurrStep, Gsp < 0 ? -1 : 1);

			_dLambda = sqrt( Gsp ) * _sControl.CtrlFactor;

			if( MathVecDot( _iNumEq, _pdSumDx, _pdDx1 ) < 0 ) _dLambda *= -1;

			MathVecZero( _iNumEq, _pdSumDx );
		}
	} else {
		_dLambda = - MathVecDot( _iNumEq, _pdDx1_i1, _pdDx2 ) /
					MathVecDot( _iNumEq, _pdDx1_i1, _pdDx1 );
	}

	MathVecAdd2( _iNumEq, _pdSumDx, _dLambda, _pdDx1, _pdDx2 );
}
