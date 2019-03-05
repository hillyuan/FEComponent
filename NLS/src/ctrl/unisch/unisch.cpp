// --------------------------------------------------------------------------
//
// UniSch.cpp - This file contains base class routines which are common among
//              a number of different control algorithm solution classes.
//
// --------------------------------------------------------------------------
//
// Format equation:
//
//   [a](i)(j-1) {dx}(i)(j) = lambda(i)(j) {r} + {u}(i)(j-1)
//
//   where:
//     [a]     -  Tangent matrix
//     {dx}    -  Incremental state vector
//     lambda  -  Incremental reference vector parameter
//     {r}     -  Reference vector
//     {u}     -  Unbalance vector
//     (i)(j)  -  (Step)(Iteraction)
//
// Reference:
//   "Solution Method for Nonlinear Problems with Multiple
//    Critical Points", Yeong-Bin Yang and Ming-Shan Shieh,
//    AIAA Journal, Vol 28, No 12, 1990.
//
// --------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "ctrl/unisch/unisch.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "linsys/linearsystem.h"
#include "model/model.h"

#define TIMING

// --------------------------------------------------------------------------
// Public functions:

// ======================== cUnifiedSchemes =================================

cUnifiedSchemes :: cUnifiedSchemes ( cModel *pcModel, sControl *psControl,
                   cLinSys *pcLinSys ) :
                   cControl        ( pcModel, psControl, pcLinSys )
{
 _pdDx1 = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdDx2 = (double *) MemAlloc( _iNumEq, sizeof(double) );

 _pdUe = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdUi = (double *) MemAlloc( _iNumEq, sizeof(double) );
 _pdR = (double *) MemAlloc( _iNumEq, sizeof(double) );

 _pdX = (double *) MemAlloc( _iNumEq, sizeof(double) );

 _dTotFactor = 0.0;

 _dNormReference = MathVecNorm( _iNumEq, _pdReference );
}

// ======================= ~cUnifiedSchemes =================================

cUnifiedSchemes :: ~cUnifiedSchemes ( void )
{
 MemFree( _pdDx1 );
 MemFree( _pdDx2 );

 MemFree( _pdUe );
 MemFree( _pdUi );
 MemFree( _pdR );

 MemFree( _pdX );

}

// ============================== Solver ====================================

void cUnifiedSchemes :: Solver ( void )
{
	//comments added by SL on 031010
	//step - outer loop _iCurrStep
	//iteration - inner loop _iCurrIte

	#ifdef TIMING
		double tic = clock();
	#endif

	int  ite;

	//initialize current step
	_iCurrStep = 1;
	
	//Outer Loop - Continues until max number of steps is reached or until diverges
	do {
		//initialize current iteration and set convergence to false	
		_iConvFlag = 0;
		_iCurrIte  = 1;

		//Inner loop - Continues until convergence or max permitted iterations is exceeded (ie divergence)
		do {
			//Compute tangent matrix if first iteration or if this is a standard update (ie compute K at every iteration
			//If this is modified then only computer tangent matrix on first iteration
			if(_iCurrIte == 1 || _sControl.UpdateType == STANDARD) {
				//get tangent matrix
				_pcModel->TangentMatrix(_pdX, _pcLinSys);
				//_pcLinSys solves Ax=b where A is _paA stored in _pcLinSys, and b is _pdDx1, _pdReference is copied to 
				//_pdDx1, the solver modifies b to be x, so result is _pdDx1 which is x
				//solving K*du_I = p, so _pdDx1 is du_I
				_pcLinSys->Solve(_pdReference, _pdDx1);
				//_pcLinSys->Solve2x2NonSym(_pdReference, _pdDx1);
			}

			if(_iCurrIte == 1) {
				//zero out _pdDx2
				MathVecZero(_iNumEq, _pdDx2);
			} else {
				//_pcLinSys solves Ax=b where A is _paA stored in _pcLinSys, and b is _pdDx2, _pdR is copied to 
				//_pdDx1, the solver modifies b to be x, so result is _pdDx2 which is x
				//solving K*du_II = r, so _pdDx2 is du_II
				_pcLinSys->Solve( _pdR, _pdDx2 );
				//_pcLinSys->Solve2x2NonSym( _pdR, _pdDx2 );
			}

			//compute iterative load factor
			Lambda( );

			//update total load factor for this step (outer loop)
			_dTotFactor += _dLambda;

			//compute _pdUe = _dLambda*_pdReference - external load vector for this iteration
			MathVecAdd1( _iNumEq, _pdUe, _dLambda, _pdReference   );

			//compute _pdX = _dLambda*_pdDx1 + _pdDx2, ie u = d_lambda*du_I + du_II
			MathVecAdd2( _iNumEq, _pdX,  _dLambda, _pdDx1, _pdDx2 );

			//sets up internal vector, gets stored in _pdUi
			_pcModel->InternalVector( _pdX, _pdUi );

			//calculate residual
		} while( !CheckConvergence( ) && ++_iCurrIte <= _sControl.NumMaxIte );

		//print results of convergering or not
		if( _iConvFlag ) {
			_pcModel->Convergence( _dTotFactor, _pdX );
			//printf( "\n" );
		}
		else {
			printf( "\n\n\t ### Convergence not achieved (Step %d) !!! ###\n\n", _iCurrStep );
			break;
		}
	//next step
	} while( ++_iCurrStep <= _sControl.NumMaxStep );

	#ifdef TIMING
	  double toc = clock();
	  printf("Solution time=%f\n", (toc-tic)/CLOCKS_PER_SEC);
	#endif

}


// --------------------------------------------------------------------------
// Protected function:

// ========================= CheckConvergence ===============================

int cUnifiedSchemes :: CheckConvergence ( void )
{
	//assign external load vector to residual, external load vector is all of the previous
	//external loads plus the lambda*referenceLoad for this iteration
	MathVecAssign( _iNumEq, _pdR, _pdUe );

	//right hand side of "governing" equation, ie K_(j-1)*deltaU_j = p_j - f_(j-1)
	//external load at j - internal load at j-1
	//_pdUe = _pdUe-_pdUi = p_j-1+lambda*p_bar-f_j-1
	MathVecSub( _iNumEq, _pdR, _pdUi );

	double Error = MathVecNorm( _iNumEq, _pdR );

	Error /= _dNormReference;

	_iConvFlag = Error <= _sControl.Tol;

	//printf( "   Step = %-5d   Ite = %-5d   Factor = %-11f   Error = %f\n",
	//		_iCurrStep, _iCurrIte, _dTotFactor, Error );

	return( _iConvFlag );
}

