// --------------------------------------------------------------------------
//
// PCGCSR.cpp - PCG Linear System class.
//
// --------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "linsys/pcgcsr.h"
#include "mempack/mempack.h"
#include "mathpack/mathpack.h"
#include "model/model.h"
#include <linsys.h>

// --------------------------------------------------------------------------
// Public functions:

// ============================= cLinSys ===================================

cPCGCSR :: cPCGCSR (int iMaxIter, double dTol ) : _iMaxIter(iMaxIter), _dTol(dTol)
{
}

// ============================= ~cLinSys ==================================

cPCGCSR :: ~cPCGCSR ( void )
{
 MemFree(_iRowPtr);
 MemFree(_iColInd);
 MemFree(_dVal);
}

void cPCGCSR :: Init ( cModel *pcModel )
{
 int i,j,p;
 _iNumEq      = pcModel->NumEq( );
 int **piSparsity = (int **) MemAlloc( _iNumEq , sizeof(int*) );
 for (i=0; i<_iNumEq; i++) {
   piSparsity[i] = (int *) MemAlloc( _iNumEq , sizeof(int) );
   for (j=0; j<_iNumEq; j++) {
     piSparsity[i][j] = 0;
   }
 }
 pcModel->SparsityPattern( piSparsity );

 _iRowPtr = (int *) MemAlloc( _iNumEq + 1 , sizeof(int) );
 _iRowPtr[0] =0;
 for (i=0; i<_iNumEq; i++) {
   for (j=0,p=0; j<_iNumEq; j++) {
     if (piSparsity[i][j]) p++;
   }
   _iRowPtr[i+1] =_iRowPtr[i]+p;
 }
 _iColInd = (int *) MemAlloc( _iRowPtr[_iNumEq] , sizeof(int) );
 _dVal = (double *) MemAlloc( _iRowPtr[_iNumEq] , sizeof(double) );
 for (i=0; i<_iNumEq; i++) {
   for (j=0,p=0; j<_iNumEq; j++) {
     if (piSparsity[i][j]) {
       _iColInd[_iRowPtr[i]+p] = j;
       p++;
     }
   }
 }
#if 0
 printf("cPCGCSR::Init | n = %d\n", _iNumEq);
 printf("cPCGCSR::Init | nonzeros = %d\n", _iRowPtr[_iNumEq]);
 for (i=0; i<_iNumEq; i++)
   printf("cPCGCSR::Init | row_length(%d) = %d\n", i, _iRowPtr[i+1]-_iRowPtr[i]);
 exit(-1);
#endif
 MemFree(piSparsity);
}

void cPCGCSR :: Zero ( )
{
 int i,j;
 for(i=0;i<_iNumEq;i++) {
    for(j=_iRowPtr[i];j<_iRowPtr[i+1];j++) {
       _dVal[j] = 0.;
    }
 }
}

void cPCGCSR :: AddA ( int i, int j, double k )
{
 int p;
 // search
 for (p=_iRowPtr[i]; p<_iRowPtr[i+1]; p++)
   if (_iColInd[p] == j) break;

 if (p == _iRowPtr[i+1]) {
   printf("cPCGCSR::AddA element(%d,%d) is not defined.\n", i,j);
   return;
 }
 _dVal[p]+=k;
}

void cPCGCSR :: Solve ( double *b, double *x )
{
  int iter;
  PCG_CSR(_iRowPtr, _iColInd, _dVal,
              _iNumEq, b, x, &iter, _iMaxIter, _dTol);
  //printf("_PCG_CSR - iterations = %d \n", iter);
}


void cPCGCSRSym :: AddA ( int i, int j, double k )
{
 int p;
 if (_iRowPtr[i] == _iRowPtr[i+1]) {
   printf("cPCGCSRSym::AddA(%d,%d) : empty line.\n", i,j);
   return;
 }

 // search
 for (p=_iRowPtr[i]; p<_iRowPtr[i+1]; p++)
   if (_iColInd[p] == j) break;

 if (p == _iRowPtr[i+1]) {
   printf("cPCGCSRSym::AddA(%d,%d) : search doesn't find the column.\n", i,j);
   return;
 }
 _dVal[p]+=k;

 // skip if it's a main diagonal entry
 if ( i == j ) return;

 // symmetric part
 for (p=_iRowPtr[j]; p<_iRowPtr[j+1]; p++)
  if ( (_iColInd[p] == i) && (_iColInd[p] != j) ) break;

 if (p == _iRowPtr[i+1]) {
   printf("cPCGCSRSym::AddA element(%d,%d) is not defined.\n", i,j);
   return;
 }
 _dVal[p]+=k;
}
