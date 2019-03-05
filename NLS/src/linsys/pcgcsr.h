// -------------------------------------------------------------------------
//
// PCGCSR.h - This header file contains the public data structure
//          definitions for the PCG linear system class.
//
// -------------------------------------------------------------------------

#ifndef _PCGCSR_LINSYS_H_
#define _PCGCSR_LINSYS_H_

#include "linearsystem.h"

// -------------------------------------------------------------------------
// Crout Linear System class:

class cPCGCSR : public cLinSys
{
 protected:
  int      _iMaxIter;
  double   _dTol;
  int     *_iRowPtr;
  int     *_iColInd;
  double  *_dVal;

 public:
           cPCGCSR    ( int iMaxIter, double dTol );
  virtual ~cPCGCSR    ( void );
  virtual  void Init  ( cModel * );
  virtual  void Zero  ( void );
  virtual  void AddA  (int, int, double );
  virtual  void Solve  ( double *b, double *x );
};

class cPCGCSRSym : public cPCGCSR
{
 protected:
 public:
           cPCGCSRSym    ( int iMaxIter, double dTol ): cPCGCSR(iMaxIter,dTol) {}
  virtual ~cPCGCSRSym    ( void ) {}
  virtual  void AddA  (int, int, double );
};

#endif
