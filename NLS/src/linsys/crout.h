// -------------------------------------------------------------------------
//
// Crout.h - This header file contains the public data structure
//          definitions for the crout linear system class.
//
// -------------------------------------------------------------------------

#ifndef _CROUT_LINSYS_H_
#define _CROUT_LINSYS_H_

#include "linearsystem.h"

class cModel;

// -------------------------------------------------------------------------
// Crout Linear System class:

class cCroutProfile : public cLinSys
{
 protected:
  int      *_piProfile;
  double   **_paA;
  int      _isFactorized;

 public:
           cCroutProfile    ( void );
  virtual ~cCroutProfile    ( void );
  virtual  void Init  ( cModel * );
  virtual  void Zero  ( void );
  virtual  void AddA  (int, int, double );
  virtual  void FillA ( int, int, double);
  virtual  void Solve  ( double *x );
  virtual  void Solve  ( double *b, double *x );
  virtual  void Solve2x2NonSym( double *b, double *x );
};

#endif
