// -------------------------------------------------------------------------
//
// LinearSystem.h - This header file contains the public data structure
//          definitions for the linear system class.
//
// -------------------------------------------------------------------------

#ifndef _LINEAR_SYSTEM_H_
#define _LINEAR_SYSTEM_H_

class cModel;

// -------------------------------------------------------------------------
// Linear System class:

class cLinSys
{
 protected:
  int       _iNumEq;

 public:
           cLinSys    ( void );
  virtual ~cLinSys    ( void );
  virtual  void Init  ( cModel * ) = 0;
  virtual  void Zero  ( void ) = 0;
  virtual  void AddA  (int, int, double ) = 0;
  virtual  void FillA ( int, int, double) = 0;
  virtual  void Solve  ( double *b );
  virtual  void Solve  ( double *b, double *x ) = 0;
  virtual  void Solve2x2NonSym( double *b, double *_pdDx ) = 0;
};

#endif
