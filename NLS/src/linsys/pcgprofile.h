// -------------------------------------------------------------------------
//
// PCGProfile.h - This header file contains the public data structure
//          definitions for the PCG linear system class.
//
// -------------------------------------------------------------------------

#ifndef _PCGPROFILE_LINSYS_H_
#define _PCGPROFILE_LINSYS_H_

#include "crout.h"

// -------------------------------------------------------------------------
// Crout Linear System class:

class cPCGProfile : public cCroutProfile
{
 protected:
  int      _iMaxIter;
  double   _dTol;

 public:
           cPCGProfile    ( int iMaxIter, double dTol );
  virtual ~cPCGProfile    ( void );
  virtual  void Solve  ( double *x );
  virtual  void Solve  ( double *b, double *x );
};

#endif
