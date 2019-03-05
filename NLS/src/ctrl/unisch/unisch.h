// -------------------------------------------------------------------------
//
// UniSch.h - This header file contains the public data structure 
//            definitions for the control problem class.
//
// -------------------------------------------------------------------------

#ifndef _UNISCH_H_
#define _UNISCH_H_

#include "ctrl/ctrl.h"

class cLinSys;

// -------------------------------------------------------------------------
// Unified Schemes class:

class cUnifiedSchemes : public cControl
{
 protected:
  double  *_pdDx1;           // Incremental state vector due to reference vector
  double  *_pdDx2;           // Incremental state vector due to unbalance vector
  double   _dLambda;         // Incremental reference vector parameter
  double   _dNormReference;  // Euclidean norm of reference vector
  double  *_pdUe;            // External contribution for unbalance vector
  double  *_pdUi;            // Internal contribution for unbalance vector
  double  *_pdX;             // Total state variable vector - total displacement
  double  *_pdR;             // Total state variable vector

 public:
           cUnifiedSchemes ( cModel *, sControl *, cLinSys * );
  virtual ~cUnifiedSchemes ( void );
  void     Solver          ( void );

 protected:
          int  CheckConvergence ( void );
  virtual void Lambda           ( void ) = 0;
};

#endif
