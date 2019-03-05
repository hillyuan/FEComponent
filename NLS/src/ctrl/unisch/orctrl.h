// -------------------------------------------------------------------------
//
// ORCtrl.h - This header file contains the public data structure
//            definitions  for  the  Orthogonal  Residual  method.
//
// -------------------------------------------------------------------------

#ifndef _ORCTRL_H_
#define _ORCTRL_H_

#include "unisch.h"
class cLinSys;

// -------------------------------------------------------------------------
// Orthogonal Residual Control class:

class cOrthResidualControl : public cUnifiedSchemes
{
 protected:
  double  _dCurrDx;
  double  _dMaxDx;
  double *_pdDx;
  double *_pdSumDx;

 public:
        cOrthResidualControl ( cModel *, sControl *, cLinSys * );
       ~cOrthResidualControl ( void );
  void  Lambda               ( void );
};

#endif
