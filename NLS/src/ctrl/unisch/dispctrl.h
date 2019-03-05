// -------------------------------------------------------------------------
//
// DispCtrl.h - This header file contains the public data structure
//              definitions for the displacement control method.
//
// -------------------------------------------------------------------------

#ifndef _DISPCTRL_H_
#define _DISPCTRL_H_

#include "unisch.h"
class cLinSys;

// -------------------------------------------------------------------------
// Displacement Control class:

class cDisplacementControl : public cUnifiedSchemes
{
 protected:
  int     _iCtrlEq;
  double  _dCtrlFactor;
  double *_pdL;
  double *_pdXPrv;

 public:
        cDisplacementControl ( cModel *, sControl *, cLinSys * );
       ~cDisplacementControl ( void );
  void  Lambda               ( void );
};

#endif
