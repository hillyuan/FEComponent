// -------------------------------------------------------------------------
//
// StrnCtrl.h - This header file contains the public data structure
//              definitions for the strain control method.
//
// -------------------------------------------------------------------------

#ifndef _STRNCTRL_H_
#define _STRNCTRL_H_

#include "unisch.h"

class cLinSys;

// -------------------------------------------------------------------------
// Strain Control class:

class cStrainControl : public cUnifiedSchemes
{
 protected:
  int     _iCtrlEq;
  double  _dCtrlFactor;
  double *_pdDe1;
  double *_pdDe2;
  double *_pdL;
  double *_pdEpsEq;
  double *_pdEpsPrv;

 public:
        cStrainControl ( cModel *, sControl *, cLinSys * );
       ~cStrainControl ( void );
  void  Lambda         ( void );
};

#endif
