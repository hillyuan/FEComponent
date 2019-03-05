// -------------------------------------------------------------------------
//
// SRCtrl.h - This header file contains the public data structure
//            definitions for the strain ration control method.
//
// -------------------------------------------------------------------------

#ifndef _SRCTRL_H_
#define _SRCTRL_H_

#include "unisch.h"

class cLinSys;

// -------------------------------------------------------------------------
// Strain Ratio Control class:

class cStrainRatioControl : public cUnifiedSchemes
{
 protected:
  int     _iNumEpsEq;
  double *_pdDe1;
  double *_pdDe2;
  double *_pdEpsEq;
  double *_pdEpsCurr;
  double *_pdSumDx;
  int    *_piEpsActv;
  double  _dEpsZero;

 public:
        cStrainRatioControl ( cModel *, sControl *, cLinSys * );
       ~cStrainRatioControl ( void );
  void  Lambda              ( void );
};

#endif
