// -------------------------------------------------------------------------
//
// ARCtrl.h - This header file contains the public data structure
//            definitions for the arc-length control method.
//
// -------------------------------------------------------------------------

#ifndef _ALCTRL_H_
#define _ALCTRL_H_

#include "unisch.h"

class cLinSys;

// -------------------------------------------------------------------------
// Arc-Length Control class:

class cArcLengthControl : public cUnifiedSchemes
{
 protected:
  eCtrlType  _eCtrlType;
  int        _iDir;
  double    *_pdDx1_i1;
  double    *_pdDx1_i_11;
  double    *_pdSumDx;
  double     _dSumLambda;

 public:
        cArcLengthControl ( cModel *, sControl *, cLinSys * );
       ~cArcLengthControl ( void );
  void  Lambda            ( void );
};

#endif
