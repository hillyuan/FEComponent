// -------------------------------------------------------------------------
//
// OORCtrl.h - This header file contains the public  data  structure
//             definitions for the orthogonal residual control class.
//
// -------------------------------------------------------------------------

#ifndef _OORCTRL_H_
#define _OORCTRL_H_

#include "ctrl/ctrl.h"

class cLinSys;

// -------------------------------------------------------------------------
// Orthogonal Residual class:

class cOldOrthResidualControl : public cControl
{
 public:
           cOldOrthResidualControl ( cModel *, sControl *, cLinSys * );
  virtual ~cOldOrthResidualControl ( void ) { }
  void     Solver                  ( void );
};

#endif
