// -------------------------------------------------------------------------
//
// NRCtrl.h - This header file contains the public data structure
//            definitions for the standart Newton-Raphson method.
//
// -------------------------------------------------------------------------

#ifndef _NRCTRL_H_
#define _NRCTRL_H_

#include "unisch.h"

class cLinSys;

// -------------------------------------------------------------------------
// Newton-Raphson control class:

class cNewtonRaphson : public cUnifiedSchemes
{
 public: 
        cNewtonRaphson ( cModel *, sControl *, cLinSys * );
       ~cNewtonRaphson ( void ) { }
  void  Lambda         ( void );
};

#endif
