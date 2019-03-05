// -------------------------------------------------------------------------
//
// WorkCtrl.h - This header file contains the public data structure
//              definitions for the work control method.
//
// -------------------------------------------------------------------------

#ifndef _WORKCTRL_H_
#define _WORKCTRL_H_

#include "unisch.h"

// -------------------------------------------------------------------------
// Work Control class:

class cWorkControl : public cUnifiedSchemes
{
 public:
        cWorkControl ( cModel *, sControl *, cLinSys * );
       ~cWorkControl ( void ) { }
  void  Lambda       ( void );
};

#endif
