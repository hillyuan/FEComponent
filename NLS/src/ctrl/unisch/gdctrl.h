// -------------------------------------------------------------------------
//
// GDCtrl.h - This header file contains the public data  structure
//            definitions for the generalized displacement control
//            method.
// -------------------------------------------------------------------------

#ifndef _GDCTRL_H_
#define _GDCTRL_H_

#include "unisch.h"
class cLinSys;

// -------------------------------------------------------------------------
// Generalized Displacement control class:

class cGenDisplacementControl : public cUnifiedSchemes
{
 protected:
  double  _dN2_Dx1_11;
  double *_pdDx1_i1;
  double *_pdDx1_i_11;
  double *_pdSumDx;

 public:
        cGenDisplacementControl ( cModel *, sControl *, cLinSys * );
       ~cGenDisplacementControl ( void );
  void  Lambda                  ( void );
};

#endif
