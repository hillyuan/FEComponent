// -------------------------------------------------------------------------
//
// MGDCtrl.h - This header file contains the public data  structure
//				definitions for the modified generalized displacement 
//				control method.
//
//				Prof. Eduardo Nobre
// -------------------------------------------------------------------------

#ifndef _MGDCTRL_H_
#define _MGDCTRL_H_

#include "unisch.h"
class cLinSys;

// -------------------------------------------------------------------------
// Generalized Displacement control class:

class cModGenDisplacementControl : public cUnifiedSchemes
{
 protected:
  double  _dN2_Dx1_11;
  double *_pdDx1_i1;
  double *_pdDx1_i_11;
  double *_pdSumDx;

 public:
        cModGenDisplacementControl ( cModel *, sControl *, cLinSys * );
       ~cModGenDisplacementControl ( void );
  void  Lambda                  ( void );
};

#endif
