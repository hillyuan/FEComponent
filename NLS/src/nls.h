/**
  NLS
  
  NLS is a software developed by Tecgraf/PUC-Rio & UIUC.
  It is requested that the NLS users provide the appropriate
  credits and references to the work.
  
  References:
  1 -  Leon S, Paulino GH, Pereira A, Menezes IFM and Lages EN (2012)
       A Unified Library of Nonlinear Solution Schemes,
       Appl. Mech. Rev. 64(4), 040803 DOI:10.1115/1.4006992

**/

#ifndef _NLS_H
#define _NLS_H

#include "ctrl/ctrl.h"
#include "ctrl/orctrl/oorctrl.h"
#include "ctrl/unisch/alctrl.h"
#include "ctrl/unisch/dispctrl.h"
#include "ctrl/unisch/gdctrl.h"
#include "ctrl/unisch/mgdctrl.h"
#include "ctrl/unisch/nrctrl.h"
#include "ctrl/unisch/orctrl.h"
#include "ctrl/unisch/srctrl.h"
#include "ctrl/unisch/strnctrl.h"
#include "ctrl/unisch/unisch.h"
#include "ctrl/unisch/workctrl.h"
#include "linsys/crout.h"
#include "linsys/linearsystem.h"
#include "linsys/pcgcsr.h"
#include "linsys/pcgprofile.h"
#include "mathpack/mathpack.h"
#include "mempack/mempack.h"
#include "model/model.h"

#endif
