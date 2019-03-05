// ---------------------------------------------------------------------
//
// ModBeam2.h - Generic bidimensional beam model.
//
// ---------------------------------------------------------------------

#ifndef _MODBEAM2_H_
#define _MODBEAM2_H_

#include <nls.h>

class cLinSys;

class cModelBeam2D : public cModel
{
 protected:
  int          _iNumNodes;	// Number of nodes
  int          _iNumElms;	// Number of elements
  int        **_paInc;		// List of incidences (initial and final node)
  double     **_paProp;		// List of properties (EA, GA and EI)
  double     **_paCoord;	// List of coordinates (x and y)
  int        **_paDof;		// Degree of freedom (u, v, rz)
  double     **_paLoad;		// List of loads (px, py and mz)
  int        **_paGra;		// List of plots (node and dir)

 public:
       cModelBeam2D      ( char *filename );
      ~cModelBeam2D      ( void );

  void Init              ( void );
  void Convergence       ( double,   double  * );
  void InternalVector    ( double *, double  * );
  void Reference         ( double * );
  void TangentMatrix     ( double *, cLinSys * );
  int  Profile ( int * );
  int  SparsityPattern   ( int ** );
};

#endif
