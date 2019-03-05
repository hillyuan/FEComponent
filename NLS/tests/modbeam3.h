// ---------------------------------------------------------------------
//
// ModBeam3.h - Generic three-dimensional beam model.
//
// ---------------------------------------------------------------------

#ifndef _MODBEAM3_H_
#define _MODBEAM3_H_

#include <nls.h>

class cLinSys;

class cModelBeam3D : public cModel
{
 protected:
  int      _iNumNodes;  // Number of nodes
  int      _iNumElms;   // Number of elements
  int      **_paInc;    // List of incidences (initial and final node)
  double   **_paProp;   // List of properties (EA, GAx, GAy, EIx, EIy, and EIz)
  double   **_paCoord;  // List of coordinates (x, y, and z)
  int      **_paDof;    // Degrees of freedom (u, v, w, rx, ry, and rz)
  double   **_paLoad;   // List of loads (px, py, pz, mx, my, and mz)
  int      **_paGra;    // List of plots (node and dir)

 public:
       cModelBeam3D      ( char *filename );
      ~cModelBeam3D      ( void );

  void Init              ( void );
  void Convergence       ( double,   double  * );
  void InternalVector    ( double *, double  * );
  void Reference         ( double * );
  void TangentMatrix     ( double *, cLinSys * );
  int  Profile ( int * );
  int  SparsityPattern   ( int ** );
};

#endif
