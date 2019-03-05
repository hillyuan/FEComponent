// ---------------------------------------------------------------------
//
// ModST.h - Space Truss.
//
// ---------------------------------------------------------------------

#ifndef _MODST_H_
#define _MODST_H_

#include <nls.h>

class cLinSys;

class cModelSpaceTruss : public cModel
{
 protected:
  int      _iNumNodes;	// Number of nodes
  int      _iNumDofsOfNodes;
  int      _iNumElms;	// Number of elements
  double **_paCoord;	// List of coordinates (x, y and z)
  double **_paLoad;		// List of loads (px, py and pz)
  int    **_paDof;		// List of degree of freedom (u, v and w)
  int    **_paInc;		// List of incidences (initial and final)
  double **_paProp;		// List of properties (EA, N0 and L0)
  int    **_paGra;		// List of plots (node and direction)

 public:
       cModelSpaceTruss  ( char *filename );
      ~cModelSpaceTruss  ( void );

  void Init              ( void );
  void Convergence       ( double,   double  * );
  void InternalVector    ( double *, double  * );
  void Reference         ( double * );
  void TangentMatrix     ( double *, cLinSys * );
  void DeltaStrainVector ( double *, double  *, double * );
  void StrainVector      ( double *, double  * );

 protected:
  double Ldef ( int, double * );
};

#endif
