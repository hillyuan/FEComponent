// ---------------------------------------------------------------------
//
// Model.h - Descricao da classe de modelos.
//
// ---------------------------------------------------------------------

#ifndef _MODEL_H_
#define _MODEL_H_

#include <stdio.h>

class cLinSys;

class cModel
{
 protected:
  int      _iNumEpsEq;
  int      _iNumEq;
  int      _iNumGra;
  FILE   **_afOut;

 public:
           cModel   ( void ) { }; //already implemented to do nothing
  virtual ~cModel   ( void ); //required, implemented in model.cpp

  int     NumEpsEq  ( void ); //required, implemented in model.cpp
  int     NumEq     ( void ); //required, implemented in model.cpp
  virtual int Profile   ( int * ); //required, implemented in model.cpp
  virtual int SparsityPattern  ( int ** ); //required, implemented in model.cpp

  virtual void Convergence ( double, double * ); //required, implemented in model.cpp

  virtual void Init              ( void )                         = 0; //optional, implemented in children of model.cpp
  virtual void InternalVector    ( double *, double  * )          = 0; //optional, implemented in children of model.cpp
  virtual void Reference         ( double * )                     = 0; //optional, implemented in children of model.cpp
  virtual void TangentMatrix     ( double *, cLinSys * )          = 0; //optional, implemented in children of model.cpp
  virtual void StrainVector      ( double *, double  * )           { } //already implemented to do nothing
  virtual void DeltaStrainVector ( double *, double  *, double * ) { } //already implemented to do nothing

 protected:
  void InitFile ( void ); //required, implemented in model.cpp
};

#endif

