// ---------------------------------------------------------------------
//
// ModFunc.h - Descricao da classe modelo da funcao do paper.
//
// ---------------------------------------------------------------------

#ifndef _MODFUNC_H_
#define _MODFUNC_H_

#include <nls.h>

class cLinSys;

class cModelFunction : public cModel
{
 public:
       cModelFunction ( char *filename );
      ~cModelFunction ( void );

  void Init           ( void );
  void InternalVector ( double *, double  * );
  void Reference      ( double * );
  void TangentMatrix  ( double *, cLinSys * );
};

#endif

