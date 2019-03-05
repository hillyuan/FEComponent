// ---------------------------------------------------------------------
//
// modfunc_mt.h - Header file for Function_MT class
//
// ---------------------------------------------------------------------

#include <nls.h>

class cLinSys;

class cModelFunction_MT : public cModel
{
 public:
       cModelFunction_MT ( char *filename );
      ~cModelFunction_MT ( void );

  void Init           ( void );
  void InternalVector ( double *, double  * );
  void Reference      ( double * );
  void TangentMatrix  ( double *, cLinSys * );
};

