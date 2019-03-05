// ---------------------------------------------------------------------
//
// ModFunc_nobre.h 
//
// ---------------------------------------------------------------------

#include <nls.h>

class cLinSys;

class cModelFunction_Nobre : public cModel
{
 public:
       cModelFunction_Nobre ( char *filename );
      ~cModelFunction_Nobre ( void );

  void Init           ( void );
  void InternalVector ( double *, double  * );
  void Reference      ( double * );
  void TangentMatrix  ( double *, cLinSys * );
};

