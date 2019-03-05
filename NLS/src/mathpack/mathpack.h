// --------------------------------------------------------------------------
//
// MathPack.h - This header file contains the prototypes of the functions
//              used by "MathPack.cpp" file.
//
// --------------------------------------------------------------------------

#ifndef _MATHPACK_H_
#define _MATHPACK_H_

void   MathMatCrout   ( int, double **, double *, int, int * );
void   MathVecAdd     ( int, double *, double * );
void   MathVecAdd1    ( int, double *, double, double * );
void   MathVecAdd2    ( int, double *, double, double *, double * );
void   MathVecAssign  ( int, double *, double * );
void   MathVecAssign1 ( int, double *, double, double * );
void   MathVecAssign2 ( int, double *, double, double *, double * );
double MathVecDot     ( int, double *, double * );
double MathVecMaxNorm ( int, double * );
double MathVecNorm    ( int, double * );
void   MathVecScale   ( int, double, double * );
void   MathVecSub     ( int, double *, double * );
void   MathVecZero    ( int, double * );

#endif
