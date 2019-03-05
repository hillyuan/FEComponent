// --------------------------------------------------------------------------
//
// MemPack.h - This header file contains the prototypes of the functions
//             used by "MemPack.cpp" file.
//
// --------------------------------------------------------------------------

#ifndef _MEMPACK_H_
#define _MEMPACK_H_

void*    MemAlloc        ( int, int );
void     MemFree         ( void * );
double** MemProfileAlloc ( int, int * );
void     MemProfileFree  ( double **, int, int * );

#endif
