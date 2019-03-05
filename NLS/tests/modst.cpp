// ----------------------------------------------------------------------
//
// ModST.cpp - Space truss.
//
// ----------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include "modst.h"
#include "utl.h"

cModelSpaceTruss :: cModelSpaceTruss ( char *filename ) : _iNumDofsOfNodes(3), cModel ( )
{
 int   i, j;
 int   id;
 int   num_load, num_sup;
 char mod_type[80];
 FILE *in = NULL;

 // Read data file

 in = fopen( filename, "r" );

 if( !in )
 {
  printf( "\n\n ### %s not defined !!! ###\n\n", filename );
  exit( -1 );
 }

 read_string( in, mod_type );

 fscanf( in, "%d", &_iNumNodes );

 _paCoord = (double**) MemAlloc( _iNumNodes, sizeof(double*) );
 _paLoad  = (double**) MemAlloc( _iNumNodes, sizeof(double*) );
 _paDof   = (int**)    MemAlloc( _iNumNodes, sizeof(int*) );

 for( i = 0; i < _iNumNodes; i++ )
 {
	 _paCoord[i] = (double*) MemAlloc( 3, sizeof(double) );
	 for( j = 0; j < 3; j++ ) fscanf( in, "%lf", &_paCoord[i][j] );

	 _paLoad[i] = (double*) MemAlloc( 3, sizeof(double) );
	 _paDof[i]  = (int*)    MemAlloc( 3, sizeof(int) );

	 for( j = 0; j < 3; j++ )
	 {
		 _paLoad[i][j] = 0.;
		 _paDof[i][j]  = 0;
	 }
 }

 fscanf( in, "%d", &num_sup );

 for( i = 0; i < num_sup; i++ )
 {
	 fscanf( in, "%d", &id );
	 for( j = 0; j < 3; j++ ) fscanf( in, "%d", &_paDof[id][j] );
 }

 fscanf( in, "%d", &num_load );

 for( i = 0; i < num_load; i++ )
 {
	 fscanf( in, "%d", &id );
	 for( j = 0; j < 3; j++ ) fscanf( in, "%lf", &_paLoad[id][j] );
 }

 fscanf( in, "%d", &_iNumElms );

 _paInc  = (int**)    MemAlloc( _iNumElms, sizeof(int*) );
 _paProp = (double**) MemAlloc( _iNumElms, sizeof(double*) );

 for( i = 0; i < _iNumElms; i++ ) 
 {
	 _paInc[i]  = (int*)    MemAlloc( 2, sizeof(int) );
	 _paProp[i] = (double*) MemAlloc( 3, sizeof(double) );

	 for( j = 0; j < 2; j++ ) fscanf( in, "%d",  &_paInc[i][j] );
	 for( j = 0; j < 2; j++ ) fscanf( in, "%lf", &_paProp[i][j] );
 }

 fscanf( in, "%d", &_iNumGra );

 _paGra = (int**) MemAlloc( _iNumGra, sizeof(int*) );

 for( i = 0; i < _iNumGra; i++ )
 {
	 _paGra[i] = (int*) MemAlloc( 2, sizeof(int) );
	 for( j = 0; j < 2; j++ ) fscanf( in, "%d", &_paGra[i][j] );
 }

 fclose( in );

 // Compute initial length of the bars

 for( i = 0; i < _iNumElms; i++ )
 {
	 _paProp[i][2] = 0.;
	 for( j = 0; j < 3; j++ ) _paProp[i][2] += pow( _paCoord[_paInc[i][1]][j] - _paCoord[_paInc[i][0]][j], 2.0 );
	 _paProp[i][2] = sqrt( _paProp[i][2] );
 }

 // Generate d.o.f.

 _iNumEq = 0;

 for( i = 0; i < _iNumNodes; i++ )
	 for( j = 0; j < 3; j++ )
	     _paDof[i][j] = ( _paDof[i][j] == 0 ) ? _iNumEq++ : -1;

 // Define number of strain components

 _iNumEpsEq = _iNumElms;
}

cModelSpaceTruss :: ~cModelSpaceTruss ( void )
{
	int i;

	for( i = 0; i < _iNumNodes; i++ )
	{
		MemFree( _paCoord[i] );
		MemFree( _paLoad[i] );
		MemFree( _paDof[i] );
	}
	MemFree( _paCoord );
	MemFree( _paLoad );
	MemFree( _paDof );

	for( i = 0; i < _iNumElms; i++ )
	{
		MemFree( _paInc[i] );
		MemFree( _paProp[i] );
	}
	MemFree( _paInc );
	MemFree( _paProp );

	for( i = 0; i < _iNumGra; i++ ) MemFree( _paGra[i] );
	MemFree( _paGra );
}

void cModelSpaceTruss :: Init ( void )
{
 char  fn[50];
 char *dir[3] = { "u", "v", "w" };

 _afOut = (FILE **) MemAlloc( _iNumGra, sizeof(FILE *) );
 for( int i = 0; i < _iNumGra; i++ )
 {
  sprintf( fn, "node%d%s.out", _paGra[i][0], dir[_paGra[i][1]] );
  _afOut[i] = fopen( fn, "w" );
  fprintf( _afOut[i], "%f   %f\n", 0.0, 0.0 );
 }
}

void cModelSpaceTruss :: Convergence ( double dFactor, double *pdSol )
{
 for( int i = 0; i < _iNumGra; i++ )
 {
  int    dof = _paDof[_paGra[i][0]][_paGra[i][1]];
  double u   = (dof >= 0) ? pdSol[dof] : 0.0;
  fprintf( _afOut[i], "%f   %f\n", u, dFactor );
 }
}

void cModelSpaceTruss :: StrainVector ( double *u, double *e )
{
 for( int i = 0; i < _iNumElms; i++ ) 
  e[i] = (Ldef( i, u ) - _paProp[i][2]) / _paProp[i][2];
}

void cModelSpaceTruss :: DeltaStrainVector ( double *u, double *du, double *de )
{
 int    i, j, k;
 int    edof[6];
 double cdef[2][3];
 double duelm[2][3];

 for( i = 0; i < _iNumElms; i++ )
 {
  double ldef = Ldef( i, u );

  for( j = 0; j < 2; j++ )
   for( k = 0; k < 3; k++ )
   {
    edof[3*j+k] = _paDof[_paInc[i][j]][k];
    cdef[j][k]  = _paCoord[_paInc[i][j]][k];
    duelm[j][k] = 0.0;
    if( edof[3*j+k] >= 0 ) 
    {
     cdef[j][k] +=  u[edof[3*j+k]];
     duelm[j][k] = du[edof[3*j+k]];
    }
   }

  de[i] = 0.0;

  for( j = 0; j < 3; j++ )
   de[i] += (cdef[1][j] - cdef[0][j]) * (duelm[1][j] - duelm[0][j]);

  de[i] *= 1.0 / (ldef * _paProp[i][2]);
 }
}

void cModelSpaceTruss :: InternalVector ( double *u, double *f )
{
 int    i, j, k;
 int    edof[6];
 double fe[6];
 double cdef[2][3];

 MathVecZero( _iNumEq, f );

 for( i = 0; i < _iNumElms; i++ )
 {
  double ldef = Ldef( i, u );
  double n    = _paProp[i][0] / _paProp[i][2] * (ldef - _paProp[i][2]) + _paProp[i][1];
  double coef = n / ldef;

  for( j = 0; j < 2; j++ )
   for( k = 0; k < 3; k++ )
   {
    edof[3*j+k] = _paDof[_paInc[i][j]][k];
    cdef[j][k]  = _paCoord[_paInc[i][j]][k];
    if( edof[3*j+k] >= 0 ) cdef[j][k] += u[edof[3*j+k]];
   }

  for( j = 0; j < 3; j++ )
  {
   fe[j]   = - coef * (cdef[1][j] - cdef[0][j]);
   fe[j+3] =   coef * (cdef[1][j] - cdef[0][j]);
  }

  for( j = 0; j < 6; j++ ) if( edof[j] >= 0 ) f[edof[j]] += fe[j];
 }
}

void cModelSpaceTruss :: Reference ( double *pdReference )
{
 int i,j;
 for( i = 0; i < _iNumNodes; i++ )
	 for( j = 0; j < _iNumDofsOfNodes; j++ )
         if( _paDof[i][j] >= 0 ) pdReference[_paDof[i][j]] = _paLoad[i][j];
}

void cModelSpaceTruss :: TangentMatrix ( double *u, cLinSys *kt )
{
	int    i, j, k;
	int    edof[6];
	double ke[6][6];
	double cdef[2][3];

	kt->Zero();

	for( i = 0; i < _iNumElms; i++ ) {
		//updated length of the element
		double ldef = Ldef( i, u );
		//stress in element (Hooke's Law)
		double n    = _paProp[i][0] / _paProp[i][2] * (ldef - _paProp[i][2]) + _paProp[i][1];

		//same code as Ldef
		for( j = 0; j < 2; j++ ) {
			for( k = 0; k < 3; k++ ) {
				edof[3*j+k] = _paDof[_paInc[i][j]][k];
				cdef[j][k]  = _paCoord[_paInc[i][j]][k];
				if( edof[3*j+k] >= 0 ) cdef[j][k] += u[edof[3*j+k]];
			}
		}

		//element tangent matrix
		for( j = 0; j < 3; j++ ) {
			for( k = 0; k < 3; k++ ) {
				double coef = (_paProp[i][0] - _paProp[i][1]) / (ldef * ldef * ldef) *
							(cdef[1][j] - cdef[0][j]) * (cdef[1][k] - cdef[0][k]);

				ke[j][k] = coef;
				if( j == k ) ke[j][k] += n / ldef;

				ke[j+3][k] = - coef;
				if( j == k ) ke[j+3][k] -= n / ldef;

				ke[j][k+3] = - coef;
				if( j == k ) ke[j][k+3] -= n / ldef;

				ke[j+3][k+3] = coef;
				if( j == k ) ke[j+3][k+3] += n / ldef;
			}
		}

		//what about coordinate transformation??? 
		//not needed for this problem, but is it available?

		//add any new contribution to total tangent matrix - only gives upper diagonal
		for( j = 0; j < 6; j++ ) {
			if( edof[j] >= 0 ) {
				for( k = 0; k < 6; k++ ) {
					if( edof[k] >= 0 && edof[j] >= edof[k] ) {
						kt->AddA(edof[j],edof[k],ke[j][k]);
					}
				}
			}
		}
	}
}

double cModelSpaceTruss :: Ldef ( int elm, double *u )
{
	int    i, j;
	int    edof[6];
	double cdef[2][3];

	//very clever!
	for( i = 0; i < 2; i++ ) {
		for( j = 0; j < 3; j++ ) {
			//getting element dofs for this element
			edof[3*i+j] = _paDof[_paInc[elm][i]][j];
			//getting node locations for this element 
			cdef[i][j]  = _paCoord[_paInc[elm][i]][j];
			//if its a free dof then edof[3*i+j] != 0, add the 
				//previous displacement associated with that dof
			if( edof[3*i+j] >= 0 ) cdef[i][j] += u[edof[3*i+j]];
		}
	}

	double ldef = 0.0;

	for( i = 0; i < 3; i++ ) ldef += pow( cdef[1][i] - cdef[0][i], 2.0 );

	//returns new length of the element - deformed length!!
	return( sqrt( ldef ) );
}