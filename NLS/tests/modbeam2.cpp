// ----------------------------------------------------------------------
//
// ModBeam2.cpp - Generic bidimensional beam model.
//
// ----------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include "modbeam2.h"
#include "utl.h"

/* Constructor for beam element object */
cModelBeam2D :: cModelBeam2D ( char *filename ) : cModel ( )
{
 int   id, i, j;
 int   num_load, num_sup;
 char  mod_type[80];
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

 //initialize arrays for coordinates, dof constraints and load
 _paCoord = (double**) MemAlloc( _iNumNodes, sizeof(double*) );
 _paDof   = (int**)    MemAlloc( _iNumNodes, sizeof(int*) );
 _paLoad  = (double**) MemAlloc( _iNumNodes, sizeof(double*) );
 
 for( i = 0; i < _iNumNodes; i++ )
 {
   //initalize 2nd dimension of dynamic arrays
   _paCoord[i] = (double*) MemAlloc( 2, sizeof(double) );
   _paDof[i]   = (int*)    MemAlloc( 3, sizeof(int) );
   _paLoad[i]  = (double*) MemAlloc( 3, sizeof(double) );

   //read in nodal coordinates
   fscanf( in, "%lf%lf", &_paCoord[i][0], &_paCoord[i][1] );
   for( j = 0; j < 3; j++ )
   {
     //initialze constraint to 0 and load to 0.0 for dof x, y, z (even though this is only 2d)
     _paDof[i][j]  = 0;
     _paLoad[i][j] = 0.0;
   }
 }

 //read in number of supports
 fscanf( in, "%d", &num_sup );

 for( i = 0; i < num_sup; i++ )
 {
  //read in constraints for dof x, y, z 
  fscanf( in, "%d", &id );
  for( j = 0; j < 3; j++ ) fscanf( in, "%d", &_paDof[id][j] );
 }

 //read in number of loads
 fscanf( in, "%d", &num_load );

 for( i = 0; i < num_load; i++ )
 {
  //read in loads for dof x, y, z
  fscanf( in, "%d", &id );
  for( j = 0; j < 3; j++ ) fscanf( in, "%lf", &_paLoad[id][j] );
 }

 //read in number of elements
 fscanf( in, "%d", &_iNumElms );

 //initalize dynamic array for element numbers and corresponsing nodes
 _paInc  = (int**)    MemAlloc( _iNumElms, sizeof(int*) );
 //initialize dynamic array for element properties
 _paProp = (double**) MemAlloc( _iNumElms, sizeof(double*) );

 for( i = 0; i < _iNumElms; i++ )
 {
   //initialize second dimesion of element arrays
   _paInc[i]  = (int*)    MemAlloc( 2, sizeof(int) );
   _paProp[i] = (double*) MemAlloc( 3, sizeof(double) );

   //read in beginning and ending nodes for each element
   for( j = 0; j < 2; j++ ) fscanf( in, "%d",  &_paInc[i][j] );
   //read in properties for each element: EA, GJ, EI
   for( j = 0; j < 3; j++ ) fscanf( in, "%lf", &_paProp[i][j] );
 }

 //read in number of displacement curves
 fscanf( in, "%d", &_iNumGra );

 //initialize dynamic array for displacement nodes and dofs
 _paGra = (int**) MemAlloc( _iNumGra, sizeof(int*) );

 for( i = 0; i < _iNumGra; i++ ) 
 {
   //initalize 2nd dimension of 2D array
   _paGra[i] = (int*) MemAlloc( 2, sizeof(int) );

   //read in node and degree of freedom
   for( j = 0; j < 2; j++ ) fscanf( in, "%d", &_paGra[i][j] );
 }

 fclose( in );


 _iNumEq = 0;

 for( i = 0; i < _iNumNodes; i++ )
  for( j = 0; j < 3; j++ )
   _paDof[i][j] = ( _paDof[i][j] == 0 ) ? _iNumEq++ : -1;


 _iNumEpsEq = _iNumElms;
}


/* Destructor for beam element object */
cModelBeam2D :: ~cModelBeam2D ( void )
{
  int i;

  for( i = 0; i < _iNumNodes; i++ )
  {
    MemFree( _paCoord[i] );
    MemFree( _paDof[i] );
    MemFree( _paLoad[i] );
  }
  MemFree( _paCoord );
  MemFree( _paDof );
  MemFree( _paLoad );
  
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

/* Initializes out file "node<number><dof (ie u, v, rz)>.out */
void cModelBeam2D :: Init ( void )
{
 char  fn[50];
 char *dir[3] = { "u", "v", "rz" };
 
 _afOut = (FILE **) MemAlloc( _iNumGra, sizeof(FILE *) );

 for( int i = 0; i < _iNumGra; i++ )
 {
  sprintf( fn, "node%d%s.out", _paGra[i][0], dir[_paGra[i][1]] );
  _afOut[i] = fopen( fn, "w" );
  fprintf( _afOut[i], "%f   %f\n", 0.0, 0.0 );
 }
}

/* Prints convergence */
void cModelBeam2D :: Convergence ( double dFactor, double *pdSol )
{
 for( int i = 0; i < _iNumGra; i++ )
 {
  double u = (_paDof[_paGra[i][0]][_paGra[i][1]] >= 0) ? pdSol[_paDof[_paGra[i][0]][_paGra[i][1]]] : 0.;
  fprintf( _afOut[i], "%f   %f\n", u, dFactor );
 }
}

/* Compute interval force vector for a beam element */
void cModelBeam2D :: InternalVector ( double *u, double *f )
{
 int    i, j, k;
 double a[6][6];
 double uelm[6], uloc[6];
 double felm[6], floc[6];

 MathVecZero( _iNumEq, f );

 for( i = 0; i < 6; i++ )
   for( j = 0; j < 6; j++ ) a[i][j] = 0.0;

 for( i = 0; i < _iNumElms; i++ )
 {
   double l = sqrt( pow( _paCoord[_paInc[i][1]][0] - _paCoord[_paInc[i][0]][0], 2.0 ) +
                  pow( _paCoord[_paInc[i][1]][1] - _paCoord[_paInc[i][0]][1], 2.0 ) );
   double c = (_paCoord[_paInc[i][1]][0] - _paCoord[_paInc[i][0]][0]) / l;
   double s = (_paCoord[_paInc[i][1]][1] - _paCoord[_paInc[i][0]][1]) / l;

   a[0][0] = a[1][1] = a[3][3] = a[4][4] = c;
   a[0][1] = a[3][4] =  s;
   a[1][0] = a[4][3] = -s;
   a[2][2] = a[5][5] = 1.0;

   for( j = 0; j < 2; j++ )
     for( k = 0; k < 3; k++ )
         uelm[3*j+k] = _paDof[_paInc[i][j]][k] >= 0 ? u[_paDof[_paInc[i][j]][k]] : 0.0;

   for( j = 0; j < 6; j++ )
   {
     uloc[j] = 0.0;
     for( k = 0; k < 6; k++ ) uloc[j] += a[j][k] * uelm[k];
   }

      floc[0] = -_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)+
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2);
      floc[1] = -_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)-
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2);
      floc[2] = l*_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc
[5]/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*(-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)
/2)+l*_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]
-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*(-(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc
[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/2)-1/l*_paProp[i]
[2]*(uloc[5]-uloc[2]);
      floc[3] = _paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)-
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2);
      floc[4] = _paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)+
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2);
      floc[5] = l*_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc
[5]/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*(-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)
/2)+l*_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]
-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*(-(uloc[4]-uloc[1])/l*sin(uloc[2]/2+uloc
[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/2)+1/l*_paProp[i]
[2]*(uloc[5]-uloc[2]);

      
   for( j = 0; j < 6; j++ )
   {
     felm[j] = 0.0;
     for( k = 0; k < 6; k++ ) felm[j] += a[k][j] * floc[k];
   }

    for( j = 0; j < 2; j++ )
     for( k = 0; k < 3; k++ )
         if( _paDof[_paInc[i][j]][k] >= 0 )
           f[_paDof[_paInc[i][j]][k]] += felm[3*j+k];
 }
}

/* Compute reference load vector */
void cModelBeam2D :: Reference ( double *pdReference )
{
 int i,j;
 for( i = 0; i < _iNumNodes; i++ )
  for( j = 0; j < 3; j++ )
   if( _paDof[i][j] >= 0 ) pdReference[_paDof[i][j]] = _paLoad[i][j];
}

/* Compute tangent matrix for beam element */
void cModelBeam2D :: TangentMatrix ( double *u, cLinSys *kt )
{
 int    i, j, k, m, n;
 double s1, uloc[6], uelm[6];
 double a[6][6];
 double kloc[6][6], kelm[6][6];

 kt->Zero();

 for( i = 0; i < 6; i++ )
   for( j = 0; j < 6; j++ ) a[i][j] = 0.0;

 for( i = 0; i < _iNumElms; i++ )
 {
   double l = sqrt( pow( _paCoord[_paInc[i][1]][0] - _paCoord[_paInc[i][0]][0], 2.0 ) +
                  pow( _paCoord[_paInc[i][1]][1] - _paCoord[_paInc[i][0]][1], 2.0 ) );
   double c = (_paCoord[_paInc[i][1]][0] - _paCoord[_paInc[i][0]][0]) / l;
   double s = (_paCoord[_paInc[i][1]][1] - _paCoord[_paInc[i][0]][1]) / l;

   a[0][0] = a[1][1] = a[3][3] = a[4][4] = c;
   a[0][1] = a[3][4] =  s;
   a[1][0] = a[4][3] = -s;
   a[2][2] = a[5][5] = 1.0;

   for( j = 0; j < 2; j++ )
     for( k = 0; k < 3; k++ )
         uelm[3*j+k] = _paDof[_paInc[i][j]][k] >= 0 ? u[_paDof[_paInc[i][j]][k]] : 0.0;

   for( j = 0; j < 6; j++ )
   {
     uloc[j] = 0.0;
     for( k = 0; k < 6; k++ ) uloc[j] += a[j][k] * uelm[k];
   }

     kloc[0][0] = _paProp[i][0]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0)+_paProp[i]
[1]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0);
     kloc[0][1] = kloc[1][0] = _paProp[i][0]/l*cos(uloc[2]/2+uloc[5]/2)*sin(uloc[2]/2+uloc
[5]/2)-_paProp[i][1]/l*sin(uloc[2]/2+uloc[5]/2)*cos(uloc[2]/2+uloc[5]/2);
     s1 = -_paProp[i][0]*cos(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)+
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[0][2] = kloc[2][0] = s1+_paProp[i][1]*sin(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)+_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[0][3] = kloc[3][0] = -_paProp[i][0]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0)-_paProp
[i][1]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0);
     kloc[0][4] = kloc[4][0] = -_paProp[i][0]/l*cos(uloc[2]/2+uloc[5]/2)*sin(uloc[2]/2+uloc
[5]/2)+_paProp[i][1]/l*sin(uloc[2]/2+uloc[5]/2)*cos(uloc[2]/2+uloc[5]/2);
     s1 = -_paProp[i][0]*cos(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)+
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[0][5] = kloc[5][0] = s1+_paProp[i][1]*sin(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)+_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[1][1] = _paProp[i][0]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0)+_paProp[i]
[1]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0);
     s1 = -_paProp[i][0]*sin(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)-
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[1][2] = kloc[2][1] = s1-_paProp[i][1]*cos(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)+_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[1][3] = kloc[3][1] = -_paProp[i][0]/l*cos(uloc[2]/2+uloc[5]/2)*sin(uloc[2]/2+uloc
[5]/2)+_paProp[i][1]/l*sin(uloc[2]/2+uloc[5]/2)*cos(uloc[2]/2+uloc[5]/2);
     kloc[1][4] = kloc[4][1] = -_paProp[i][0]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0)-_paProp
[i][1]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0);
     s1 = -_paProp[i][0]*sin(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)-
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[1][5] = kloc[5][1] = s1-_paProp[i][1]*cos(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)+_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2)/2;
     s1 = l*_paProp[i][0]*pow(-(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]
/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*_paProp[i][0]*((1.0
+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]
/2+uloc[5]/2)-1.0)*(-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/4-(uloc
[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)/4);
     kloc[2][2] = s1+l*_paProp[i][1]*pow(-(uloc[4]-uloc[1])/l*sin(uloc[2]/2+
uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*(-(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2
)/4+(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2)/4)+1/l*_paProp[i][2];
     s1 = _paProp[i][0]*cos(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)-
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[2][3] = kloc[3][2] = s1-_paProp[i][1]*sin(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)-_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2)/2;
     s1 = _paProp[i][0]*sin(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)+
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[2][4] = kloc[4][2] = s1+_paProp[i][1]*cos(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)-_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2)/2;
     s1 = l*_paProp[i][0]*pow(-(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]
/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*_paProp[i][0]*((1.0
+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]
/2+uloc[5]/2)-1.0)*(-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/4-(uloc
[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)/4);
     kloc[2][5] = kloc[5][2] = s1+l*_paProp[i][1]*pow(-(uloc[4]-uloc[1])/l*sin(uloc[2]/2+
uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*(-(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2
)/4+(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2)/4)-1/l*_paProp[i][2];
     s1 = _paProp[i][0]*cos(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)-
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[3][3] = _paProp[i][0]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0)+_paProp[i]
[1]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0);
     kloc[3][4] = kloc[4][3] = _paProp[i][0]/l*cos(uloc[2]/2+uloc[5]/2)*sin(uloc[2]/2+uloc
[5]/2)-_paProp[i][1]/l*sin(uloc[2]/2+uloc[5]/2)*cos(uloc[2]/2+uloc[5]/2);
     s1 = _paProp[i][0]*cos(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)-
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*sin(uloc[2]/2+uloc[5]/2)/2;
     kloc[3][5] = kloc[5][3] = s1-_paProp[i][1]*sin(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)-_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[4][4] = _paProp[i][0]/l*pow(sin(uloc[2]/2+uloc[5]/2),2.0)+_paProp[i]
[1]/l*pow(cos(uloc[2]/2+uloc[5]/2),2.0);
     s1 = _paProp[i][0]*sin(uloc[2]/2+uloc[5]/2)*(-(1.0+(uloc[3]-uloc[0])/l)*
sin(uloc[2]/2+uloc[5]/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2)+
_paProp[i][0]*((1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc
[1])/l*sin(uloc[2]/2+uloc[5]/2)-1.0)*cos(uloc[2]/2+uloc[5]/2)/2;
     kloc[4][5] = kloc[5][4] = s1+_paProp[i][1]*cos(uloc[2]/2+uloc[5]/2)*(-(uloc[4]-uloc[1]
)/l*sin(uloc[2]/2+uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]
/2)/2)-_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc
[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2))*sin(uloc[2]/2+uloc[5]/2)/2;
     s1 = l*_paProp[i][0]*pow(-(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]
/2)/2+(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*_paProp[i][0]*((1.0
+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)+(uloc[4]-uloc[1])/l*sin(uloc[2]
/2+uloc[5]/2)-1.0)*(-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/4-(uloc
[4]-uloc[1])/l*sin(uloc[2]/2+uloc[5]/2)/4);
     kloc[5][5] = s1+l*_paProp[i][1]*pow(-(uloc[4]-uloc[1])/l*sin(uloc[2]/2+
uloc[5]/2)/2-(1.0+(uloc[3]-uloc[0])/l)*cos(uloc[2]/2+uloc[5]/2)/2,2.0)+l*
_paProp[i][1]*((uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2)-(1.0+(uloc[3]-uloc
[0])/l)*sin(uloc[2]/2+uloc[5]/2))*(-(uloc[4]-uloc[1])/l*cos(uloc[2]/2+uloc[5]/2
)/4+(1.0+(uloc[3]-uloc[0])/l)*sin(uloc[2]/2+uloc[5]/2)/4)+1/l*_paProp[i][2];

   for( j = 0; j < 6; j++ )
       for( k = 0; k < 6; k++ )
       {
          kelm[j][k] = 0.0;
        for( m = 0; m < 6; m++ )
          for( n = 0; n < 6; n++ )
              kelm[j][k] += a[m][j] * kloc[m][n] * a[n][k];
       }

      for( j = 0; j < 2; j++ )
     for( k = 0; k < 3; k++ )
         if( _paDof[_paInc[i][j]][k] >= 0 )
          for( m = 0; m < 2; m++ )
            for( n = 0; n < 3; n++ )
              if( _paDof[_paInc[i][m]][n] >= 0 && _paDof[_paInc[i][m]][n] <= _paDof[_paInc[i][j]][k]  )
               kt->AddA(_paDof[_paInc[i][j]][k],_paDof[_paInc[i][m]][n], kelm[3*j+k][3*m+n]);
 }
}

int cModelBeam2D :: Profile ( int *piProfile )
{
#if 1
 int i,j;
 int **piSparsity = (int **) MemAlloc( _iNumEq , sizeof(int*) );
 for (i=0; i<_iNumEq; i++) {
   piSparsity[i] = (int *) MemAlloc( _iNumEq , sizeof(int) );
   for (j=0; j<_iNumEq; j++) {
     piSparsity[i][j] = 0;
   }
 }
 SparsityPattern ( piSparsity );
 for (i=0; i<_iNumEq; i++) {
   for (j=0; j<_iNumEq; j++) {
     if (piSparsity[i][j]==1) {
       piProfile[i]=j;
       break;
     }
   }
 }
 MemFree(piSparsity);
#else
 for( int i = 0; i < _iNumEq; i++ ) piProfile[i] = 0;
#endif
 return 1;
}

int cModelBeam2D :: SparsityPattern ( int **piSparsity )
{
 if (piSparsity==0) return( 0 );
 // Generate sparticity matrix
 int i,j; 
 for ( i = 0; i < _iNumElms; i++ )
 {
   for ( j = 0; j < 2; j++ )
   {
      for ( int k = 0; k < 3; k++ )
      {
         if ( _paDof[_paInc[i][j]][k] >= 0 )
         {
            for ( int m = 0; m < 2; m++ )
            {
               for ( int n = 0; n < 3; n++ )
               {
                  if ( ( _paDof[_paInc[i][m]][n] >= 0 ) && 
                       ( _paDof[_paInc[i][m]][n] <= _paDof[_paInc[i][j]][k] ) ) 
                  {
                     piSparsity[_paDof[_paInc[i][j]][k]][_paDof[_paInc[i][m]][n]] = 1;
                  }
               }
            }
         }
      }
   }
 }
 for ( i = 0; i < _iNumEq; i++ ) {
   for ( j = 0; j < i; j++ )
     piSparsity[j][i] = piSparsity[i][j];
 }
 return( 1 );
}

