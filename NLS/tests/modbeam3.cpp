// ----------------------------------------------------------------------
//
// ModBeam3.cpp - Generic three-dimensional beam model.
//
// ----------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include "modbeam3.h"
#include "utl.h"

#define     TOL      1.0e-10

/* Prototypes of Local Functions */
static void  ComputeRotMat       (double[3],double[3],double[3][3]);
static void  NormalizeVector     (double[3]);
static void  BeamTL_LocalIntLoad (double,double,double,double,double,
                                  double,double,double[12],double[12]);
static void  BeamTL_LocalStiff   (double,double,double,double,double,double,
                                  double,double[12],double[12][12]);


cModelBeam3D :: cModelBeam3D ( char *filename ) : cModel ( )
{
   int   id, i, j;
   int   num_load, num_sup;
   char mod_type[80];
   FILE *in = NULL;

   // Read data file

   in = fopen( filename, "r" );

   if ( !in )
   {
      printf ( "\n\n ### %s not defined !!! ###\n\n", filename );
      exit ( -1 );
   }

   read_string( in, mod_type );

   fscanf ( in, "%d", &_iNumNodes );     
   //printf ( "%d\n", _iNumNodes );     

   _paCoord = ( double** ) MemAlloc ( _iNumNodes, sizeof ( double* ) );
   _paDof   = ( int** )    MemAlloc ( _iNumNodes, sizeof ( int* ) );
   _paLoad  = ( double** ) MemAlloc ( _iNumNodes, sizeof ( double* ) );
   
   for ( i = 0; i < _iNumNodes; i++ )
   {
      _paCoord[i] = ( double* ) MemAlloc ( 3, sizeof ( double ) );
      _paDof[i]   = ( int* )    MemAlloc ( 6, sizeof ( int ) );
      _paLoad[i]  = ( double* ) MemAlloc ( 6, sizeof ( double ) );

      fscanf ( in, "%lf%lf%lf", &_paCoord[i][0], 
                                &_paCoord[i][1],  
                                &_paCoord[i][2] );
      for ( j = 0; j < 6; j++ )
      {
         _paDof[i][j]  = 0;
         _paLoad[i][j] = 0.0;
      }
      //printf ("%lf\t%lf\t%lf\n", _paCoord[i][0], _paCoord[i][1],  _paCoord[i][2] );
   }

   fscanf ( in, "%d", &num_sup );

   for ( i = 0; i < num_sup; i++ )
   {
      fscanf ( in, "%d", &id );
      for ( j = 0; j < 6; j++ ) fscanf ( in, "%d", &_paDof[id][j] );
   }

   fscanf ( in, "%d", &num_load );

   for ( i = 0; i < num_load; i++ )
   {
      fscanf( in, "%d", &id );
      for( j = 0; j < 6; j++ ) fscanf ( in, "%lf", &_paLoad[id][j] );
   }

   fscanf ( in, "%d", &_iNumElms );

   _paInc  = ( int** )    MemAlloc ( _iNumElms, sizeof ( int* ) );
   _paProp = ( double** ) MemAlloc ( _iNumElms, sizeof ( double* ) );

   for ( i = 0; i < _iNumElms; i++ )
   {
      _paInc[i]  = ( int* )    MemAlloc( 2, sizeof ( int ) );
      _paProp[i] = ( double* ) MemAlloc( 9, sizeof ( double ) );

      for ( j = 0; j < 2; j++ ) fscanf( in, "%d",  &_paInc[i][j] );
      for ( j = 0; j < 9; j++ ) fscanf( in, "%lf", &_paProp[i][j] );
   }

   fscanf( in, "%d", &_iNumGra );

   _paGra = ( int** ) MemAlloc( _iNumGra, sizeof(int*) );

   for ( i = 0; i < _iNumGra; i++ ) 
   {
      _paGra[i] = ( int* ) MemAlloc ( 2, sizeof ( int ) );

      for( j = 0; j < 2; j++ ) fscanf( in, "%d", &_paGra[i][j] );
   }

   fclose ( in );

   // Generate d.o.f.

   _iNumEq = 0;

   for ( i = 0; i < _iNumNodes; i++ )
      for ( j = 0; j < 6; j++ )
         _paDof[i][j] = ( _paDof[i][j] == 0 ) ? _iNumEq++ : -1;

   // Define number of strain components

   _iNumEpsEq = _iNumElms;
}

cModelBeam3D :: ~cModelBeam3D ( void )
{
   int i;

   for ( i = 0; i < _iNumNodes; i++ )
   {
      MemFree ( _paCoord[i] );
      MemFree ( _paDof[i] );
      MemFree ( _paLoad[i] );
   }
   MemFree ( _paCoord );
   MemFree ( _paDof );
   MemFree ( _paLoad );
   
   for ( i = 0; i < _iNumElms; i++ )
   {
      MemFree ( _paInc[i] );
      MemFree ( _paProp[i] );
   }
   MemFree ( _paInc );
   MemFree ( _paProp );

   for ( i = 0; i < _iNumGra; i++ ) MemFree ( _paGra[i] );
   MemFree ( _paGra );
}

void cModelBeam3D :: Init ( void )
{
   char  fn[50];
   char *dir[6] = { "u", "v", "w", "rx", "ry", "rz" };
   
   _afOut = ( FILE ** ) MemAlloc ( _iNumGra, sizeof ( FILE * ) );
   for ( int i = 0; i < _iNumGra; i++ )
   {
      sprintf( fn, "node%d%s.out", _paGra[i][0], dir[_paGra[i][1]] );
      _afOut[i] = fopen( fn, "w" );
      fprintf( _afOut[i], "%f   %f\n", 0.0, 0.0 );
   }
}

void cModelBeam3D :: Convergence ( double dFactor, double *pdSol )
{
   for ( int i = 0; i < _iNumGra; i++ )
   {
      double u = (_paDof[_paGra[i][0]][_paGra[i][1]] >= 0) ? 
                 pdSol[_paDof[_paGra[i][0]][_paGra[i][1]]] : 0.;
      fprintf( _afOut[i], "%f   %f\n", u, dFactor );
   }
}

void cModelBeam3D :: InternalVector ( double *u, double *f )
{
   int     i, j, k, l;
   double  Lo, X1, X2, Y1, Y2, Z1, Z2, EA, EIx, EIy, EIz, GAy, GAz,
           v1[3], v3[3], r[3][3], uloc[12], uelm[12], floc[12], felm[12];

   MathVecZero( _iNumEq, f );

   for( i = 0; i < _iNumElms; i++ )
   {

      X1 = _paCoord[_paInc[i][0]][0];
      X2 = _paCoord[_paInc[i][1]][0];
      Y1 = _paCoord[_paInc[i][0]][1];
      Y2 = _paCoord[_paInc[i][1]][1];
      Z1 = _paCoord[_paInc[i][0]][2];
      Z2 = _paCoord[_paInc[i][1]][2];

      Lo = sqrt ( ( X2 - X1 ) * ( X2 - X1 ) +
                  ( Y2 - Y1 ) * ( Y2 - Y1 ) +
                  ( Z2 - Z1 ) * ( Z2 - Z1 ) );

      for ( j = 0; j < 12; j++ ) 
      {
         uelm[j] = felm[j] = floc[j] = 0.0;
      }

      EA  = _paProp[i][0];
      EIx = _paProp[i][1];
      EIy = _paProp[i][2];
      EIz = _paProp[i][3];
      GAy = _paProp[i][4];
      GAz = _paProp[i][5];

      v1[0] = ( X2 - X1 ) / Lo;
      v1[1] = ( Y2 - Y1 ) / Lo;
      v1[2] = ( Z2 - Z1 ) / Lo;
      v3[0] = _paProp[i][6];
      v3[1] = _paProp[i][7];
      v3[2] = _paProp[i][8];

      /* Compute Rotation Matrix */
      ComputeRotMat( v1, v3, r );

      /* Get Element Displacements in Global Coordinates */
      for ( j = 0; j < 2; j++ )
        for ( k = 0; k < 6; k++ )
           uelm[6*j+k] = _paDof[_paInc[i][j]][k] >= 0 ? 
                         u[_paDof[_paInc[i][j]][k]] : 0.0;

      /* Transform Displacements from Global to Local Coordinates */
      for ( k = 0; k < 4; k++ ) {
         for ( l = 0; l < 3; l++ ) {
            uloc[3*k+l] = 0.0;
            for ( j = 0; j < 3; j++) {
               uloc[3*k+l] += r[l][j] * uelm[3*k+j];
            }
         }
      }

      /* Compute Local Internal Forces for Current Element */
      BeamTL_LocalIntLoad ( Lo, EA, GAy, GAz, EIx, EIy, EIz, uloc, floc );

      /* Transform Forces from Local to Global Coordinates */
      for ( k = 0; k < 4; k++ ) {
         for ( l = 0; l < 3; l++ ) {
            felm[3*k+l] = 0.0;
            for ( j = 0; j < 3; j++ )  {
               felm[3*k+l] += r[j][l] * floc[3*k+j];
            }
         }
      }

      /* Add Local Internal Forces to the Internal Force Vector */
      for ( j = 0; j < 2; j++ )
         for ( k = 0; k < 6; k++ )
            if ( _paDof[_paInc[i][j]][k] >= 0 )
               f[_paDof[_paInc[i][j]][k]] += felm[6*j+k];
   }
}

void cModelBeam3D :: Reference ( double *pdReference )
{
   int i,j;
   for ( i = 0; i < _iNumNodes; i++ )
      for ( j = 0; j < 6; j++ )
         if ( _paDof[i][j] >= 0 ) pdReference[_paDof[i][j]] = _paLoad[i][j];
}

void cModelBeam3D :: TangentMatrix ( double *u, cLinSys *kt )
{
   int     i, j, k, l, m, n;
   double  Lo, X1, X2, Y1, Y2, Z1, Z2, EA, EIx, EIy, EIz, GAy, GAz,
           v1[3], v3[3], r[3][3], uloc[12], uelm[12], kloc[12][12], 
           a[12][12], kelm[12][12];

   kt->Zero();

   for( i = 0; i < _iNumElms; i++ )
   {

      X1 = _paCoord[_paInc[i][0]][0];
      X2 = _paCoord[_paInc[i][1]][0];
      Y1 = _paCoord[_paInc[i][0]][1];
      Y2 = _paCoord[_paInc[i][1]][1];
      Z1 = _paCoord[_paInc[i][0]][2];
      Z2 = _paCoord[_paInc[i][1]][2];

      Lo = sqrt ( ( X2 - X1 ) * ( X2 - X1 ) +
                  ( Y2 - Y1 ) * ( Y2 - Y1 ) +
                  ( Z2 - Z1 ) * ( Z2 - Z1 ) );

      for ( j = 0; j < 12; j++ ) 
      {
         uelm[j] = 0.0;
         for ( k = 0; k < 12; k++ ) 
         {
            kloc[j][k] = 0.0;
            kelm[j][k] = 0.0;
         }
      }

      EA  = _paProp[i][0];
      EIx = _paProp[i][1];
      EIy = _paProp[i][2];
      EIz = _paProp[i][3];
      GAy = _paProp[i][4];
      GAz = _paProp[i][5];

      v1[0] = ( X2 - X1 ) / Lo;
      v1[1] = ( Y2 - Y1 ) / Lo;
      v1[2] = ( Z2 - Z1 ) / Lo;
      v3[0] = _paProp[i][6];
      v3[1] = _paProp[i][7];
      v3[2] = _paProp[i][8];

      /* Compute Rotation Matrix */
      ComputeRotMat( v1, v3, r );

      /* Get Element Displacements in Global Coordinates */
      for ( j = 0; j < 2; j++ )
        for ( k = 0; k < 6; k++ )
           uelm[6*j+k] = _paDof[_paInc[i][j]][k] >= 0 ? 
                         u[_paDof[_paInc[i][j]][k]] : 0.0;

      /* Transform Displacements from Global to Local Coordinates */
      for ( k = 0; k < 4; k++ ) {
         for ( l = 0; l < 3; l++ ) {
            uloc[3*k+l] = 0.0;
            for ( j = 0; j < 3; j++) {
               uloc[3*k+l] += r[l][j] * uelm[3*k+j];
            }
         }
      }

      /* Compute Local Stiffness Matrix for Current Element */
      BeamTL_LocalStiff(  Lo, EA, GAy, GAz, EIx, EIy, EIz, uloc, kloc );

      /* Transform Stiffness Matrix from Local to Global Coordinates */
      for ( m = 0; m < 12; m++ ) {
         for ( k = 0; k < 4; k++ ) {
            for ( j = 0; j < 3; j++ ) {
               a[m][3*k+j] = 0.0;
               for ( l = 0; l < 3; l++ ) {
                  a[m][3*k+j] += kloc[m][3*k+l] * r[l][j];
               }
            }
         }
      }
      for ( m = 0; m < 12; m++ ) {
         for ( k = 0; k < 4; k++ ) {
            for ( j = 0; j < 3; j++ ) {
               kelm[m][3*k+j] = 0.0;
               for ( l = 0 ; l < 3; l++ ) {
                  kelm[m][3*k+j] += r[l][j] * a[3*k+l][m];
               }
            }
         }
      }

      /* Assemble Element Matrix into Global Stiffness Matrix */
      for ( j = 0; j < 2; j++ )
      {
         for ( k = 0; k < 6; k++ )
         {
            if ( _paDof[_paInc[i][j]][k] >= 0 )
            {
               for ( m = 0; m < 2; m++ )
               {
                  for ( n = 0; n < 6; n++ )
                  {
                     if ( ( _paDof[_paInc[i][m]][n] >= 0 ) && 
                          ( _paDof[_paInc[i][m]][n] <= _paDof[_paInc[i][j]][k] ) ) 
                     {
                        kt->AddA(_paDof[_paInc[i][j]][k],_paDof[_paInc[i][m]][n],
                          kelm[6*j+k][6*m+n]);
                     }
                  }
               }
            }
         }
      }
   }

}


/* ======================= ComputeRotMat =============================== */
static void ComputeRotMat( double   u[3],
                           double   w[3],
                           double   r[3][3] )
{
   int     i, j;
   double  v[3];
   
   for(i=0;i<3;i++) {
      for(j=0;j<3;j++) {
         r[i][j] = 0.0;
      }
   }
   
   v[0] = w[1] * u[2] - w[2] * u[1];
   v[1] = w[2] * u[0] - w[0] * u[2];
   v[2] = w[0] * u[1] - w[1] * u[0];

   w[0] = u[1] * v[2] - u[2] * v[1];
   w[1] = u[2] * v[0] - u[0] * v[2];
   w[2] = u[0] * v[1] - u[1] * v[0];

   NormalizeVector( u );
   NormalizeVector( v );
   NormalizeVector( w );

   for(i=0;i<3;i++) {
      r[0][i] = u[i];
      r[1][i] = v[i];
      r[2][i] = w[i];
   }

}  /* End of ComputeRotMat  */


/* ========================== NormalizeVector ========================== */
static void NormalizeVector( double  v[3] )
{
   int     i;
   double  norm = 0.0;

   for ( i = 0; i < 3; i++ ) norm += v[i] * v[i];
   norm = sqrt ( norm );
   if ( norm > TOL ) {
      for ( i = 0; i < 3; i++ ) v[i] /= norm;
   }

}  /* End of NormalizeVector */


/* ==================== BeamTL_LocalIntLoad ============================ */
static void  BeamTL_LocalIntLoad( double   L,
                                  double   EA,
                                  double   GAy,
                                  double   GAz,
                                  double   EIx,
                                  double   EIy,
                                  double   EIz,
                                  double   U[12],
                                  double   F[12] )
{
   double t1,   t10,  t103, t11,  t112, t115, t116, t117, t118, t12, 
          t121, t122, t123, t129, t130, t132, t133, t135, t136, t138, 
          t14,  t140, t141, t145, t147, t148, t15,  t150, t152, t153, 
          t157, t16,  t160, t166, t167, t17,  t171, t173, t174, t176, 
          t177, t178, t179, t18,  t181, t185, t19,  t190, t199, t2,   
          t20,  t200, t201, t203, t206, t208, t209, t21,  t210, t213, 
          t215, t216, t217, t22,  t220, t223, t224, t225, t226, t227, 
          t23,  t230, t232, t235, t237, t238, t239, t240, t241, t242, 
          t243, t245, t247, t248, t249, t250, t251, t254, t256, t258, 
          t26,  t261, t264, t266, t267, t27,  t270, t272, t273, t275, 
          t276, t278, t279, t280, t282, t284, t286, t289, t291, t294, 
          t297, t299, t300, t302, t304, t305, t307, t32,  t321, t323, 
          t325, t327, t331, t333, t337, t34,  t347, t349, t351, t355, 
          t36,  t368, t369, t37,  t370, t371, t372, t375, t378, t379, 
          t38,  t380, t382, t384, t386, t387, t388, t390, t392, t394, 
          t398, t4,   t40,  t403, t404, t405, t406, t408, t41,  t412, 
          t414, t415, t417, t418, t420, t434, t436, t440, t442, t446, 
          t457, t459, t46,  t463, t476, t477, t478, t479, t48,  t480, 
          t483, t486, t487, t489, t49,  t491, t493, t494, t496, t497, 
          t5,   t501, t504, t505, t507, t508, t51,  t512, t514, t519, 
          t52,  t520, t57,  t6,   t60,  t63,  t65,  t68,  t7,   t70,  
          t71,  t73,  t74,  t77,  t80,  t82,  t84,  t87,  t9,   t90, 
          t94;

   t1 = U[3];
   t2 = U[9];
   t4 = t1/2.0+t2/2.0;
   t5 = t4*t4;
   t6 = U[4];
   t7 = U[10];
   t9 = t6/2.0+t7/2.0;
   t10 = t9*t9;
   t11 = U[5];
   t12 = U[11];
   t14 = t11/2.0+t12/2.0;
   t15 = t14*t14;
   t16 = t5+t10+t15;
   if( t16 < TOL ) t16 = TOL;
   t17 = sqrt(t16);
   t18 = t17/2.0;
   t19 = sin(t18);
   t20 = t19*t19;
   t21 = 1/t16;
   t22 = t20*t21;
   t23 = t10+t15;
   t26 = 1.0-2.0*t22*t23;
   t27 = 1/L;
   t32 = 1.0-t27*U[0]+t27*U[6];
   t34 = sin(t17);
   t36 = t34/t17;
   t37 = t36*t14;
   t38 = t4*t9;
   t40 = 2.0*t22*t38;
   t41 = t37+t40;
   t46 = -t27*U[1]+t27*U[7];
   t48 = t36*t9;
   t49 = t4*t14;
   t51 = 2.0*t22*t49;
   t52 = -t48+t51;
   t57 = -t27*U[2]+t27*U[8];
   t60 = EA*(t26*t32+t41*t46+t52*t57-1.0);
   t63 = -t37+t40;
   t65 = t5+t15;
   t68 = 1.0-2.0*t22*t65;
   t70 = t36*t4;
   t71 = t9*t14;
   t73 = 2.0*t22*t71;
   t74 = t70+t73;
   t77 = GAy*(t63*t32+t68*t46+t74*t57);
   t80 = t48+t51;
   t82 = -t70+t73;
   t84 = t5+t10;
   t87 = 1.0-2.0*t22*t84;
   t90 = GAz*(t80*t32+t82*t46+t87*t57);
   t94 = -2.0*t60*t26*t27-2.0*t77*t63*t27-2.0*t90*t80*t27;
   F[0] = L*t94/2.0;
   t103 = -2.0*t60*t41*t27-2.0*t77*t68*t27-2.0*t90*t82*t27;
   F[1] = L*t103/2.0;
   t112 = -2.0*t60*t52*t27-2.0*t77*t74*t27-2.0*t90*t87*t27;
   F[2] = L*t112/2.0;
   t115 = 1/t17/t16;
   t116 = t19*t115;
   t117 = cos(t18);
   t118 = t23*t117;
   t121 = t16*t16;
   t122 = 1/t121;
   t123 = t20*t122;
   t129 = cos(t17);
   t130 = t129*t21;
   t132 = t130*t49/2.0;
   t133 = t34*t115;
   t135 = t133*t49/2.0;
   t136 = t5*t9;
   t138 = t116*t136*t117;
   t140 = 2.0*t123*t136;
   t141 = t22*t9;
   t145 = t130*t38/2.0;
   t147 = t133*t38/2.0;
   t148 = t5*t14;
   t150 = t116*t148*t117;
   t152 = 2.0*t123*t148;
   t153 = t22*t14;
   t157 = t60*((-t116*t118*t4+2.0*t123*t23*t4)*t32+
          (t132-t135+t138-t140+t141)*t46+
          (-t145+t147+t150-t152+t153)*t57);
   t160 = t65*t117;
   t166 = t22*t4;
   t167 = 2.0*t166;
   t171 = t130*t5/2.0;
   t173 = t133*t5/2.0;
   t174 = t36/2.0;
   t176 = t14*t117;
   t177 = t176*t4;
   t178 = t116*t9*t177;
   t179 = t71*t4;
   t181 = 2.0*t123*t179;
   t185 = t77*((-t132+t135+t138-t140+t141)*t32+
          (-t116*t160*t4+2.0*t123*t65*t4-t167)*t46+
          (t171-t173+t174+t178-t181)*t57);
   t190 = t84*t117;
   t199 = t90*((t145-t147+t150-t152+t153)*t32+
          (-t171+t173-t174+t178-t181)*t46+
          (-t116*t190*t4+2.0*t123*t84*t4-t167)*t57);
   t200 = 1.0-t36;
   t201 = t200*t21;
   t203 = t36+t201*t5;
   t206 = -t27*t1+t27*t2;
   t208 = t201*t38;
   t209 = 2.0*t153;
   t210 = t208+t209;
   t213 = -t27*t6+t27*t7;
   t215 = t201*t49;
   t216 = 2.0*t141;
   t217 = t215-t216;
   t220 = -t27*t11+t27*t12;
   t223 = EIx*(t203*t206+t210*t213+t217*t220);
   t224 = t130*t4;
   t225 = t224/2.0;
   t226 = t133*t4;
   t227 = t226/2.0;
   t230 = (-t224/2.0+t226/2.0)*t21;
   t232 = t200*t122;
   t235 = t201*t4;
   t237 = (t225-t227+t230*t5-t232*t5*t4+t235)*t206;
   t238 = t203*t27;
   t239 = t230*t38;
   t240 = t232*t136;
   t241 = t201*t9;
   t242 = t241/2.0;
   t243 = t116*t177;
   t245 = 2.0*t123*t49;
   t247 = (t239-t240+t242+t243-t245)*t213;
   t248 = t230*t49;
   t249 = t232*t148;
   t250 = t201*t14;
   t251 = t250/2.0;
   t254 = t116*t9*t117*t4;
   t256 = 2.0*t123*t38;
   t258 = (t248-t249+t251-t254+t256)*t220;
   t261 = t208-t209;
   t264 = t36+t201*t10;
   t266 = t201*t71;
   t267 = t266+t167;
   t270 = EIy*(t261*t206+t264*t213+t267*t220);
   t272 = (t239-t240+t242-t243+t245)*t206;
   t273 = t261*t27;
   t275 = t10*t4;
   t276 = t232*t275;
   t278 = (t225-t227+t230*t10-t276)*t213;
   t279 = t230*t71;
   t280 = t232*t179;
   t282 = t116*t5*t117;
   t284 = 2.0*t123*t5;
   t286 = (t279-t280+t282-t284+t22)*t220;
   t289 = t215+t216;
   t291 = t266-t167;
   t294 = t36+t201*t15;
   t297 = EIz*(t289*t206+t291*t213+t294*t220);
   t299 = (t248-t249+t251+t254-t256)*t206;
   t300 = t289*t27;
   t302 = (t279-t280-t282+t284-t22)*t213;
   t304 = t15*t4;
   t305 = t232*t304;
   t307 = (t225-t227+t230*t15-t305)*t220;
   F[3] = L*(2.0*t157+2.0*t185+2.0*t199+2.0*t223*
          (t237-t238+t247+t258)+2.0*t270*
          (t272-t273+t278+t286)+2.0*t297*
          (t299-t300+t302+t307))/2.0;
   t321 = t130*t71/2.0;
   t323 = t133*t71/2.0;
   t325 = t116*t275*t117;
   t327 = 2.0*t123*t275;
   t331 = t130*t10/2.0;
   t333 = t133*t10/2.0;
   t337 = t60*((-t116*t118*t9+2.0*t123*t23*t9-t216)*t32+
          (t321-t323+t325-t327+t166)*t46+
          (-t331+t333-t174+t178-t181)*t57);
   t347 = t10*t14;
   t349 = t116*t347*t117;
   t351 = 2.0*t123*t347;
   t355 = t77*((-t321+t323+t325-t327+t166)*t32+
          (-t116*t160*t9+2.0*t123*t65*t9)*t46+
          (t145-t147+t349-t351+t153)*t57);
   t368 = t90*((t331-t333+t174+t178-t181)*t32+
          (-t145+t147+t349-t351+t153)*t46+
          (-t116*t190*t9+2.0*t123*t84*t9-t216)*t57);
   t369 = t130*t9;
   t370 = t369/2.0;
   t371 = t133*t9;
   t372 = t371/2.0;
   t375 = (-t369/2.0+t371/2.0)*t21;
   t378 = (t370-t372+t375*t5-t240)*t206;
   t379 = t375*t38;
   t380 = t235/2.0;
   t382 = t116*t176*t9;
   t384 = 2.0*t123*t71;
   t386 = (t379-t276+t380+t382-t384)*t213;
   t387 = t210*t27;
   t388 = t375*t49;
   t390 = t116*t10*t117;
   t392 = 2.0*t123*t10;
   t394 = (t388-t280-t390+t392-t22)*t220;
   t398 = (t379-t276+t380-t382+t384)*t206;
   t403 = (t370-t372+t375*t10-t232*t10*t9+t241)*t213;
   t404 = t264*t27;
   t405 = t375*t71;
   t406 = t232*t347;
   t408 = (t405-t406+t251+t254-t256)*t220;
   t412 = (t388-t280+t390-t392+t22)*t206;
   t414 = (t405-t406+t251-t254+t256)*t213;
   t415 = t291*t27;
   t417 = t15*t9;
   t418 = t232*t417;
   t420 = (t370-t372+t375*t15-t418)*t220;
   F[4] = L*(2.0*t337+2.0*t355+2.0*t368+2.0*t223*
          (t378+t386-t387+t394)+2.0*t270*
          (t398+t403-t404+t408)+2.0*t297*
          (t412+t414-t415+t420))/2.0;
   t434 = t130*t15/2.0;
   t436 = t133*t15/2.0;
   t440 = t116*t304*t117;
   t442 = 2.0*t123*t304;
   t446 = t60*((-t116*t118*t14+2.0*t123*t23*t14-t209)*t32+
          (t434-t436+t174+t178-t181)*t46+
          (-t321+t323+t440-t442+t166)*t57);
   t457 = t116*t417*t117;
   t459 = 2.0*t123*t417;
   t463 = t77*((-t434+t436-t174+t178-t181)*t32+
          (-t116*t160*t14+2.0*t123*t65*t14-t209)*t46+
          (t132-t135+t457-t459+t141)*t57);
   t476 = t90*((t321-t323+t440-t442+t166)*t32+
          (-t132+t135+t457-t459+t141)*t46+
          (-t116*t190*t14+2.0*t123*t84*t14)*t57);
   t477 = t130*t14;
   t478 = t477/2.0;
   t479 = t133*t14;
   t480 = t479/2.0;
   t483 = (-t477/2.0+t479/2.0)*t21;
   t486 = (t478-t480+t483*t5-t249)*t206;
   t487 = t483*t38;
   t489 = t116*t15*t117;
   t491 = 2.0*t123*t15;
   t493 = (t487-t280+t489-t491+t22)*t213;
   t494 = t483*t49;
   t496 = (t494-t305+t380-t382+t384)*t220;
   t497 = t217*t27;
   t501 = (t487-t280-t489+t491-t22)*t206;
   t504 = (t478-t480+t483*t10-t406)*t213;
   t505 = t483*t71;
   t507 = (t505-t418+t242+t243-t245)*t220;
   t508 = t267*t27;
   t512 = (t494-t305+t380+t382-t384)*t206;
   t514 = (t505-t418+t242-t243+t245)*t213;
   t519 = (t478-t480+t483*t15-t232*t15*t14+t250)*t220;
   t520 = t294*t27;
   F[5] = L*(2.0*t446+2.0*t463+2.0*t476+2.0*t223*
          (t486+t493+t496-t497)+2.0*t270*
          (t501+t504+t507-t508)+2.0*t297*
          (t512+t514+t519-t520))/2.0;
   F[6] = -L*t94/2.0;
   F[7] = -L*t103/2.0;
   F[8] = -L*t112/2.0;
   F[9] = L*(2.0*t157+2.0*t185+2.0*t199+2.0*t223*
          (t237+t238+t247+t258)+2.0*t270*
          (t272+t273+t278+t286)+2.0*t297*
          (t299+t300+t302+t307))/2.0;
   F[10] = L*(2.0*t337+2.0*t355+2.0*t368+2.0*t223*
           (t378+t386+t387+t394)+2.0*t270*
           (t398+t403+t404+t408)+2.0*t297*
           (t412+t414+t415+t420))/2.0;
   F[11] = L*(2.0*t446+2.0*t463+2.0*t476+2.0*t223*
           (t486+t493+t496+t497)+2.0*t270*
           (t501+t504+t507+t508)+2.0*t297*
           (t512+t514+t519+t520))/2.0;

   return;

}  /* End of BeamTL_LocalIntLoad */


/* ======================= BeamTL_LocalStiff =========================== */
static void  BeamTL_LocalStiff( double   L,
                                double   EA,
                                double   GAy,
                                double   GAz,
                                double   EIx,
                                double   EIy,
                                double   EIz,
                                double   U[12],
                                double   K[12][12] )
{
   double t1,    t10,   t100,  t1000, t1001, t1002, t1003, t1005, t1008, 
          t1009, t101,  t1011, t1012, t1013, t1014, t1015, t1016, t1017, 
          t1018, t102,  t1020, t1022, t1023, t1025, t1026, t1027, t1028, 
          t1031, t1032, t1033, t1035, t1037, t1038, t1039, t1040, t1041, 
          t1042, t1043, t1044, t1046, t1047, t1048, t105,  t1050, t1052, 
          t1054, t1056, t1057, t1058, t1059, t1060, t1062, t1063, t1064, 
          t1065, t1066, t1069, t1070, t1071, t1072, t1073, t1074, t1076, 
          t1077, t1078, t1079, t1081, t1082, t1083, t1084, t1085, t1086, 
          t1088, t1089, t1090, t1092, t1094, t1097, t11,   t110,  t1100, 
          t1101, t1113, t1115, t1117, t1118, t112,  t1121, t1124, t1126,
          t1127, t1129, t113,  t1131, t1132, t1133, t1136, t1137, t1138, 
          t1140, t115,  t1154, t1156, t1158, t116,  t1160, t1161, t1162, 
          t1165, t1166, t1167, t1169, t117,  t1171, t118,  t1183, t1184, 
          t1185, t1187, t1188, t119,  t1190, t1191, t1192, t1193, t1194, 
          t1195, t1196, t1197, t1198, t1199, t12,   t120,  t1200, t1202, 
          t1204, t1206, t1208, t121,  t1210, t1212, t1214, t1215, t1216, 
          t1217, t1218, t1219, t122,  t1221, t1223, t1224, t1226, t1227, 
          t1228, t1229, t123,  t1230, t1231, t1232, t1233, t1234, t1235, 
          t1238, t1239, t124,  t1241, t1242, t1243, t1244, t1245, t1246, 
          t1247, t1249, t125,  t1250, t1251, t1253, t1255, t1256, t1257,
          t1258, t1259, t1260, t1261, t1262, t1265, t1266, t1267, t1268, 
          t1270, t1272, t1273, t1274, t1275, t1277, t1278, t1279, t1280, 
          t1281, t1283, t1285, t1287, t1289, t1290, t1293, t130,  t1302, 
          t1306, t1310, t1314, t1317, t1321, t1325, t1329, t133,  t1332, 
          t1336, t134,  t1340, t1344, t1347, t1348, t1349, t135,  t1355, 
          t136,  t1360, t1364, t1366, t1368, t137,  t1370, t1373, t1375, 
          t1377, t1378, t138,  t1381, t1383, t1384, t1386, t1387, t1388, 
          t139,  t1391, t1392, t1393, t1394, t1396, t14,   t140,  t1408, 
          t141,  t1410, t1413, t1414, t1416, t1418, t1419, t142,  t1420, 
          t1423, t1424, t1425, t1426, t1428, t1430, t1443, t1444, t1446,
          t1448, t1450, t1452, t1456, t1457, t1458, t1460, t1462, t1463, 
          t1465, t1466, t1467, t1469, t147,  t1470, t1471, t1473, t1476, 
          t1477, t1479, t1481, t1482, t1484, t1487, t1489, t149,  t1490, 
          t1495, t1499, t15,   t150,  t1500, t1502, t1503, t1504, t1505, 
          t1507, t1508, t1509, t1510, t1513, t1516, t1517, t1518, t1520, 
          t1522, t1524, t1526, t1528, t1531, t1535, t1536, t1550, t1552, 
          t1554, t1555, t1557, t156,  t1560, t1562, t1563, t1565, t1577, 
          t1580, t1582, t1584, t1585, t1588, t159,  t1590, t1591, t1593, 
          t1595, t16,   t1607, t1608, t161,  t1610, t1612, t1614, t1616, 
          t1619, t1620, t1621, t1623, t1625, t1626, t1628, t1629, t1630,
          t1631, t1632, t1633, t1634, t1635, t1636, t1639, t164,  t1641, 
          t1642, t1646, t1647, t1648, t1649, t1650, t1651, t1652, t1655, 
          t1657, t1658, t1659, t1660, t1661, t1663, t1665, t1667, t1669, 
          t167,  t1670, t1673, t168,  t169,  t1691, t17,   t1703, t1715, 
          t1718, t1719, t172,  t1720, t1726, t173,  t1731, t1735, t1737, 
          t1738, t174,  t1740, t1741, t1742, t1745, t1748, t175,  t1750, 
          t1752, t1753, t1756, t1757, t1758, t1759, t1761, t177,  t1774, 
          t1777, t1779, t178,  t1781, t1782, t1785, t1786, t1787, t1788, 
          t179,  t1790, t1792, t18,   t180,  t1805, t1806, t1808, t181, 
          t1810, t1812, t1814, t1818, t1819, t182,  t1820, t1822, t1825,
          t1826, t1828, t183,  t1830, t1831, t1833, t1834, t1835, t1836, 
          t1837, t1839, t1842, t1845, t1849, t185,  t1850, t1851, t1852, 
          t1853, t1855, t1858, t186,  t1860, t1861, t1862, t1863, t1868, 
          t1872, t1873, t1875, t1878, t1887, t1891, t1895, t1899, t19, 
          t1911, t192,  t1923, t1926, t1930, t1934, t1938, t1941, t1945, 
          t1949, t195,  t1953, t1965, t1968, t197,  t1972, t1976, t1980, 
          t199,  t1995, t1998, t2,    t20,   t2002, t2006, t2010, t202, 
          t205,  t207,  t208,  t21,   t214,  t218,  t22,   t222,  t225, 
          t226,  t229,  t23,   t230,  t231,  t232,  t233,  t234,  t235, 
          t236,  t237,  t240,  t241,  t242,  t243,  t245,  t246,  t250, 
          t254,  t257,  t259,  t26,   t260,  t261,  t262,  t263,  t264, 
          t266,  t267,  t27,   t271,  t273,  t277,  t280,  t282,  t283, 
          t288,  t29,   t295,  t296,  t299,  t30,   t300,  t301,  t302, 
          t304,  t305,  t306,  t307,  t308,  t309,  t311,  t312,  t316, 
          t32,   t323,  t325,  t326,  t327,  t328,  t329,  t330,  t332, 
          t333,  t337,  t339,  t34,   t346,  t348,  t349,  t35,   t354, 
          t36,   t362,  t365,  t368,  t372,  t374,  t376,  t378,  t38, 
          t381,  t39,   t393,  t4,    t40,   t405,  t417,  t423,  t426, 
          t429,  t43,   t433,  t435,  t439,  t44,   t443,  t448,  t46, 
          t460,  t47,   t472,  t476,  t477,  t478,  t479,  t48,   t480, 
          t484,  t485,  t492,  t494,  t495,  t498,  t5,    t502,  t503, 
          t505,  t506,  t507,  t508,  t510,  t511,  t512,  t513,  t514, 
          t516,  t519,  t52,   t520,  t521,  t522,  t523,  t525,  t527, 
          t528,  t529,  t530,  t533,  t535,  t536,  t537,  t539,  t54, 
          t540,  t541,  t542,  t544,  t547,  t548,  t549,  t55,   t551, 
          t553,  t554,  t555,  t556,  t559,  t560,  t561,  t562,  t564, 
          t571,  t572,  t576,  t579,  t58,   t580,  t581,  t585,  t587, 
          t588,  t589,  t59,   t591,  t592,  t593,  t594,  t596,  t597, 
          t599,  t6,    t601,  t602,  t603,  t604,  t606,  t607,  t608, 
          t611,  t612,  t613,  t614,  t616,  t618,  t62,   t627,  t630, 
          t634,  t635,  t636,  t638,  t639,  t641,  t642,  t644,  t645, 
          t646,  t649,  t65,   t650,  t652,  t653,  t654,  t655,  t656, 
          t657,  t658,  t659,  t66,   t662,  t663,  t664,  t665,  t666, 
          t667,  t668,  t669,  t67,   t672,  t673,  t674,  t675,  t678, 
          t679,  t681,  t682,  t685,  t686,  t688,  t689,  t69,   t691, 
          t692,  t694,  t696,  t699,  t7,    t70,   t701,  t702,  t705, 
          t707,  t708,  t709,  t711,  t712,  t713,  t714,  t715,  t717, 
          t718,  t719,  t721,  t723,  t724,  t725,  t726,  t728,  t729, 
          t730,  t731,  t732,  t733,  t734,  t735,  t737,  t738,  t739, 
          t74,   t741,  t743,  t744,  t745,  t746,  t748,  t749,  t750, 
          t751,  t754,  t755,  t756,  t757,  t759,  t76,   t760,  t761, 
          t762,  t763,  t764,  t765,  t766,  t767,  t768,  t77,   t772, 
          t774,  t775,  t778,  t779,  t780,  t782,  t784,  t786,  t788, 
          t789,  t79,   t790,  t792,  t793,  t794,  t795,  t797,  t798, 
          t799,  t80,   t801,  t804,  t806,  t807,  t809,  t811,  t812, 
          t813,  t815,  t818,  t819,  t82,   t820,  t821,  t822,  t823, 
          t825,  t826,  t827,  t828,  t829,  t832,  t835,  t838,  t839, 
          t840,  t842,  t844,  t846,  t848,  t85,   t850,  t851,  t852, 
          t854,  t857,  t86,   t860,  t861,  t862,  t865,  t875,  t877, 
          t879,  t881,  t884,  t886,  t887,  t889,  t89,   t891,  t892, 
          t893,  t894,  t895,  t898,  t9,    t900,  t902,  t903,  t904, 
          t905,  t908,  t909,  t91,   t910,  t911,  t913,  t916,  t925, 
          t927,  t93,   t930,  t932,  t934,  t935,  t936,  t939,  t94, 
          t940,  t941,  t942,  t944,  t946,  t949,  t95,   t96,   t961, 
          t962,  t963,  t964,  t966,  t967,  t969,  t970,  t971,  t972, 
          t973,  t974,  t975,  t976,  t977,  t978,  t979,  t980,  t981, 
          t983,  t985,  t987,  t989,  t99,   t991,  t993,  t995,  t996, 
          t997,  t998,  t999;

   t1 = U[3];
   t2 = U[9];
   t4 = t1/2.0+t2/2.0;
   t5 = t4*t4;
   t6 = U[4];
   t7 = U[10];
   t9 = t6/2.0+t7/2.0;
   t10 = t9*t9;
   t11 = U[5];
   t12 = U[11];
   t14 = t11/2.0+t12/2.0;
   t15 = t14*t14;
   t16 = t5+t10+t15;
   if( t16 < TOL ) t16 = TOL;
   t17 = sqrt(t16);
   t18 = t17/2.0;
   t19 = sin(t18);
   t20 = t19*t19;
   t21 = 1/t16;
   t22 = t20*t21;
   t23 = t10+t15;
   t26 = 1.0-2.0*t22*t23;
   t27 = t26*t26;
   t29 = L*L;
   t30 = 1/t29;
   t32 = sin(t17);
   t34 = t32/t17;
   t35 = t34*t14;
   t36 = t4*t9;
   t38 = 2.0*t22*t36;
   t39 = -t35+t38;
   t40 = t39*t39;
   t43 = t34*t9;
   t44 = t4*t14;
   t46 = 2.0*t22*t44;
   t47 = t43+t46;
   t48 = t47*t47;
   t52 = 2.0*EA*t27*t30+2.0*GAy*t40*t30+2.0*GAz*t48*t30;
   K[0][0] = L*t52/2.0;
   t54 = EA*t26;
   t55 = t35+t38;
   t58 = GAy*t39;
   t59 = t5+t15;
   t62 = 1.0-2.0*t22*t59;
   t65 = GAz*t47;
   t66 = t34*t4;
   t67 = t9*t14;
   t69 = 2.0*t22*t67;
   t70 = -t66+t69;
   t74 = 2.0*t54*t30*t55+2.0*t58*t30*t62+2.0*t65*t30*t70;
   K[0][1] = L*t74/2.0;
   t76 = -t43+t46;
   t77 = t30*t76;
   t79 = t66+t69;
   t80 = t30*t79;
   t82 = t5+t10;
   t85 = 1.0-2.0*t22*t82;
   t86 = t30*t85;
   t89 = 2.0*t54*t77+2.0*t58*t80+2.0*t65*t86;
   K[0][2] = L*t89/2.0;
   t91 = 1/L;
   t93 = 1/t17/t16;
   t94 = t19*t93;
   t95 = cos(t18);
   t96 = t23*t95;
   t99 = t16*t16;
   t100 = 1/t99;
   t101 = t20*t100;
   t102 = t23*t4;
   t105 = -t94*t96*t4+2.0*t101*t102;
   t110 = 1.0-t91*U[0]+t91*U[6];
   t112 = cos(t17);
   t113 = t112*t21;
   t115 = t113*t44/2.0;
   t116 = t32*t93;
   t117 = t116*t44;
   t118 = t117/2.0;
   t119 = t5*t9;
   t120 = t119*t95;
   t121 = t94*t120;
   t122 = t101*t119;
   t123 = 2.0*t122;
   t124 = t22*t9;
   t125 = t115-t118+t121-t123+t124;
   t130 = -t91*U[1]+t91*U[7];
   t133 = t113*t36/2.0;
   t134 = t116*t36;
   t135 = t134/2.0;
   t136 = t5*t14;
   t137 = t136*t95;
   t138 = t94*t137;
   t139 = t101*t136;
   t140 = 2.0*t139;
   t141 = t22*t14;
   t142 = -t133+t135+t138-t140+t141;
   t147 = -t91*U[2]+t91*U[8];
   t149 = t105*t110+t125*t130+t142*t147;
   t150 = t91*t149;
   t156 = EA*(t26*t110+t55*t130+t76*t147-1.0);
   t159 = -t115+t118+t121-t123+t124;
   t161 = t59*t95;
   t164 = t59*t4;
   t167 = t22*t4;
   t168 = 2.0*t167;
   t169 = -t94*t161*t4+2.0*t101*t164-t168;
   t172 = t113*t5/2.0;
   t173 = t116*t5;
   t174 = t173/2.0;
   t175 = t34/2.0;
   t177 = t14*t95;
   t178 = t177*t4;
   t179 = t94*t9*t178;
   t180 = t67*t4;
   t181 = t101*t180;
   t182 = 2.0*t181;
   t183 = t172-t174+t175+t179-t182;
   t185 = t159*t110+t169*t130+t183*t147;
   t186 = t91*t185;
   t192 = GAy*(t39*t110+t62*t130+t79*t147);
   t195 = t133-t135+t138-t140+t141;
   t197 = -t172+t174-t175+t179-t182;
   t199 = t82*t95;
   t202 = t82*t4;
   t205 = -t94*t199*t4+2.0*t101*t202-t168;
   t207 = t195*t110+t197*t130+t205*t147;
   t208 = t91*t207;
   t214 = GAz*(t47*t110+t70*t130+t85*t147);
   t218 = -2.0*t54*t150-2.0*t156*t105*t91-2.0*t58*t186-
          2.0*t192*t159*t91-2.0*t65*t208-2.0*t214*t195*t91;
   K[0][3] = L*t218/2.0;
   t222 = t23*t9;
   t225 = 2.0*t124;
   t226 = -t94*t96*t9+2.0*t101*t222-t225;
   t229 = t113*t67/2.0;
   t230 = t116*t67;
   t231 = t230/2.0;
   t232 = t10*t4;
   t233 = t232*t95;
   t234 = t94*t233;
   t235 = t101*t232;
   t236 = 2.0*t235;
   t237 = t229-t231+t234-t236+t167;
   t240 = t113*t10/2.0;
   t241 = t116*t10;
   t242 = t241/2.0;
   t243 = -t240+t242-t175+t179-t182;
   t245 = t226*t110+t237*t130+t243*t147;
   t246 = t91*t245;
   t250 = -t229+t231+t234-t236+t167;
   t254 = t59*t9;
   t257 = -t94*t161*t9+2.0*t101*t254;
   t259 = t10*t14;
   t260 = t259*t95;
   t261 = t94*t260;
   t262 = t101*t259;
   t263 = 2.0*t262;
   t264 = t133-t135+t261-t263+t141;
   t266 = t250*t110+t257*t130+t264*t147;
   t267 = t91*t266;
   t271 = t240-t242+t175+t179-t182;
   t273 = -t133+t135+t261-t263+t141;
   t277 = t82*t9;
   t280 = -t94*t199*t9+2.0*t101*t277-t225;
   t282 = t271*t110+t273*t130+t280*t147;
   t283 = t91*t282;
   t288 = -2.0*t54*t246-2.0*t156*t226*t91-2.0*t58*t267-
          2.0*t192*t250*t91-2.0*t65*t283-2.0*t214*t271*t91;
   K[0][4] = L*t288/2.0;
   t295 = 2.0*t141;
   t296 = -t94*t96*t14+2.0*t101*t23*t14-t295;
   t299 = t113*t15/2.0;
   t300 = t116*t15;
   t301 = t300/2.0;
   t302 = t299-t301+t175+t179-t182;
   t304 = t15*t4;
   t305 = t304*t95;
   t306 = t94*t305;
   t307 = t101*t304;
   t308 = 2.0*t307;
   t309 = -t229+t231+t306-t308+t167;
   t311 = t296*t110+t302*t130+t309*t147;
   t312 = t91*t311;
   t316 = -t299+t301-t175+t179-t182;
   t323 = -t94*t161*t14+2.0*t101*t59*t14-t295;
   t325 = t15*t9;
   t326 = t325*t95;
   t327 = t94*t326;
   t328 = t101*t325;
   t329 = 2.0*t328;
   t330 = t115-t118+t327-t329+t124;
   t332 = t316*t110+t323*t130+t330*t147;
   t333 = t91*t332;
   t337 = t229-t231+t306-t308+t167;
   t339 = -t115+t118+t327-t329+t124;
   t346 = -t94*t199*t14+2.0*t101*t82*t14;
   t348 = t337*t110+t339*t130+t346*t147;
   t349 = t91*t348;
   t354 = -2.0*t54*t312-2.0*t156*t296*t91-2.0*t58*t333-
          2.0*t192*t316*t91-2.0*t65*t349-2.0*t214*t337*t91;
   K[0][5] = L*t354/2.0;
   K[0][6] = -L*t52/2.0;
   K[0][7] = -L*t74/2.0;
   K[0][8] = -L*t89/2.0;
   K[0][9] = K[0][3];
   K[0][10] = K[0][4];
   K[0][11] = K[0][5];
   K[1][0] = K[0][1];
   t362 = t55*t55;
   t365 = t62*t62;
   t368 = t70*t70;
   t372 = 2.0*EA*t362*t30+2.0*GAy*t365*t30+2.0*GAz*t368*t30;
   K[1][1] = L*t372/2.0;
   t374 = EA*t55;
   t376 = GAy*t62;
   t378 = GAz*t70;
   t381 = 2.0*t374*t77+2.0*t376*t80+2.0*t378*t86;
   K[1][2] = L*t381/2.0;
   t393 = -2.0*t374*t150-2.0*t156*t125*t91-2.0*t376*t186-
          2.0*t192*t169*t91-2.0*t378*t208-2.0*t214*t197*t91;
   K[1][3] = L*t393/2.0;
   t405 = -2.0*t374*t246-2.0*t156*t237*t91-2.0*t376*t267-
          2.0*t192*t257*t91-2.0*t378*t283-2.0*t214*t273*t91;
   K[1][4] = L*t405/2.0;
   t417 = -2.0*t374*t312-2.0*t156*t302*t91-2.0*t376*t333-
          2.0*t192*t323*t91-2.0*t378*t349-2.0*t214*t339*t91;
   K[1][5] = L*t417/2.0;
   K[1][6] = K[0][7];
   K[1][7] = -L*t372/2.0;
   K[1][8] = -L*t381/2.0;
   K[1][9] = K[1][3];
   K[1][10] = K[1][4];
   K[1][11] = K[1][5];
   K[2][0] = K[0][2];
   K[2][1] = K[1][2];
   t423 = t76*t76;
   t426 = t79*t79;
   t429 = t85*t85;
   t433 = 2.0*EA*t423*t30+2.0*GAy*t426*t30+2.0*GAz*t429*t30;
   K[2][2] = L*t433/2.0;
   t435 = EA*t76;
   t439 = GAy*t79;
   t443 = GAz*t85;
   t448 = -2.0*t435*t150-2.0*t156*t142*t91-2.0*t439*t186-
          2.0*t192*t183*t91-2.0*t443*t208-2.0*t214*t205*t91;
   K[2][3] = L*t448/2.0;
   t460 = -2.0*t435*t246-2.0*t156*t243*t91-2.0*t439*t267-
          2.0*t192*t264*t91-2.0*t443*t283-2.0*t214*t280*t91;
   K[2][4] = L*t460/2.0;
   t472 = -2.0*t435*t312-2.0*t156*t309*t91-2.0*t439*t333-
          2.0*t192*t330*t91-2.0*t443*t349-2.0*t214*t346*t91;
   K[2][5] = L*t472/2.0;
   K[2][6] = K[0][8];
   K[2][7] = K[1][8];
   K[2][8] = -L*t433/2.0;
   K[2][9] = K[2][3];
   K[2][10] = K[2][4];
   K[2][11] = K[2][5];
   K[3][0] = K[0][9];
   K[3][1] = K[1][9];
   K[3][2] = K[2][9];
   t476 = t149*t149;
   t477 = EA*t476;
   t478 = t95*t95;
   t479 = t478*t100;
   t480 = t5*t23;
   t484 = 1/t17/t99;
   t485 = t19*t484;
   t492 = t94*t96/2.0;
   t494 = 1/t99/t16;
   t495 = t20*t494;
   t498 = t101*t23;
   t502 = t116*t136/4.0;
   t503 = t112*t100;
   t505 = 3.0/4.0*t503*t136;
   t506 = t113*t14;
   t507 = t506/4.0;
   t508 = t32*t484;
   t510 = 3.0/4.0*t508*t136;
   t511 = t116*t14;
   t512 = t511/4.0;
   t513 = t5*t4;
   t514 = t513*t9;
   t516 = t479*t514/4.0;
   t519 = 5.0/2.0*t485*t514*t95;
   t520 = t9*t95;
   t521 = t520*t4;
   t522 = t94*t521;
   t523 = 3.0/2.0*t522;
   t525 = t101*t514/4.0;
   t527 = 4.0*t495*t514;
   t528 = t101*t36;
   t529 = 3.0*t528;
   t530 = -t502-t505+t507+t510-t512+t516-
          t519+t523-t525+t527-t529;
   t533 = t116*t119/4.0;
   t535 = 3.0/4.0*t503*t119;
   t536 = t113*t9;
   t537 = t536/4.0;
   t539 = 3.0/4.0*t508*t119;
   t540 = t116*t9;
   t541 = t540/4.0;
   t542 = t513*t14;
   t544 = t479*t542/4.0;
   t547 = 5.0/2.0*t485*t542*t95;
   t548 = t94*t178;
   t549 = 3.0/2.0*t548;
   t551 = t101*t542/4.0;
   t553 = 4.0*t495*t542;
   t554 = t101*t44;
   t555 = 3.0*t554;
   t556 = t533+t535-t537-t539+t541+t544-
          t547+t549-t551+t553-t555;
   t559 = t156*((-t479*t480/4.0+5.0/2.0*t485*t96*t5+t101*
          t480/4.0-t492-4.0*t495*t480+t498)*t110+t530*
          t130+t556*t147);
   t560 = t185*t185;
   t561 = GAy*t560;
   t562 = t502+t505-t507-t510+t512+t516-
          t519+t523-t525+t527-t529;
   t564 = t5*t59;
   t571 = t94*t5*t95;
   t572 = 2.0*t571;
   t576 = t94*t161/2.0;
   t579 = t101*t5;
   t580 = 4.0*t579;
   t581 = t101*t59;
   t585 = t116*t513/4.0;
   t587 = 3.0/4.0*t503*t513;
   t588 = t113*t4;
   t589 = 3.0/4.0*t588;
   t591 = 3.0/4.0*t508*t513;
   t592 = t116*t4;
   t593 = 3.0/4.0*t592;
   t594 = t119*t14;
   t596 = t479*t594/4.0;
   t597 = t485*t9;
   t599 = 5.0/2.0*t597*t137;
   t601 = t101*t594/4.0;
   t602 = t177*t9;
   t603 = t94*t602;
   t604 = t603/2.0;
   t606 = 4.0*t495*t594;
   t607 = t101*t67;
   t608 = -t585-t587+t589+t591-t593+t596-
          t599-t601+t604+t606-t607;
   t611 = t192*(t562*t110+(-t479*t564/4.0+5.0/2.0*
          t485*t161*t5-t572+t101*t564/4.0-t576-
          4.0*t495*t564+t580+t581-t22)*t130+t608*t147);
   t612 = t207*t207;
   t613 = GAz*t612;
   t614 = -t533-t535+t537+t539-t541+t544-
          t547+t549-t551+t553-t555;
   t616 = t585+t587-t589-t591+t593+t596-
          t599-t601+t604+t606-t607;
   t618 = t5*t82;
   t627 = t94*t199/2.0;
   t630 = t101*t82;
   t634 = t214*(t614*t110+t616*t130+(-t479*t618/4.0+5.0/2.0*
          t485*t199*t5-t572+t101*t618/4.0-t627-4.0*t495*
          t618+t580+t630-t22)*t147);
   t635 = t588/2.0;
   t636 = t592/2.0;
   t638 = -t588/2.0+t592/2.0;
   t639 = t638*t21;
   t641 = 1.0-t34;
   t642 = t641*t100;
   t644 = t641*t21;
   t645 = t644*t4;
   t646 = t635-t636+t639*t5-t642*t513+t645;
   t649 = -t91*t1+t91*t2;
   t650 = t646*t649;
   t652 = t34+t644*t5;
   t653 = t652*t91;
   t654 = t639*t36;
   t655 = t642*t119;
   t656 = t644*t9;
   t657 = t656/2.0;
   t658 = 2.0*t554;
   t659 = t654-t655+t657+t548-t658;
   t662 = -t91*t6+t91*t7;
   t663 = t659*t662;
   t664 = t639*t44;
   t665 = t642*t136;
   t666 = t644*t14;
   t667 = t666/2.0;
   t668 = 2.0*t528;
   t669 = t664-t665+t667-t522+t668;
   t672 = -t91*t11+t91*t12;
   t673 = t669*t672;
   t674 = t650-t653+t663+t673;
   t675 = t674*t674;
   t678 = t644*t36;
   t679 = t678+t295;
   t681 = t644*t44;
   t682 = t681-t225;
   t685 = EIx*(t652*t649+t679*t662+t682*t672);
   t686 = t173/4.0;
   t688 = 3.0/4.0*t503*t5;
   t689 = t113/4.0;
   t691 = 3.0/4.0*t508*t5;
   t692 = t116/4.0;
   t694 = (t686+t688-t689-t691+t692)*t21;
   t696 = t638*t100;
   t699 = t639*t4;
   t701 = t641*t494;
   t702 = t5*t5;
   t705 = t642*t5;
   t707 = t644/2.0;
   t708 = -t686-t688+t689+t691-t692+t694*t5-2.0*t696*
          t513+2.0*t699+2.0*t701*t702-5.0/2.0*t705+t707;
   t709 = t708*t649;
   t711 = 2.0*t646*t91;
   t712 = t694*t36;
   t713 = t696*t119;
   t714 = 2.0*t713;
   t715 = t639*t9;
   t717 = 2.0*t701*t514;
   t718 = t642*t36;
   t719 = 3.0/2.0*t718;
   t721 = t479*t136/4.0;
   t723 = 5.0/2.0*t485*t137;
   t724 = t139/4.0;
   t725 = t94*t177;
   t726 = t725/2.0;
   t728 = 4.0*t495*t136;
   t729 = t101*t14;
   t730 = t712-t714+t715+t717-t719+t721-
          t723-t724+t726+t728-t729;
   t731 = t730*t662;
   t732 = t694*t44;
   t733 = t696*t136;
   t734 = 2.0*t733;
   t735 = t639*t14;
   t737 = 2.0*t701*t542;
   t738 = t642*t44;
   t739 = 3.0/2.0*t738;
   t741 = t479*t119/4.0;
   t743 = 5.0/2.0*t485*t120;
   t744 = t122/4.0;
   t745 = t94*t520;
   t746 = t745/2.0;
   t748 = 4.0*t495*t119;
   t749 = t101*t9;
   t750 = t732-t734+t735+t737-t739-t741+t743+
          t744-t746-t748+t749;
   t751 = t750*t672;
   t754 = t654-t655+t657-t548+t658;
   t755 = t754*t649;
   t756 = t678-t295;
   t757 = t756*t91;
   t759 = t642*t232;
   t760 = t635-t636+t639*t10-t759;
   t761 = t760*t662;
   t762 = t639*t67;
   t763 = t642*t180;
   t764 = 2.0*t579;
   t765 = t762-t763+t571-t764+t22;
   t766 = t765*t672;
   t767 = t755-t757+t761+t766;
   t768 = t767*t767;
   t772 = t34+t644*t10;
   t774 = t644*t67;
   t775 = t774+t168;
   t778 = EIy*(t756*t649+t772*t662+t775*t672);
   t779 = t712-t714+t715+t717-t719-t721+t723+
          t724-t726-t728+t729;
   t780 = t779*t649;
   t782 = 2.0*t754*t91;
   t784 = t696*t232;
   t786 = t10*t5;
   t788 = 2.0*t701*t786;
   t789 = t642*t10;
   t790 = t789/2.0;
   t792 = (-t686-t688+t689+t691-t692+t694*t10-
          2.0*t784+t788-t790)*t662;
   t793 = t694*t67;
   t794 = t696*t180;
   t795 = 2.0*t794;
   t797 = 2.0*t701*t594;
   t798 = t642*t67;
   t799 = t798/2.0;
   t801 = t479*t513/4.0;
   t804 = 5.0/2.0*t485*t513*t95;
   t806 = t94*t4*t95;
   t807 = 3.0/2.0*t806;
   t809 = t101*t513/4.0;
   t811 = 4.0*t495*t513;
   t812 = t101*t4;
   t813 = 3.0*t812;
   t815 = (t793-t795+t797-t799+t801-t804+t807-
          t809+t811-t813)*t672;
   t818 = t664-t665+t667+t522-t668;
   t819 = t818*t649;
   t820 = t681+t225;
   t821 = t820*t91;
   t822 = t762-t763-t571+t764-t22;
   t823 = t822*t662;
   t825 = t642*t304;
   t826 = t635-t636+t639*t15-t825;
   t827 = t826*t672;
   t828 = t819-t821+t823+t827;
   t829 = t828*t828;
   t832 = t774-t168;
   t835 = t34+t644*t15;
   t838 = EIz*(t820*t649+t832*t662+t835*t672);
   t839 = t732-t734+t735+t737-t739+t741-
          t743-t744+t746+t748-t749;
   t840 = t839*t649;
   t842 = 2.0*t818*t91;
   t844 = (t793-t795+t797-t799-t801+t804-
          t807+t809-t811+t813)*t662;
   t846 = t696*t304;
   t848 = t15*t5;
   t850 = 2.0*t701*t848;
   t851 = t642*t15;
   t852 = t851/2.0;
   t854 = (-t686-t688+t689+t691-t692+t694*t15-
          2.0*t846+t850-t852)*t672;
   t857 = t477+t559+t561+t611+t613+t634+EIx*t675+
          t685*(t709-t711+t731+t751)+EIy*t768+t778*
          (t780-t782+t792+t815)+EIz*t829+t838*
          (t840-t842+t844+t854);
   K[3][3] = L*t857;
   t860 = EA*t149;
   t861 = t860*t245;
   t862 = t102*t9;
   t865 = t485*t23;
   t875 = t116*t180/4.0;
   t877 = 3.0/4.0*t503*t180;
   t879 = 3.0/4.0*t508*t180;
   t881 = t479*t786/4.0;
   t884 = 5.0/2.0*t485*t786*t95;
   t886 = t94*t10*t95;
   t887 = t886/2.0;
   t889 = t101*t786/4.0;
   t891 = 4.0*t495*t786;
   t892 = t101*t10;
   t893 = t571/2.0;
   t894 = t22/2.0;
   t895 = -t875-t877+t879+t881-t884+t887-t889+t891-
          t892+t893-t579+t894;
   t898 = t116*t232/4.0;
   t900 = 3.0/4.0*t503*t232;
   t902 = 3.0/4.0*t508*t232;
   t903 = t588/4.0;
   t904 = t592/4.0;
   t905 = t898+t900-t902-t903+t904+t596-
          t599-t601+t604+t606-t607;
   t908 = t156*((-t479*t862/4.0+5.0/2.0*t865*t521+
          t101*t862/4.0-4.0*t495*t862-t522+t668)*
          t110+t895*t130+t905*t147);
   t909 = GAy*t185;
   t910 = t909*t266;
   t911 = t875+t877-t879+t881-t884+t887-t889+
          t891-t892+t893-t579+t894;
   t913 = t164*t9;
   t916 = t485*t59;
   t925 = t232*t14;
   t927 = t479*t925/4.0;
   t930 = 5.0/2.0*t485*t10*t178;
   t932 = t101*t925/4.0;
   t934 = 4.0*t495*t925;
   t935 = t548/2.0;
   t936 = -t533-t535+t537+t539-t541+t927-
          t930-t932+t934+t935-t554;
   t939 = t192*(t911*t110+(-t479*t913/4.0+5.0/2.0*
          t916*t521-t522+t101*t913/4.0-4.0*t495*
          t913+t668)*t130+t936*t147);
   t940 = GAz*t207;
   t941 = t940*t282;
   t942 = -t898-t900+t902+t903-t904+t596-t599-
          t601+t604+t606-t607;
   t944 = t533+t535-t537-t539+t541+t927-t930-
          t932+t934+t935-t554;
   t946 = t202*t9;
   t949 = t485*t82;
   t961 = t214*(t942*t110+t944*t130+
          (-t479*t946/4.0+5.0/2.0*t949*t521-2.0*t522+t101*
          t946/4.0-4.0*t495*t946+4.0*t528)*t147);
   t962 = EIx*t674;
   t963 = t536/2.0;
   t964 = t540/2.0;
   t966 = -t536/2.0+t540/2.0;
   t967 = t966*t21;
   t969 = t963-t964+t967*t5-t655;
   t970 = t969*t649;
   t971 = t967*t36;
   t972 = t645/2.0;
   t973 = 2.0*t607;
   t974 = t971-t759+t972+t603-t973;
   t975 = t974*t662;
   t976 = t679*t91;
   t977 = t967*t44;
   t978 = 2.0*t892;
   t979 = t977-t763-t886+t978-t22;
   t980 = t979*t672;
   t981 = t970+t975-t976+t980;
   t983 = t134/4.0;
   t985 = 3.0/4.0*t503*t36;
   t987 = 3.0/4.0*t508*t36;
   t989 = (t983+t985-t987)*t21;
   t991 = t966*t100;
   t993 = t967*t4;
   t995 = (-t983-t985+t987+t989*t5-t991*t513+
          t993-t713+t717-t718)*t649;
   t996 = t969*t91;
   t997 = t989*t36;
   t998 = t991*t119;
   t999 = t967*t9;
   t1000 = t999/2.0;
   t1001 = t699/2.0;
   t1002 = t705/2.0;
   t1003 = t644/4.0;
   t1005 = t479*t180/4.0;
   t1008 = 5.0/2.0*t485*t14*t521;
   t1009 = t181/4.0;
   t1011 = 4.0*t495*t180;
   t1012 = t997-t998+t1000-t784+t788-t790+t1001-
           t1002+t1003+t1005-t1008-t1009+t1011;
   t1013 = t1012*t662;
   t1014 = t659*t91;
   t1015 = t989*t44;
   t1016 = t991*t136;
   t1017 = t967*t14;
   t1018 = t1017/2.0;
   t1020 = t479*t232/4.0;
   t1022 = 5.0/2.0*t485*t233;
   t1023 = t235/4.0;
   t1025 = 4.0*t495*t232;
   t1026 = t806/2.0;
   t1027 = t1015-t1016+t1018-t794+t797-t799-t1020+
           t1022+t1023-t1025-t1026+t812;
   t1028 = t1027*t672;
   t1031 = EIy*t767;
   t1032 = t971-t759+t972-t603+t973;
   t1033 = t1032*t649;
   t1035 = t10*t9;
   t1037 = t963-t964+t967*t10-t642*t1035+t656;
   t1038 = t1037*t662;
   t1039 = t772*t91;
   t1040 = t967*t67;
   t1041 = t642*t259;
   t1042 = t1040-t1041+t667+t522-t668;
   t1043 = t1042*t672;
   t1044 = t1033+t1038-t1039+t1043;
   t1046 = t997-t998+t1000-t784+t788-t790+t1001-
           t1002+t1003-t1005+t1008+t1009-t1011;
   t1047 = t1046*t649;
   t1048 = t1032*t91;
   t1050 = t991*t232;
   t1052 = t1035*t4;
   t1054 = 2.0*t701*t1052;
   t1056 = (-t983-t985+t987+t989*t10-t1050-t696*t1035+
           t1054+t715-t718)*t662;
   t1057 = t760*t91;
   t1058 = t989*t67;
   t1059 = t991*t180;
   t1060 = t696*t259;
   t1062 = 2.0*t701*t925;
   t1063 = t735/2.0;
   t1064 = t738/2.0;
   t1065 = t1058-t1059-t1060+t1062+t1063-t1064+t741-
           t743-t744+t746+t748-t749;
   t1066 = t1065*t672;
   t1069 = EIz*t828;
   t1070 = t977-t763+t886-t978+t22;
   t1071 = t1070*t649;
   t1072 = t1040-t1041+t667-t522+t668;
   t1073 = t1072*t662;
   t1074 = t832*t91;
   t1076 = t642*t325;
   t1077 = t963-t964+t967*t15-t1076;
   t1078 = t1077*t672;
   t1079 = t1071+t1073-t1074+t1078;
   t1081 = t1015-t1016+t1018-t794+t797-t799+t1020-
           t1022-t1023+t1025+t1026-t812;
   t1082 = t1081*t649;
   t1083 = t1070*t91;
   t1084 = t1058-t1059-t1060+t1062+t1063-t1064-t741+
           t743+t744-t746-t748+t749;
   t1085 = t1084*t662;
   t1086 = t822*t91;
   t1088 = t991*t304;
   t1089 = t696*t325;
   t1090 = t325*t4;
   t1092 = 2.0*t701*t1090;
   t1094 = (-t983-t985+t987+t989*t15-
           t1088-t1089+t1092)*t672;
   t1097 = t861+t908+t910+t939+t941+t961+t962*
           t981+t685*(t995-t996+t1013-t1014+t1028)+
           t1031*t1044+t778*(t1047-t1048+t1056-
           t1057+t1066)+t1069*t1079+t838*(t1082-
           t1083+t1085-t1086+t1094);
   K[3][4] = L*t1097;
   t1100 = t860*t311;
   t1101 = t102*t14;
   t1113 = t116*t304/4.0;
   t1115 = 3.0/4.0*t503*t304;
   t1117 = 3.0/4.0*t508*t304;
   t1118 = -t1113-t1115+t1117+t903-t904+t596-
           t599-t601+t604+t606-t607;
   t1121 = t479*t848/4.0;
   t1124 = 5.0/2.0*t485*t848*t95;
   t1126 = t94*t15*t95;
   t1127 = t1126/2.0;
   t1129 = t101*t848/4.0;
   t1131 = 4.0*t495*t848;
   t1132 = t101*t15;
   t1133 = t875+t877-t879+t1121-t1124+t1127-
           t1129+t1131-t1132+t893-t579+t894;
   t1136 = t156*((-t479*t1101/4.0+5.0/2.0*t865*
           t178+t101*t1101/4.0-4.0*t495*t1101-
           t548+t658)*t110+t1118*t130+t1133*t147);
   t1137 = t909*t332;
   t1138 = t1113+t1115-t1117-t903+t904+t596-
           t599-t601+t604+t606-t607;
   t1140 = t164*t14;
   t1154 = t479*t1090/4.0;
   t1156 = 5.0/2.0*t597*t305;
   t1158 = t101*t1090/4.0;
   t1160 = 4.0*t495*t1090;
   t1161 = t522/2.0;
   t1162 = -t502-t505+t507+t510-t512+t1154-
           t1156-t1158+t1160+t1161-t528;
   t1165 = t192*(t1138*t110+(-t479*t1140/4.0+5.0/2.0*
           t916*t178-2.0*t548+t101*t1140/4.0-4.0*t495*
           t1140+4.0*t554)*t130+t1162*t147);
   t1166 = t940*t348;
   t1167 = -t875-t877+t879+t1121-t1124+t1127-
           t1129+t1131-t1132+t893-t579+t894;
   t1169 = t502+t505-t507-t510+t512+t1154-
           t1156-t1158+t1160+t1161-t528;
   t1171 = t202*t14;
   t1183 = t214*(t1167*t110+t1169*t130+
           (-t479*t1171/4.0+5.0/2.0*t949*t178-t548+
           t101*t1171/4.0-4.0*t495*t1171+t658)*t147);
   t1184 = t506/2.0;
   t1185 = t511/2.0;
   t1187 = -t506/2.0+t511/2.0;
   t1188 = t1187*t21;
   t1190 = t1184-t1185+t1188*t5-t665;
   t1191 = t1190*t649;
   t1192 = t1188*t36;
   t1193 = 2.0*t1132;
   t1194 = t1192-t763+t1126-t1193+t22;
   t1195 = t1194*t662;
   t1196 = t1188*t44;
   t1197 = t1196-t825+t972-t603+t973;
   t1198 = t1197*t672;
   t1199 = t682*t91;
   t1200 = t1191+t1195+t1198-t1199;
   t1202 = t117/4.0;
   t1204 = 3.0/4.0*t503*t44;
   t1206 = 3.0/4.0*t508*t44;
   t1208 = (t1202+t1204-t1206)*t21;
   t1210 = t1187*t100;
   t1212 = t1188*t4;
   t1214 = (-t1202-t1204+t1206+t1208*t5-t1210*
           t513+t1212-t733+t737-t738)*t649;
   t1215 = t1190*t91;
   t1216 = t1208*t36;
   t1217 = t1210*t119;
   t1218 = t1188*t9;
   t1219 = t1218/2.0;
   t1221 = t479*t304/4.0;
   t1223 = 5.0/2.0*t485*t305;
   t1224 = t307/4.0;
   t1226 = 4.0*t495*t304;
   t1227 = t1216-t1217+t1219-t794+t797-t799+t1221-
           t1223-t1224+t1226+t1026-t812;
   t1228 = t1227*t662;
   t1229 = t1208*t44;
   t1230 = t1210*t136;
   t1231 = t1188*t14;
   t1232 = t1231/2.0;
   t1233 = t1229-t1230+t1232-t846+t850-t852+t1001-
           t1002+t1003-t1005+t1008+t1009-t1011;
   t1234 = t1233*t672;
   t1235 = t669*t91;
   t1238 = t1192-t763-t1126+t1193-t22;
   t1239 = t1238*t649;
   t1241 = t1184-t1185+t1188*t10-t1041;
   t1242 = t1241*t662;
   t1243 = t1188*t67;
   t1244 = t1243-t1076+t657+t548-t658;
   t1245 = t1244*t672;
   t1246 = t775*t91;
   t1247 = t1239+t1242+t1245-t1246;
   t1249 = t1216-t1217+t1219-t794+t797-t799-t1221+
           t1223+t1224-t1226-t1026+t812;
   t1250 = t1249*t649;
   t1251 = t1238*t91;
   t1253 = t1210*t232;
   t1255 = (-t1202-t1204+t1206+t1208*t10-
           t1253-t1060+t1062)*t662;
   t1256 = t1208*t67;
   t1257 = t1210*t180;
   t1258 = t715/2.0;
   t1259 = t718/2.0;
   t1260 = t1256-t1257-t1089+t1092+t1258-t1259+
           t721-t723-t724+t726+t728-t729;
   t1261 = t1260*t672;
   t1262 = t765*t91;
   t1265 = t1196-t825+t972+t603-t973;
   t1266 = t1265*t649;
   t1267 = t1243-t1076+t657-t548+t658;
   t1268 = t1267*t662;
   t1270 = t15*t14;
   t1272 = t1184-t1185+t1188*t15-t642*t1270+t666;
   t1273 = t1272*t672;
   t1274 = t835*t91;
   t1275 = t1266+t1268+t1273-t1274;
   t1277 = t1229-t1230+t1232-t846+t850-t852+t1001-
           t1002+t1003+t1005-t1008-t1009+t1011;
   t1278 = t1277*t649;
   t1279 = t1265*t91;
   t1280 = t1256-t1257-t1089+t1092+t1258-t1259-
           t721+t723+t724-t726-t728+t729;
   t1281 = t1280*t662;
   t1283 = t1210*t304;
   t1285 = t1270*t4;
   t1287 = 2.0*t701*t1285;
   t1289 = (-t1202-t1204+t1206+t1208*t15-t1283-
           t696*t1270+t1287+t735-t738)*t672;
   t1290 = t826*t91;
   t1293 = t1100+t1136+t1137+t1165+t1166+t1183+
           t962*t1200+t685*(t1214-t1215+t1228+t1234-
           t1235)+t1031*t1247+t778*(t1250-t1251+
           t1255+t1261-t1262)+t1069*t1275+t838*(t1278-
           t1279+t1281+t1289-t1290);
   K[3][5] = L*t1293;
   K[3][6] = -L*t218/2.0;
   K[3][7] = -L*t393/2.0;
   K[3][8] = -L*t448/2.0;
   t1302 = t650+t653+t663+t673;
   t1306 = t755+t757+t761+t766;
   t1310 = t819+t821+t823+t827;
   t1314 = t477+t559+t561+t611+t613+t634+t962*
           t1302+t685*(t709+t731+t751)+t1031*t1306+
           t778*(t780+t792+t815)+t1069*t1310+t838*
           (t840+t844+t854);
   K[3][9] = L*t1314;
   t1317 = t970+t975+t976+t980;
   t1321 = t1033+t1038+t1039+t1043;
   t1325 = t1071+t1073+t1074+t1078;
   t1329 = t861+t908+t910+t939+t941+t961+t962*
           t1317+t685*(t995-t996+t1013+t1014+t1028)+
           t1031*t1321+t778*(t1047-t1048+t1056+t1057+
           t1066)+t1069*t1325+t838*(t1082-t1083+
           t1085+t1086+t1094);
   K[3][10] = L*t1329;
   t1332 = t1191+t1195+t1198+t1199;
   t1336 = t1239+t1242+t1245+t1246;
   t1340 = t1266+t1268+t1273+t1274;
   t1344 = t1100+t1136+t1137+t1165+t1166+t1183+t962*
           t1332+t685*(t1214-t1215+t1228+t1234+t1235)+
           t1031*t1336+t778*(t1250-t1251+t1255+t1261+
           t1262)+t1069*t1340+t838*(t1278-t1279+t1281+
           t1289+t1290);
   K[3][11] = L*t1344;
   K[4][0] = K[0][10];
   K[4][1] = K[1][10];
   K[4][2] = K[2][10];
   K[4][3] = K[3][4];
   t1347 = t245*t245;
   t1348 = EA*t1347;
   t1349 = t10*t23;
   t1355 = 2.0*t886;
   t1360 = 4.0*t892;
   t1364 = t116*t259/4.0;
   t1366 = 3.0/4.0*t503*t259;
   t1368 = 3.0/4.0*t508*t259;
   t1370 = t479*t1052/4.0;
   t1373 = 5.0/2.0*t485*t1052*t95;
   t1375 = t101*t1052/4.0;
   t1377 = 4.0*t495*t1052;
   t1378 = -t1364-t1366+t507+t1368-t512+t1370-
           t1373+t523-t1375+t1377-t529;
   t1381 = t116*t1035/4.0;
   t1383 = 3.0/4.0*t503*t1035;
   t1384 = 3.0/4.0*t536;
   t1386 = 3.0/4.0*t508*t1035;
   t1387 = 3.0/4.0*t540;
   t1388 = t1381+t1383-t1384-t1386+t1387+t927-
           t930+t935-t932+t934-t554;
   t1391 = t156*((-t479*t1349/4.0+5.0/2.0*t485*t96*
           t10-t1355+t101*t1349/4.0-t492-4.0*t495*
           t1349+t1360+t498-t22)*t110+t1378*
           t130+t1388*t147);
   t1392 = t266*t266;
   t1393 = GAy*t1392;
   t1394 = t1364+t1366-t507-t1368+t512+t1370-
           t1373+t523-t1375+t1377-t529;
   t1396 = t10*t59;
   t1408 = t1035*t14;
   t1410 = t479*t1408/4.0;
   t1413 = 5.0/2.0*t485*t1408*t95;
   t1414 = 3.0/2.0*t603;
   t1416 = t101*t1408/4.0;
   t1418 = 4.0*t495*t1408;
   t1419 = 3.0*t607;
   t1420 = -t898-t900+t903+t902-t904+t1410-
           t1413+t1414-t1416+t1418-t1419;
   t1423 = t192*(t1394*t110+(-t479*t1396/4.0+5.0/2.0*
           t485*t161*t10+t101*t1396/4.0-t576-4.0*t495*
           t1396+t581)*t130+t1420*t147);
   t1424 = t282*t282;
   t1425 = GAz*t1424;
   t1426 = -t1381-t1383+t1384+t1386-t1387+t927-
           t930+t935-t932+t934-t554;
   t1428 = t898+t900-t903-t902+t904+t1410-
           t1413+t1414-t1416+t1418-t1419;
   t1430 = t10*t82;
   t1443 = t214*(t1426*t110+t1428*t130+(-t479*
           t1430/4.0+5.0/2.0*t485*t199*t10-t1355+t101*
           t1430/4.0-t627-4.0*t495*t1430+
           t1360+t630-t22)*t147);
   t1444 = t981*t981;
   t1446 = t241/4.0;
   t1448 = 3.0/4.0*t503*t10;
   t1450 = 3.0/4.0*t508*t10;
   t1452 = (t1446+t1448-t689-t1450+t692)*t21;
   t1456 = (-t1446-t1448+t689+t1450-t692+t1452*t5-
           2.0*t998+t788-t1002)*t649;
   t1457 = t1452*t36;
   t1458 = 2.0*t1050;
   t1460 = t479*t259/4.0;
   t1462 = 5.0/2.0*t485*t260;
   t1463 = t262/4.0;
   t1465 = 4.0*t495*t259;
   t1466 = t1457-t1458+t993+t1054-t719+t1460-
           t1462-t1463+t726+t1465-t729;
   t1467 = t1466*t662;
   t1469 = 2.0*t974*t91;
   t1470 = t1452*t44;
   t1471 = 2.0*t1059;
   t1473 = t479*t1035/4.0;
   t1476 = 5.0/2.0*t485*t1035*t95;
   t1477 = 3.0/2.0*t745;
   t1479 = t101*t1035/4.0;
   t1481 = 4.0*t495*t1035;
   t1482 = 3.0*t749;
   t1484 = (t1470-t1471+t1062-t1064-t1473+t1476-
           t1477+t1479-t1481+t1482)*t672;
   t1487 = t1044*t1044;
   t1489 = t1457-t1458+t993+t1054-t719-t1460+
           t1462+t1463-t726-t1465+t729;
   t1490 = t1489*t649;
   t1495 = t10*t10;
   t1499 = -t1446-t1448+t689+t1450-t692+t1452*t10-
           2.0*t991*t1035+2.0*t999+2.0*t701*t1495-
           5.0/2.0*t789+t707;
   t1500 = t1499*t662;
   t1502 = 2.0*t1037*t91;
   t1503 = t1452*t67;
   t1504 = t991*t259;
   t1505 = 2.0*t1504;
   t1507 = 2.0*t701*t1408;
   t1508 = 3.0/2.0*t798;
   t1509 = t1503-t1505+t1017+t1507-t1508+t1020-
           t1022+t1026-t1023+t1025-t812;
   t1510 = t1509*t672;
   t1513 = t1079*t1079;
   t1516 = (t1470-t1471+t1062-t1064+t1473-t1476+
           t1477-t1479+t1481-t1482)*t649;
   t1517 = t1503-t1505+t1017+t1507-t1508-
           t1020+t1022-t1026+t1023-t1025+t812;
   t1518 = t1517*t662;
   t1520 = 2.0*t1072*t91;
   t1522 = t991*t325;
   t1524 = t15*t10;
   t1526 = 2.0*t701*t1524;
   t1528 = (-t1446-t1448+t689+t1450-t692+t1452*
           t15-2.0*t1522+t1526-t852)*t672;
   t1531 = t1348+t1391+t1393+t1423+t1425+t1443+
           EIx*t1444+t685*(t1456+t1467-t1469+t1484)+
           EIy*t1487+t778*(t1490+t1500-t1502+t1510)+
           EIz*t1513+t838*(t1516+t1518-t1520+t1528);
   K[4][4] = L*t1531;
   t1535 = EA*t245*t311;
   t1536 = t222*t14;
   t1550 = t116*t325/4.0;
   t1552 = 3.0/4.0*t503*t325;
   t1554 = 3.0/4.0*t508*t325;
   t1555 = -t1550-t1552+t1554+t537-t541+t927-
           t930+t935-t932+t934-t554;
   t1557 = t1364+t1366-t507-t1368+t512+t1154-
           t1156-t1158+t1160+t1161-t528;
   t1560 = t156*((-t479*t1536/4.0+5.0/2.0*t865*t602-
           2.0*t603+t101*t1536/4.0-4.0*t495*t1536+4.0*t607)*
           t110+t1555*t130+t1557*t147);
   t1562 = GAy*t266*t332;
   t1563 = t1550+t1552-t1554-t537+t541+t927-
           t930+t935-t932+t934-t554;
   t1565 = t254*t14;
   t1577 = t479*t1524/4.0;
   t1580 = 5.0/2.0*t485*t1524*t95;
   t1582 = t101*t1524/4.0;
   t1584 = 4.0*t495*t1524;
   t1585 = -t875-t877+t879+t1577-t1580+t1127-
           t1582+t1584-t1132+t887-t892+t894;
   t1588 = t192*(t1563*t110+(-t479*t1565/4.0+5.0/2.0*t916*
           t602+t101*t1565/4.0-4.0*t495*t1565-t603+t973)*
           t130+t1585*t147);
   t1590 = GAz*t282*t348;
   t1591 = -t1364-t1366+t507+t1368-t512+t1154-
           t1156-t1158+t1160+t1161-t528;
   t1593 = t875+t877-t879+t1577-t1580+t1127-
           t1582+t1584-t1132+t887-t892+t894;
   t1595 = t277*t14;
   t1607 = t214*(t1591*t110+t1593*t130+(-t479*
           t1595/4.0+5.0/2.0*t949*t602-t603+t101*t1595/4.0-
           4.0*t495*t1595+t973)*t147);
   t1608 = EIx*t981;
   t1610 = t230/4.0;
   t1612 = 3.0/4.0*t503*t67;
   t1614 = 3.0/4.0*t508*t67;
   t1616 = (t1610+t1612-t1614)*t21;
   t1619 = (-t1610-t1612+t1614+t1616*t5-
           t1217-t1016+t797)*t649;
   t1620 = t1616*t36;
   t1621 = t1212/2.0;
   t1623 = t479*t325/4.0;
   t1625 = 5.0/2.0*t485*t326;
   t1626 = t328/4.0;
   t1628 = 4.0*t495*t325;
   t1629 = t1620-t1253+t1621-t1059+t1062-t1064+t1623-
           t1625-t1626+t1628+t746-t749;
   t1630 = t1629*t662;
   t1631 = t1194*t91;
   t1632 = t1616*t44;
   t1633 = t993/2.0;
   t1634 = t1632-t1257-t1088+t1092+t1633-t1259-t1460+
           t1462+t1463-t726-t1465+t729;
   t1635 = t1634*t672;
   t1636 = t979*t91;
   t1639 = EIy*t1044;
   t1641 = t1620-t1253+t1621-t1059+t1062-t1064-t1623+
           t1625+t1626-t1628-t746+t749;
   t1642 = t1641*t649;
   t1646 = (-t1610-t1612+t1614+t1616*t10-t1210*t1035+
           t1218-t1504+t1507-t798)*t662;
   t1647 = t1241*t91;
   t1648 = t1616*t67;
   t1649 = t1210*t259;
   t1650 = t1648-t1649+t1232-t1522+t1526-t852+t1000-
           t790+t1003+t1005-t1008-t1009+t1011;
   t1651 = t1650*t672;
   t1652 = t1042*t91;
   t1655 = EIz*t1079;
   t1657 = t1632-t1257-t1088+t1092+t1633-t1259+t1460-
           t1462-t1463+t726+t1465-t729;
   t1658 = t1657*t649;
   t1659 = t1648-t1649+t1232-t1522+t1526-t852+t1000-
           t790+t1003-t1005+t1008+t1009-t1011;
   t1660 = t1659*t662;
   t1661 = t1267*t91;
   t1663 = t1210*t325;
   t1665 = t1270*t9;
   t1667 = 2.0*t701*t1665;
   t1669 = (-t1610-t1612+t1614+t1616*t15-t1663-
           t991*t1270+t1667+t1017-t798)*t672;
   t1670 = t1077*t91;
   t1673 = t1535+t1560+t1562+t1588+t1590+t1607+t1608*
           t1200+t685*(t1619+t1630-t1631+t1635-t1636)+
           t1639*t1247+t778*(t1642+t1646-t1647+t1651-
           t1652)+t1655*t1275+t838*(t1658+t1660-
           t1661+t1669-t1670);
   K[4][5] = L*t1673;
   K[4][6] = -L*t288/2.0;
   K[4][7] = -L*t405/2.0;
   K[4][8] = -L*t460/2.0;
   t1691 = t861+t908+t910+t939+t941+t961+t1608*
           t1302+t685*(t995+t996+t1013-t1014+t1028)+t1639*
           t1306+t778*(t1047+t1048+t1056-t1057+t1066)+
           t1655*t1310+t838*(t1082+t1083+t1085-t1086+t1094);
   K[4][9] = L*t1691;
   t1703 = t1348+t1391+t1393+t1423+t1425+t1443+t1608*
           t1317+t685*(t1456+t1467+t1484)+t1639*t1321+t778*
           (t1490+t1500+t1510)+t1655*t1325+t838*(t1516+
           t1518+t1528);
   K[4][10] = L*t1703;
   t1715 = t1535+t1560+t1562+t1588+t1590+t1607+t1608*
           t1332+t685*(t1619+t1630-t1631+t1635+t1636)+
           t1639*t1336+t778*(t1642+t1646-t1647+t1651+
           t1652)+t1655*t1340+t838*(t1658+t1660-t1661+
           t1669+t1670);
   K[4][11] = L*t1715;
   K[5][0] = K[0][11];
   K[5][1] = K[1][11];
   K[5][2] = K[2][11];
   K[5][3] = K[3][5];
   K[5][4] = K[4][5];
   t1718 = t311*t311;
   t1719 = EA*t1718;
   t1720 = t15*t23;
   t1726 = 2.0*t1126;
   t1731 = 4.0*t1132;
   t1735 = t116*t1270/4.0;
   t1737 = 3.0/4.0*t503*t1270;
   t1738 = 3.0/4.0*t506;
   t1740 = 3.0/4.0*t508*t1270;
   t1741 = 3.0/4.0*t511;
   t1742 = -t1735-t1737+t1738+t1740-t1741+t1154-
           t1156+t1161-t1158+t1160-t528;
   t1745 = t479*t1285/4.0;
   t1748 = 5.0/2.0*t485*t1285*t95;
   t1750 = t101*t1285/4.0;
   t1752 = 4.0*t495*t1285;
   t1753 = t1550+t1552-t537-t1554+t541+t1745-
           t1748+t549-t1750+t1752-t555;
   t1756 = t156*((-t479*t1720/4.0+5.0/2.0*t485*t96*t15-
           t1726+t101*t1720/4.0-t492-4.0*t495*t1720+t1731+
           t498-t22)*t110+t1742*t130+t1753*t147);
   t1757 = t332*t332;
   t1758 = GAy*t1757;
   t1759 = t1735+t1737-t1738-t1740+t1741+t1154-
           t1156+t1161-t1158+t1160-t528;
   t1761 = t15*t59;
   t1774 = t479*t1665/4.0;
   t1777 = 5.0/2.0*t485*t1665*t95;
   t1779 = t101*t1665/4.0;
   t1781 = 4.0*t495*t1665;
   t1782 = -t1113-t1115+t903+t1117-t904+t1774-
           t1777+t1414-t1779+t1781-t1419;
   t1785 = t192*(t1759*t110+(-t479*t1761/4.0+5.0/2.0*t485*
           t161*t15-t1726+t101*t1761/4.0-t576-4.0*t495*
           t1761+t1731+t581-t22)*t130+t1782*t147);
   t1786 = t348*t348;
   t1787 = GAz*t1786;
   t1788 = -t1550-t1552+t537+t1554-t541+t1745-
           t1748+t549-t1750+t1752-t555;
   t1790 = t1113+t1115-t903-t1117+t904+t1774-
           t1777+t1414-t1779+t1781-t1419;
   t1792 = t15*t82;
   t1805 = t214*(t1788*t110+t1790*t130+
           (-t479*t1792/4.0+5.0/2.0*t485*t199*t15+t101*
           t1792/4.0-t627-4.0*t495*t1792+t630)*t147);
   t1806 = t1200*t1200;
   t1808 = t300/4.0;
   t1810 = 3.0/4.0*t503*t15;
   t1812 = 3.0/4.0*t508*t15;
   t1814 = (t1808+t1810-t689-t1812+t692)*t21;
   t1818 = (-t1808-t1810+t689+t1812-t692+t1814*t5-
           2.0*t1230+t850-t1002)*t649;
   t1819 = t1814*t36;
   t1820 = 2.0*t1257;
   t1822 = t479*t1270/4.0;
   t1825 = 5.0/2.0*t485*t1270*t95;
   t1826 = 3.0/2.0*t725;
   t1828 = t101*t1270/4.0;
   t1830 = 4.0*t495*t1270;
   t1831 = 3.0*t729;
   t1833 = (t1819-t1820+t1092-t1259+t1822-t1825+t1826-
           t1828+t1830-t1831)*t662;
   t1834 = t1814*t44;
   t1835 = 2.0*t1283;
   t1836 = t1834-t1835+t1212+t1287-t739-t1623+t1625-
           t746+t1626-t1628+t749;
   t1837 = t1836*t672;
   t1839 = 2.0*t1197*t91;
   t1842 = t1247*t1247;
   t1845 = (t1819-t1820+t1092-t1259-t1822+t1825-
           t1826+t1828-t1830+t1831)*t649;
   t1849 = (-t1808-t1810+t689+t1812-t692+t1814*t10-
           2.0*t1649+t1526-t790)*t662;
   t1850 = t1814*t67;
   t1851 = 2.0*t1663;
   t1852 = t1850-t1851+t1218+t1667-t1508+t1221-
           t1223+t1026-t1224+t1226-t812;
   t1853 = t1852*t672;
   t1855 = 2.0*t1244*t91;
   t1858 = t1275*t1275;
   t1860 = t1834-t1835+t1212+t1287-t739+t1623-
           t1625+t746-t1626+t1628-t749;
   t1861 = t1860*t649;
   t1862 = t1850-t1851+t1218+t1667-t1508-t1221+t1223-
           t1026+t1224-t1226+t812;
   t1863 = t1862*t662;
   t1868 = t15*t15;
   t1872 = -t1808-t1810+t689+t1812-t692+t1814*t15-
           2.0*t1210*t1270+2.0*t1231+2.0*t701*
           t1868-5.0/2.0*t851+t707;
   t1873 = t1872*t672;
   t1875 = 2.0*t1272*t91;
   t1878 = t1719+t1756+t1758+t1785+t1787+t1805+
           EIx*t1806+t685*(t1818+t1833+t1837-t1839)+
           EIy*t1842+t778*(t1845+t1849+t1853-t1855)+
           EIz*t1858+t838*(t1861+t1863+t1873-t1875);
   K[5][5] = L*t1878;
   K[5][6] = -L*t354/2.0;
   K[5][7] = -L*t417/2.0;
   K[5][8] = -L*t472/2.0;
   t1887 = EIx*t1200;
   t1891 = EIy*t1247;
   t1895 = EIz*t1275;
   t1899 = t1100+t1136+t1137+t1165+t1166+t1183+t1887*
           t1302+t685*(t1214+t1215+t1228+t1234-t1235)+
           t1891*t1306+t778*(t1250+t1251+t1255+t1261-
           t1262)+t1895*t1310+t838*(t1278+t1279+t1281+
           t1289-t1290);
   K[5][9] = L*t1899;
   t1911 = t1535+t1560+t1562+t1588+t1590+t1607+t1887*
           t1317+t685*(t1619+t1630+t1631+t1635-t1636)+
           t1891*t1321+t778*(t1642+t1646+t1647+t1651-
           t1652)+t1895*t1325+t838*(t1658+t1660+t1661+
           t1669-t1670);
   K[5][10] = L*t1911;
   t1923 = t1719+t1756+t1758+t1785+t1787+t1805+t1887*
           t1332+t685*(t1818+t1833+t1837)+t1891*t1336+
           t778*(t1845+t1849+t1853)+t1895*t1340+t838*
           (t1861+t1863+t1873);
   K[5][11] = L*t1923;
   K[6][0] = K[0][6];
   K[6][1] = K[1][6];
   K[6][2] = K[2][6];
   K[6][3] = K[3][6];
   K[6][4] = K[4][6];
   K[6][5] = K[5][6];
   K[6][6] = K[0][0];
   K[6][7] = K[0][1];
   K[6][8] = K[0][2];
   K[6][9] = K[3][6];
   K[6][10] = K[4][6];
   K[6][11] = K[5][6];
   K[7][0] = K[1][6];
   K[7][1] = K[1][7];
   K[7][2] = K[2][7];
   K[7][3] = K[3][7];
   K[7][4] = K[4][7];
   K[7][5] = K[5][7];
   K[7][6] = K[6][7];
   K[7][7] = K[1][1];
   K[7][8] = K[1][2];
   K[7][9] = K[3][7];
   K[7][10] = K[4][7];
   K[7][11] = K[5][7];
   K[8][0] = K[2][6];
   K[8][1] = K[2][7];
   K[8][2] = K[2][8];
   K[8][3] = K[3][8];
   K[8][4] = K[4][8];
   K[8][5] = K[5][8];
   K[8][6] = K[6][8];
   K[8][7] = K[7][8];
   K[8][8] = K[2][2];
   K[8][9] = K[3][8];
   K[8][10] = K[4][8];
   K[8][11] = K[5][8];
   K[9][0] = K[3][0];
   K[9][1] = K[3][1];
   K[9][2] = K[3][2];
   K[9][3] = K[3][9];
   K[9][4] = K[4][9];
   K[9][5] = K[5][9];
   K[9][6] = K[6][9];
   K[9][7] = K[7][9];
   K[9][8] = K[8][9];
   t1926 = t1302*t1302;
   t1930 = t1306*t1306;
   t1934 = t1310*t1310;
   t1938 = t477+t559+t561+t611+t613+t634+
           EIx*t1926+t685*(t709+t711+t731+t751)+
           EIy*t1930+t778*(t780+t782+t792+t815)+
           EIz*t1934+t838*(t840+t842+t844+t854);
   K[9][9] = L*t1938;
   t1941 = EIx*t1302;
   t1945 = EIy*t1306;
   t1949 = EIz*t1310;
   t1953 = t861+t908+t910+t939+t941+t961+t1941*
           t1317+t685*(t995+t996+t1013+t1014+
           t1028)+t1945*t1321+t778*(t1047+t1048+
           t1056+t1057+t1066)+t1949*t1325+t838*
           (t1082+t1083+t1085+t1086+t1094);
   K[9][10] = L*t1953;
   t1965 = t1100+t1136+t1137+t1165+t1166+t1183+
           t1941*t1332+t685*(t1214+t1215+t1228+
           t1234+t1235)+t1945*t1336+t778*
           (t1250+t1251+t1255+t1261+t1262)+
           t1949*t1340+t838*(t1278+t1279+t1281+
           t1289+t1290);
   K[9][11] = L*t1965;
   K[10][0] = K[4][0];
   K[10][1] = K[4][1];
   K[10][2] = K[4][2];
   K[10][3] = K[3][10];
   K[10][4] = K[4][10];
   K[10][5] = K[5][10];
   K[10][6] = K[6][10];
   K[10][7] = K[7][10];
   K[10][8] = K[8][10];
   K[10][9] = K[9][10];
   t1968 = t1317*t1317;
   t1972 = t1321*t1321;
   t1976 = t1325*t1325;
   t1980 = t1348+t1391+t1393+t1423+t1425+t1443+
           EIx*t1968+t685*(t1456+t1467+t1469+t1484)+
           EIy*t1972+t778*(t1490+t1500+t1502+t1510)+
           EIz*t1976+t838*(t1516+t1518+t1520+t1528);
   K[10][10] = L*t1980;
   t1995 = t1535+t1560+t1562+t1588+t1590+t1607+
           EIx*t1317*t1332+t685*(t1619+t1630+t1631+
           t1635+t1636)+EIy*t1321*t1336+t778*(t1642+
           t1646+t1647+t1651+t1652)+EIz*t1325*t1340+
           t838*(t1658+t1660+t1661+t1669+t1670);
   K[10][11] = L*t1995;
   K[11][0] = K[5][0];
   K[11][1] = K[5][1];
   K[11][2] = K[5][2];
   K[11][3] = K[3][11];
   K[11][4] = K[4][11];
   K[11][5] = K[5][11];
   K[11][6] = K[6][11];
   K[11][7] = K[7][11];
   K[11][8] = K[8][11];
   K[11][9] = K[9][11];
   K[11][10] = K[10][11];
   t1998 = t1332*t1332;
   t2002 = t1336*t1336;
   t2006 = t1340*t1340;
   t2010 = t1719+t1756+t1758+t1785+t1787+t1805+
           EIx*t1998+t685*(t1818+t1833+t1837+t1839)+
           EIy*t2002+t778*(t1845+t1849+t1853+t1855)+
           EIz*t2006+t838*(t1861+t1863+t1873+t1875);
   K[11][11] = L*t2010;

   return;

}  /* End of BeamTL_LocalStiff */

int cModelBeam3D :: Profile ( int *piProfile )
{
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
 return 1;
}

int cModelBeam3D :: SparsityPattern ( int **piSparsity )
{
 if (piSparsity==0) return( 0 );
 // Generate sparticity matrix
 int i,j; 
 for ( i = 0; i < _iNumElms; i++ )
 {
   for ( j = 0; j < 2; j++ )
   {
      for ( int k = 0; k < 6; k++ )
      {
         if ( _paDof[_paInc[i][j]][k] >= 0 )
         {
            for ( int m = 0; m < 2; m++ )
            {
               for ( int n = 0; n < 6; n++ )
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

