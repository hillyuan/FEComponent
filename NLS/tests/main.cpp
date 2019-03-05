// --------------------------------------------------------------------------
//
// Main.cpp - Main program for testing the NLS library.
//
// --------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//NLS library
#include <nls.h>

#include "modbeam2.h"
#include "modbeam3.h"
#include "modfunc_1d.h"
#include "modfunc_2d.h"
#include "modst.h"
#include "modstp.h"
#include "utl.h"

// --------------------------------------------------------------------------
// Main Program:

int main(int argc, char* argv[])
{
 char      sMod[80];
 int       iAlg;
 int       iLinSys;
 int	     iLinSysMaxIter;
 double    iLinSysTol;
 cControl  *pcCtrl;
 sControl  sCtrl;
 cModel    *pcModel;
 cLinSys   *pcLinSys;

 // check input arguments
 if (argc < 3) {
  fprintf(stderr, "Usage:\n   tn1s model_file algorithm_file\n");
  exit(-1);
 }

 // get space for filenames
 char *pcModFileName = new char[strlen(argv[1])+1]; 
 char *pcAlgFileName = new char[strlen(argv[2])+1];
 
 // get arguments from command line
 strcpy(pcModFileName,argv[1]); //model file name, i.e. model_lee.inp
 strcpy(pcAlgFileName,argv[2]); //algorithm file, i.e. alg_lcm.inp


  //char *pcModFileName = "D:\\Sofie's Documents\\UIUC\\Research\\NLS\\NLS Code\\bin\\vc7\\model_unidimensional.inp"; 
  //char *pcAlgFileName = "D:\\Sofie's Documents\\UIUC\\Research\\NLS\\NLS Code\\bin\\vc7\\alg_gdcm.inp";


 /*Try to open model file*/
 FILE *fpi = fopen( pcModFileName, "r" );
 if( !fpi )
 {
  printf( "\n\nModel file %s not found...\n\n", pcModFileName );
  exit( -1 );
 }

 /*reads model type from first line of model input file to sMod*/
 read_string( fpi, sMod );

 /*close model file*/
 fclose(fpi);

 /*Create a new model based on sMod*/

 if(strcmp(sMod, "function" ) == 0)
 {
   pcModel = new cModelFunction( pcModFileName );
 }
 else if(strcmp(sMod, "function_mt" ) == 0)
 {
   pcModel = new cModelFunction_MT( pcModFileName );
 }
 else if(strcmp(sMod, "function_nobre" ) == 0)
 {
   pcModel = new cModelFunction_Nobre( pcModFileName );
 }
 else if(strcmp(sMod, "space_truss") == 0) 
 { 
   printf("\n\ncreating a cModelSpaceTruss ( %s )\n", pcModFileName );
   pcModel = new cModelSpaceTruss( pcModFileName );
 }
 else if( strcmp ( sMod, "plane_frame" ) == 0 ) 
 {
   printf("\n\ncreating a cModelBeam2D ( %s )\n", pcModFileName );
   pcModel = new cModelBeam2D( pcModFileName );
 } 
 else if( strcmp ( sMod, "space_frame" ) == 0 ) 
 {
   printf("\n\ncreating a cModelBeam3D ( %s )\n", pcModFileName );
   pcModel = new cModelBeam3D( pcModFileName );
 } 
 else 
 {
   printf( "\n\nModel not available...\n\n" );
   exit( -1 );
 }

 /*Initialize model*/
 pcModel->Init( );

 /* Open algorithm file */
 fpi = fopen( pcAlgFileName, "r" );
 if( !fpi )
 {
  printf( "\n\nAlgoritm file %s not found...\n\n", pcAlgFileName );
  exit( -1 );
 }

 /* Read the linear solver */
 fscanf( fpi, "%d", &iLinSys);

 if( (iLinSys == 1) ||
     (iLinSys == 2) )
  fscanf( fpi, "%d%lf", &iLinSysMaxIter, &iLinSysTol );

 /* Creates a new LinSys (LinearSystems) */
 switch( iLinSys )
 {
  case 0:
        pcLinSys = new cCroutProfile( );
  break;

  case 1:
        pcLinSys = new cPCGProfile( iLinSysMaxIter, iLinSysTol );
  break;

  case 2:
//        pcLinSys = new cPCGCSRSym( iLinSysMaxIter, iLinSysTol );
  break;

  default:
        printf( "\n\nLinear solver not available...\n\n" );
        exit( -1 );
  break;
 }

 /* Read algorithm type */
 fscanf( fpi, "%d%d%lf", &iAlg, &sCtrl.UpdateType, &sCtrl.CtrlFactor );

 /* Read appropriate parameters based on algorithm type */
 if( iAlg == 1 ) fscanf( fpi, "%d%d", &sCtrl.CtrlEq, &sCtrl.CtrlType );

 if( iAlg == 2 ) fscanf( fpi, "%d",   &sCtrl.CtrlType );

 if( iAlg == 5 ) fscanf( fpi, "%lf",  &sCtrl.CtrlIniFactor );

 if( iAlg == 6 ) fscanf( fpi, "%lf",  &sCtrl.CtrlIniFactor );

 if( iAlg == 7 ) fscanf( fpi, "%d%d", &sCtrl.CtrlEq, &sCtrl.CtrlType );

 if( iAlg == 8 ) fscanf( fpi, "%lf",  &sCtrl.CtrlIniFactor );

 fscanf( fpi, "%d%d%lf", &sCtrl.NumMaxStep, &sCtrl.NumMaxIte, &sCtrl.Tol );

 /* Close algorithm file */
 fclose( fpi );

 /*create new Control*/
 switch( iAlg )
 {
  case 0:
        pcCtrl = new cNewtonRaphson( pcModel, &sCtrl, pcLinSys );
  break;

  case 1:
        pcCtrl = new cDisplacementControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 2:
        pcCtrl = new cArcLengthControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 3:
        pcCtrl = new cWorkControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 4:
        pcCtrl = new cGenDisplacementControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 5:
        pcCtrl = new cOrthResidualControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 6:
        pcCtrl = new cStrainRatioControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 7:
        pcCtrl = new cStrainControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 8:
        pcCtrl = new cOldOrthResidualControl( pcModel, &sCtrl, pcLinSys );
  break;

  case 9:
        pcCtrl = new cModGenDisplacementControl( pcModel, &sCtrl, pcLinSys );
  break;

  default:
        printf( "\n\nAlgorithm not available...\n\n" );
        exit( -1 );
  break;
 }

 /* run */
 pcCtrl->Solver( );

 

 delete pcCtrl;
 delete pcModel;
 return 1;
}
