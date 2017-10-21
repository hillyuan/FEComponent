/*
    Purpose: Build classes of 1D numerical integration and , as an except,
        2D Hammer integration. Those classes are used to generate quadrature
        classes used by the user and as conseqence, transprent to users.
           This is a part of Library of Generic Numerical Components *GNCLib*.
           You need only add your new derived class of QuadratureBase to
        accomplish different quadrature scheme. ID of each quadrature scheme
        is randam number generated automatically. Dangerous in generating same
        ID then?


        Date                    Programmer              Status
        ===========            ============          =============
        Oct. 08,2005             YUAN Xi               Original
*/


#ifndef _strainfunction_h_
#define _strainfunction_h_

#include <stdlib.h>
#include "physicalBase.h"


namespace GNCLib
{
  /*
     Values of dynamic matrix of all quadrature points.
     Those values can be fetched in numerical caculation with no-calling of
     above shape functions to save computation time. It is much a optimized
     way considering the fact that the memeory used by this table is generally
     disregardably small.
     It is a const singleton.
     Do not use InterpolationTable<CInterpolation, CQuadrature> in the furture
     because it generate a static interpolation table that (maybe?) useless
     it you adopt PhysicalTable here.
   */
   template<class CInterpolation, class CQuadrature>
   class PhysicalTable<Strain_Stress, CInterpolation, CQuadrature>
   {
     public:
       static const unsigned int numcomponents=
        CInterpolation::dimension*(CInterpolation::dimension-1)/2+CInterpolation::dimension;
       static const unsigned int numdofs=CInterpolation::numpoints*CInterpolation::dimension;

     private:
       unsigned int numgauss;
       double*** TPhysicalFunction;

     private:
       PhysicalTable();
       ~PhysicalTable();

     public:
       static PhysicalTable& Instance();

       void PhysicalFunction(double[][numcomponents][numdofs]) const
       {return TPhysicalFunction;}
   };


   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> :: PhysicalTable()
   {
     unsigned int i,j,k;

     unsigned int numintps=CInterpolation::numpoints;
     CQuadrature aquadrature;
     CInterpolation ainterpolator;
     numgauss = aquadrature.numQuadraturePoints();

     double shapederivative[numintps][CInterpolation::dimension];

     //-----allocate functional table
     TPhysicalFunction = (double ***)malloc(numgauss * sizeof(double *));
     if(TPhysicalFunction == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TPhysicalFunction[i] = (double **)malloc(numcomponents * sizeof(double));
         if(TPhysicalFunction[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
         for(j = 0; j < numcomponents; j++) {
             TPhysicalFunction[i][j] = (double *)malloc(numdofs * sizeof(double));
             if(TPhysicalFunction[i][j] == NULL) {
                 fprintf(stderr, "out of memory when building interpolation table\n");
                 exit(1);
             }
         }
     }

     //-----Initilize TPhysicalFunction
     for(i=0; i<numgauss; i++)
       for(j=0; j<numcomponents; j++)
         for(k=0; k<numdofs; k++)
           TPhysicalFunction[i][j][k]=0.0;

     unsigned int dimension = CInterpolation::dimension;
     unsigned int nintps = CInterpolation::numpoints;
     double coord[dimension];
     //specify the value of physical table
     for(i=0; i<numgauss; i++) {
       aquadrature.quadraturePoint(i, coord);
       ainterpolator.ShapeDerivative(coord, shapederivative[i]);
       for(j=0; j<dimension; j++)   // normal strain components
         for(k=0; k<nintps; k++)
           TPhysicalFunction[i][j][k*dimension+j]=shapederivative[k][j];
       for(j=dimension; j<numcomponents; j++) //shear strain components
         for(k=0; k<nintps; k++) {
           if(j==dimension) {  //possibly 2-d and 3-d dimensional problem
             TPhysicalFunction[i][j][k*dimension]  =shapederivative[k][1];
             TPhysicalFunction[i][j][k*dimension+1]=shapederivative[k][0];
           }
           if(j==dimension+1) {  //only 3-d dimensional problem
             TPhysicalFunction[i][j][k*dimension+1]=shapederivative[k][2];
             TPhysicalFunction[i][j][k*dimension+2]=shapederivative[k][1];
           }
           if(j==dimension+2) {  //only 3-d dimensional problem
             TPhysicalFunction[i][j][k*dimension]  =shapederivative[k][2];
             TPhysicalFunction[i][j][k*dimension+2]=shapederivative[k][0];
           }
         }
     }

   };

                         /* Destructor defined here is mostly called at the
                            end of program execution. If is needless in fact
                            ha! ha! ha!
                         */
   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> :: ~PhysicalTable()
   {
     unsigned int i,j;

     for(i = 0; i < numgauss; i++) {
         for(j =0; j < numcomponents; j++)
             free(TPhysicalFunction[i][j]);
         free(TPhysicalFunction[i]);
     }
     free(TPhysicalFunction);
   };

   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress, CInterpolation, CQuadrature>&
   PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> :: Instance()
   {
      static PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> physicaltable;
      return physicaltable;
   }


   //  -------------------------------------------------
   //    Axis-symmetric problem
   //  -------------------------------------------------

   template<class CInterpolation, class CQuadrature>
   class PhysicalTable<Strain_Stress_AxisSymmetric, CInterpolation, CQuadrature>
   {
     public:
       static const unsigned int numcomponents=4;
       static const unsigned int numdofs=CInterpolation::numpoints*CInterpolation::dimension;

     private:
       unsigned int numgauss;
       double*** TPhysicalFunction;

     private:
       PhysicalTable();
       ~PhysicalTable();

     public:
       static PhysicalTable& Instance();

       void PhysicalFunction(double[][numcomponents][numdofs]) const
       {return TPhysicalFunction;}
   };


   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_AxisSymmetric, CInterpolation, CQuadrature> :: PhysicalTable()
   {
     unsigned int i,j,k;

     unsigned int numintps=CInterpolation::numpoints;
     CQuadrature aquadrature;
     CInterpolation ainterpolator;
     numgauss = aquadrature.numQuadraturePoints();

     unsigned int dimension = CInterpolation::dimension;
     unsigned int nintps = CInterpolation::numpoints;
     double shapefunction[numintps];
     double shapederivative[numintps][dimension];


     //-----allocate functional table
     TPhysicalFunction = (double ***)malloc(numgauss * sizeof(double *));
     if(TPhysicalFunction == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TPhysicalFunction[i] = (double **)malloc(numcomponents * sizeof(double));
         if(TPhysicalFunction[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
         for(j = 0; j < numcomponents; j++) {
             TPhysicalFunction[i][j] = (double *)malloc(numdofs * sizeof(double));
             if(TPhysicalFunction[i][j] == NULL) {
                 fprintf(stderr, "out of memory when building interpolation table\n");
                 exit(1);
             }
         }
     }

     //-----Initilize TPhysicalFunction
     for(i=0; i<numgauss; i++)
       for(j=0; j<numcomponents; j++)
         for(k=0; k<numdofs; k++)
           TPhysicalFunction[i][j][k]=0.0;


     double coord[dimension];
     //specify the value of physical table
     for(i=0; i<numgauss; i++) {
       aquadrature.quadraturePoint(i, coord);
       ainterpolator.ShapeDerivative(coord, shapederivative[i]);
       ainterpolator.ShapeFunction(coord, shapefunction);
       for(k=0; k<nintps; k++) {
         for(j=0; j<dimension; j++)   // normal strain components
           TPhysicalFunction[i][j][k*dimension+j]=shapederivative[k][j];
         TPhysicalFunction[i][dimension][k*dimension]  =shapederivative[k][1];
         TPhysicalFunction[i][dimension][k*dimension+1]=shapederivative[k][0];
         TPhysicalFunction[i][dimension+1][k*dimension]=shapefunction[k];
       }

     }

   };

                         /* Destructor defined here is mostly called at the
                            end of program execution. If is needless in fact
                            ha! ha! ha!
                         */
   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_AxisSymmetric, CInterpolation, CQuadrature> :: ~PhysicalTable()
   {
     unsigned int i,j;
     for(i = 0; i < numgauss; i++) {
         for(j =0; j < numcomponents; j++)
             free(TPhysicalFunction[i][j]);
         free(TPhysicalFunction[i]);
     }
     free(TPhysicalFunction);
   };

   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_AxisSymmetric, CInterpolation, CQuadrature>&
   PhysicalTable<Strain_Stress_AxisSymmetric, CInterpolation, CQuadrature> :: Instance()
   {
      static PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> physicaltable;
      return physicaltable;
   }



   //  -------------------------------------------------
   //    Membrane problem
   //  -------------------------------------------------

   template<class CInterpolation, class CQuadrature>
   class PhysicalTable<Strain_Stress_Membrane, CInterpolation, CQuadrature>
   {
     public:
       static const unsigned int numcomponents=3;
       static const unsigned int numdofs=CInterpolation::numintps*CInterpolation::dimension;

     private:
       unsigned int numgauss;
       double*** TPhysicalFunction;

     private:
       PhysicalTable();
       ~PhysicalTable();

     public:
       static PhysicalTable& Instance();

       void PhysicalFunction(double[][numcomponents][numdofs]) const
       {return TPhysicalFunction;}
   };


   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Membrane, CInterpolation, CQuadrature> :: PhysicalTable()
   {
     unsigned int i,j,k;

     unsigned int numintps=CInterpolation::numpoints;
     unsigned int dimension = CInterpolation::dimension;
     CQuadrature aquadrature;
     CInterpolation ainterpolator;
     numgauss = aquadrature.numQuadraturePoints();

     double shapederivative[numintps][dimension];

     //-----allocate functional table
     TPhysicalFunction = (double ***)malloc(numgauss * sizeof(double *));
     if(TPhysicalFunction == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TPhysicalFunction[i] = (double **)malloc(numcomponents * sizeof(double));
         if(TPhysicalFunction[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
         for(j = 0; j < numcomponents; j++) {
             TPhysicalFunction[i][j] = (double *)malloc(numdofs * sizeof(double));
             if(TPhysicalFunction[i][j] == NULL) {
                 fprintf(stderr, "out of memory when building interpolation table\n");
                 exit(1);
             }
         }
     }

     //-----Initilize TPhysicalFunction
     for(i=0; i<numgauss; i++)
       for(j=0; j<numcomponents; j++)
         for(k=0; k<numdofs; k++)
           TPhysicalFunction[i][j][k]=0.0;


     double coord[dimension];
     //specify the value of physical table
     for(i=0; i<numgauss; i++) {
       aquadrature.quadraturePoint(i, coord);
       ainterpolator.ShapeDerivative(coord, shapederivative[i]);
       for(k=0; k<numintps; k++) {
         for(j=0; j<2; j++)   // normal strain components
           TPhysicalFunction[i][j][k*2+j]=shapederivative[k][j];
         TPhysicalFunction[i][2][k*dimension]  =shapederivative[k][1];
         TPhysicalFunction[i][2][k*dimension+1]=shapederivative[k][0];
       }

     }

   };

                         /* Destructor defined here is mostly called at the
                            end of program execution. If is needless in fact
                            ha! ha! ha!
                         */
   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Membrane, CInterpolation, CQuadrature> :: ~PhysicalTable()
   {
     unsigned int i,j;

     for(i = 0; i < numgauss; i++) {
         for(j =0; j < numcomponents; j++)
             free(TPhysicalFunction[i][j]);
         free(TPhysicalFunction[i]);
     }
     free(TPhysicalFunction);
   };

   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Membrane, CInterpolation, CQuadrature>&
   PhysicalTable<Strain_Stress_Membrane, CInterpolation, CQuadrature> :: Instance()
   {
      static PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> physicaltable;
      return physicaltable;
   }


  //  -------------------------------------------------
  //    Shell and beam problems (currently shell only)
  //  -------------------------------------------------
   template<class CInterpolation, class CQuadrature>
   class PhysicalTable<Strain_Stress_Shell, CInterpolation, CQuadrature>
   {
     public:
       static const unsigned int numcomponents=5;
       static const unsigned int numdofs=CInterpolation::numintps*6;

     private:
       unsigned int numgauss;
       double*** TPhysicalFunction;

     private:
       PhysicalTable();
       ~PhysicalTable();   //only for shell element

     public:
       static PhysicalTable& Instance();

       void PhysicalFunction(double[][numcomponents][numdofs]) const
       {return TPhysicalFunction;}
   };


   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Shell, CInterpolation, CQuadrature> :: PhysicalTable()
   {
     unsigned int i,j,k;

     unsigned int numintps=CInterpolation::numpoints;
     unsigned int dimension = CInterpolation::dimension;
     CQuadrature aquadrature;
     CInterpolation ainterpolator;
     numgauss = aquadrature.numQuadraturePoints();

     double shapefunction[numintps];
     double shapederivative[numintps][dimension];

     //-----allocate functional table
     TPhysicalFunction = (double ***)malloc(numgauss * sizeof(double *));
     if(TPhysicalFunction == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TPhysicalFunction[i] = (double **)malloc(numcomponents * sizeof(double));
         if(TPhysicalFunction[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
         for(j = 0; j < numcomponents; j++) {
             TPhysicalFunction[i][j] = (double *)malloc(numdofs * sizeof(double));
             if(TPhysicalFunction[i][j] == NULL) {
                 fprintf(stderr, "out of memory when building interpolation table\n");
                 exit(1);
             }
         }
     }

     //-----Initilize TPhysicalFunction
     for(i=0; i<numgauss; i++)
       for(j=0; j<numcomponents; j++)
         for(k=0; k<numdofs; k++)
           TPhysicalFunction[i][j][k]=0.0;

     double coord[dimension];
     //specify the value of physical table
     for(i=0; i<numgauss; i++) {
       aquadrature.quadraturePoint(i, coord);
       ainterpolator.ShapeDerivative(coord, shapederivative[i]);
       ainterpolator.ShapeFunction(coord, shapefunction);

       for(k=0; k<numintps; k++) {
         for(j=0; j<2; j++)   // normal strain components
           TPhysicalFunction[i][j][k*6+j]=shapederivative[k][j];
         TPhysicalFunction[i][2][k*6]  =shapederivative[k][1];
         TPhysicalFunction[i][2][k*6+1]=shapederivative[k][0];
         TPhysicalFunction[i][3][k*6+2]=shapederivative[k][1];
         TPhysicalFunction[i][4][k*6+2]=shapederivative[k][0];

         for(j=3; j<5; j++)   //
           TPhysicalFunction[i][j][k*6+j]=coord[2]*shapederivative[k][j-3];
         TPhysicalFunction[i][2][k*6]  =coord[2]*shapederivative[k][1];
         TPhysicalFunction[i][2][k*6+1]=coord[2]*shapederivative[k][0];
         TPhysicalFunction[i][3][k*6+1]=shapefunction[k];
         TPhysicalFunction[i][3][k*6+2]=coord[2]*shapederivative[k][1];
         TPhysicalFunction[i][4][k*6]  =shapefunction[k];
         TPhysicalFunction[i][4][k*6+2]=coord[2]*shapederivative[k][0];
       }

     }

   };

                         /* Destructor defined here is mostly called at the
                            end of program execution. If is needless in fact
                            ha! ha! ha!
                         */
   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Shell, CInterpolation, CQuadrature> :: ~PhysicalTable()
   {
     unsigned int i,j;

     for(i = 0; i < numgauss; i++) {
         for(j =0; j < numcomponents; j++)
             free(TPhysicalFunction[i][j]);
         free(TPhysicalFunction[i]);
     }
     free(TPhysicalFunction);
   };

   template<class CInterpolation, class CQuadrature>
   PhysicalTable<Strain_Stress_Shell, CInterpolation, CQuadrature>&
   PhysicalTable<Strain_Stress_Shell, CInterpolation, CQuadrature> :: Instance()
   {
      static PhysicalTable<Strain_Stress, CInterpolation, CQuadrature> physicaltable;
      return physicaltable;
   }





} // End of namespace

#endif
