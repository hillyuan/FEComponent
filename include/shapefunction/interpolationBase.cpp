#include "interpolateBase.h"
#include <iostream>

using namespace GNCLib;

  template<class CInterpolation, class CQuadrature>
  InterpolationTable<CInterpolation, CQuadrature> :: InterpolationTable()
  {
     unsigned int i,j;

     CQuadrature aquadrature;
     CInterpolation ainterpolator;

     numgauss = aquadrature.numQuadraturePoints();
     numintps = CInterpolation::numpoints;
     dimension = CInterpolation::dimension;

     //-----allocate TShapeFunction
     TShapeFunction = (double **)malloc(numgauss * sizeof(double *));
     if(TShapeFunction == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TShapeFunction[i] = (double *)malloc(numintps * sizeof(double));
         if(TShapeFunction[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
     }

     //-----allocate TShapeDerivative
     TShapeDerivative = (double ***)malloc(numgauss * sizeof(double *));
     if(TShapeDerivative == NULL) {
         fprintf(stderr, "out of memory when building interpolation table\n");
         exit(1);
     }
     for(i = 0; i < numgauss; i++) {
         TShapeDerivative[i] = (double **)malloc(numintps * sizeof(double));
         if(TShapeDerivative[i] == NULL) {
             fprintf(stderr, "out of memory when building interpolation table\n");
             exit(1);
         }
         for(j = 0; j < numintps; j++) {
             TShapeDerivative[i][j] = (double *)malloc(dimension * sizeof(double));
             if(TShapeDerivative[i][j] == NULL) {
                 fprintf(stderr, "out of memory when building interpolation table\n");
                 exit(1);
             }
         }
     }

     // Fill in TShapeFunction and TShapeDerivative
     double coord[dimension];
     for(i = 0; i < numgauss; i++) {
         aquadrature.quadraturePoint(i, coord);
         ainterpolator.ShapeFunction(coord, TShapeFunction[i]);
         ainterpolator.ShapeDerivative(coord, TShapeDerivative[i]);
     }

  };

                         /* Destructor defined here is mostly called at the
                            end of program execution. If is needless in fact
                            ha! ha! ha!
                         */
  template<class CInterpolation, class CQuadrature>
  InterpolationTable<CInterpolation, CQuadrature> :: ~InterpolationTable()
  {
     unsigned int i,j;

     for(i = 0; i < numgauss; i++)
         free(TShapeFunction[i]);
     free(TShapeFunction);

     for(i = 0; i < numgauss; i++) {
         for(j =0; j < numintps; j++)
             free(TShapeDerivative[i][j]);
         free(TShapeDerivative[i]);
     }
     free(TShapeDerivative);
  };

  template<class CInterpolation, class CQuadrature>
  InterpolationTable<CInterpolation, CQuadrature>&
  InterpolationTable<CInterpolation, CQuadrature> :: Instance()
  {
      static InterpolationTable<CInterpolation, CQuadrature> interpolation;
      return interpolation;
  }

  template<class CInterpolation, class CQuadrature>
  void
  InterpolationTable<CInterpolation, CQuadrature>
  :: ShapeFunction(unsigned int i, double shapefunction[]) const
  {
     assert(i<=numgauss);
     shapefunction=TShapeFunction[i];
  }

  template<class CInterpolation, class CQuadrature>
  void
  InterpolationTable<CInterpolation, CQuadrature>
  :: ShapeDerivative(unsigned int i, double shapederivative[][CInterpolation::dimension]) const
  {
     assert(i<=numgauss);
     shapederivative=TShapeDerivative[i];
  }



