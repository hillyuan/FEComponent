/*
   Purpose: Define the Gauss-Lengedre quadrature scheme.
            Global variable GaussLengendreRegistered defined here to guarantee
          the registration of this scheme to Quadratory factory.


    Date                    Programmer              Status
    ===========            ============          =============
    July 11,2005             YUAN Xi               Original
*/


#ifndef _qgausslegendre_h_
#define _qgausslegendre_h_

#ifndef _quadraturebase_h_
  #include "quadratureBase.h"
#endif


namespace GNCLib
{
  struct Quadrature_GaussLegendre : public QuadratureBase
  {
      Quadrature_GaussLegendre(unsigned int);
  };

  namespace
  {
      QuadratureBase* CreateGaussLegendre(unsigned int n)
      {
          return new Quadrature_GaussLegendre(n);
      }
      const bool GaussLengendreRegistered =
        QuadratureFactory::Instance().Register(GaussLegendre, CreateGaussLegendre);
      const bool AssumedStrain_ShearRegistered =
        QuadratureFactory::Instance().Register(AssumedStrain_Shear, CreateGaussLegendre);
  }


  // --------Implementation-Legendre------------
  Quadrature_GaussLegendre::Quadrature_GaussLegendre(unsigned int n)
  :QuadratureBase(n,1)
  {
     switch (n) {
     case 1:
         weights[0]   = 2.0;
         abscissas[0] = 0.0;
     break;
     case 2:
         abscissas[0]=-0.577350269189626;    weights[0]=1.0;
         abscissas[1]= 0.577350269189626;    weights[1]=1.0;
     break;
     case 3:
         abscissas[0]=-0.774596669241483;    weights[0]=0.555555555555556;
         abscissas[1]= 0.0;                  weights[1]=0.888888888888889;
         abscissas[2]= 0.774596669241483;    weights[2]=0.555555555555556;
     break;
     case 4:
         abscissas[0]=-0.861136311594053;    weights[0]=0.347854845137454;
         abscissas[1]=-0.339981043584856;    weights[1]=0.652145154862546;
         abscissas[2]=-abscissas[0];         weights[2]=weights[0];
         abscissas[3]=-abscissas[1];         weights[3]=weights[1];
     break;
     case 5:
         abscissas[0]=-0.906179845938664;    weights[0]=0.236926885056189;
         abscissas[1]=-0.538469310105683;    weights[1]=0.478628670499366;
         abscissas[2]= 0.0;                  weights[2]=0.568888888888889;
         abscissas[3]=-abscissas[0];         weights[3]=weights[0];
         abscissas[4]=-abscissas[1];         weights[4]=weights[1];
     break;
     }
  };
}  //End namespcae GNCLib

#endif
