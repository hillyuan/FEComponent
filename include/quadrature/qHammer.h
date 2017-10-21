
#ifndef _qhammer_h_
#define _qhammer_h_


#ifndef _quadraturebase_h_
  #include "qudratureBase.h"
#endif


namespace GNCLib
{

  struct Quadrature_Hammer2 : public QuadratureBase
  {
     Quadrature_Hammer2(unsigned int);
  };

  namespace
  {
     QuadratureBase* CreateHammer2(unsigned int n)
     {
         return new Quadrature_Hammer2(n);
     }
     const bool Hammer2Registered =
       QuadratureFactory::Instance().Register(Hammer2, CreateHammer2);
  }


  // --------Implementation-Hammer------------
  Quadrature_Hammer2::Quadrature_Hammer2(unsigned int n)
  :QuadratureBase(n,2)
  {
      switch (n) {
      case 1:
      default:
          abscissas[0]=1./3.;   abscissas[1]=1./3.;
          weights[0]=1.;
      break;
      case 3:
          abscissas[0]=1./2.;   abscissas[1]=1./2.;
          abscissas[2]=0.;      abscissas[3]=1./2.;
          abscissas[4]=1./2.;   abscissas[5]=0.;
          weights[0]=1./3.;     weights[1]=1./3.;
      break;
      case 4:
          abscissas[0]=1./3.;   abscissas[1]=1./3.;
          weights[0]=-27.0/48;
          abscissas[2]=11./15.; abscissas[3]=2./15.;
          weights[1]=25./48.;
          abscissas[4]=2./15.;  abscissas[5]=11./15.;
          weights[2]=25./48.;
          abscissas[6]=2./15.;  abscissas[7]=2./15.;
          weights[3]=25./48.;
      break;
      case 7:
          long double afa, beta;
          afa=0.05971587;  beta=0.47014206;
          abscissas[0] = 1./3.; abscissas[1] = 1./3.;
          abscissas[2] = afa;   abscissas[3] = beta;
          abscissas[4] = beta;  abscissas[5] = afa;
          abscissas[6] = beta;  abscissas[7] = beta;
          afa=0.79742699;  beta=0.10128651;
          abscissas[8] = afa;   abscissas[9] = beta;
          abscissas[10] = beta; abscissas[11] = afa;
          abscissas[12] = beta; abscissas[13] = beta;
          weights[0] = 0.225;
          weights[1] = 0.13239415;
          weights[2] = 0.13239415;
          weights[3] = 0.13239415;
          weights[4] = 0.12593918;
          weights[5] = 0.12593918;
          weights[6] = 0.12593918;
          break;

      }
  }



  struct Quadrature_Hammer3 : public QuadratureBase
  {
     Quadrature_Hammer3(unsigned int);
  };

  namespace
  {
     QuadratureBase* CreateHammer3(unsigned int n)
     {
         return new Quadrature_Hammer3(n);
     }
     const bool Hammer3Registered =
       QuadratureFactory::Instance().Register(Hammer3, CreateHammer3);
  }


  // --------Implementation-Hammer------------
  Quadrature_Hammer3::Quadrature_Hammer3(unsigned int n)
  :QuadratureBase(n,2)
  {
      switch (n) {
      case 1:
      default:
          abscissas[0]=1./4.;   abscissas[1]=1./4.;   abscissas[2]=1./4.;
          weights[0]=1.;
      break;
      case 4:
          long double afa, beta;
          afa=0.58541020;  beta=0.13819660;
          abscissas[0]=afa;     abscissas[1]=beta;   abscissas[2]=beta;
          abscissas[3]=beta;    abscissas[4]=afa;    abscissas[5]=beta;
          abscissas[6]=beta;    abscissas[7]=beta;   abscissas[8]=afa;
          abscissas[9]=beta;    abscissas[10]=beta;  abscissas[11]=beta;
          weights[0] = 1./4.;
          weights[1] = 1./4.;
          weights[2] = 1./4.;
          weights[3] = 1./4.;
      break;
      case 5:
          abscissas[0] = 1./4.; abscissas[1] = 1./4.; abscissas[2] = 1./4.;
          abscissas[3] = 1./2.; abscissas[4] = 1./6.; abscissas[5] = 1./6.;
          abscissas[6] = 1./6.; abscissas[7] = 1./2.; abscissas[8] = 1./6.;
          abscissas[9] = 1./6.; abscissas[10] = 1./6.;abscissas[11] = 1./2.;
          abscissas[12] = 1./6.;abscissas[13] = 1./6.;abscissas[14] = 1./6.;
          weights[0] = -4./5.;
          weights[1] = 9./20.;
          weights[2] = 9./20.;
          weights[3] = 9./20.;
          weights[4] = 9./20.;
          break;
      }

  };

} //End nampespace GNCLib

#endif
