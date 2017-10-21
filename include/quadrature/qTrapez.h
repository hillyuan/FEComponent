#ifndef _qtrapez_h_
#define _qtrapez_h_


#ifndef _quadraturebase_h_
  #include "qudratureBase.h"
#endif


namespace GNCLib
{

  struct Quadrature_Trapez : public QuadratureBase
  {
     Quadrature_Trapez(unsigned int n=2);
  };

  namespace
  {
     QuadratureBase* CreateTrapez(unsigned int n)
     {
         return new Quadrature_Trapez(n);
     }
     const bool TrapezRegistered =
       QuadratureFactory::Instance().Register(Trapez, CreateTrapez);
  }


  // --------Implementation-Hammer------------
  Quadrature_Trapez::Quadrature_Trapez(unsigned int n)
  :QuadratureBase(2,1)
  {
      abscissas[0]=0.0;   abscissas[1]=1.0;
      weights[0]=0.5;     weights[1]=0.5;
  }

}

#endif
