/*
    Purpose: Define the relation between geomatric convex, or cell, and
         quadrature scheme and all those quadrature scheme are define as
         singleton of Loki.
            Great flexibility and extendability are obtained from this
         method.


       Date                    Programmer              Status
        ===========            ============          =============
        July 16,2005             YUAN Xi               Original
*/

#include "../loki/Singleton.h"


enum ConvexType
{
   LNE2,
   LNE3,
   TRI3,
   TRI3_SHELL,
   REC4,
   REC4_SHELL,
   TET4,
   PRI6,
   HEX8
};

namespace GNCLib
{
  template<ConvexType etype>
  class ConvexQuadrature  {};

  template<>
  class ConvexQuadrature<LNE2>
  {
    typedef Loki::SingletonHolder
            < Quadrature<1,1,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<LNE3>
  {
    typedef Loki::SingletonHolder
            < Quadrature<1,2,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<TRI3>
  {
    typedef Loki::SingletonHolder
            < Quadrature<2,1,Hammer2> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<REC4>
  {
    typedef Loki::SingletonHolder
            < Quadrature<2,2,GaussLegendre,2,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<TET4>
  {
    typedef Loki::SingletonHolder
            < Quadrature<3,1,Hammer3> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<PRI6>
  {
    typedef Loki::SingletonHolder
            < Quadrature<3,1,Hammer2,2,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<HEX8>
  {
    typedef Loki::SingletonHolder
            < Quadrature<3,2,GaussLegendre,2,GaussLegendre,2,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<TRI3_SHELL>
  {
    typedef Loki::SingletonHolder
            < Quadrature<3,1,Hammer2,2,GaussLegendre> >
            QuadratureObject;
  };

  template<>
  class ConvexQuadrature<REC4_SHELL>
  {
    typedef Loki::SingletonHolder
            < Quadrature<3,2,GaussLegendre,2,GaussLegendre,2,GaussLegendre> >
            QuadratureObject;
  };


}
