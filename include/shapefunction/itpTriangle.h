/*
   Purpose: Define the Gauss-Lengedre quadrature scheme.
            Global variable GaussLengendreRegistered defined here to guarantee
          the registration of this scheme to Quadratory factory.


    Date                    Programmer              Status
    ===========            ============          =============
    Oct 06,2005             YUAN Xi               Original
*/


#ifndef _itpTriangle_h_
#define _itpTriangle_h_

#ifndef _interpolatebase_h_
  #include "interpolateBase.h"
#endif


namespace GNCLib
{
    /*
        Three nodes triangle element
                      p3
                    /  |
                  /    |
                /      |
              /        |
            p1--------p2
    */
    struct Interpolation<Triangle,1,0,0>
    {
       static const unsigned int dimension = 2;
       void ShapeFunction(double[dimension],double[]);
       void ShapeDerivative(double[dimension],double[][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);

    };

    void
    Interpolation<Triangle,1,0,0>::ShapeFunction(double p[dimension],double v[])
    {
      v[0]=1.-p[0]-p[1];
      v[1]=p[0];
      v[2]=p[1];
    }

    void
    Interpolation<Triangle,1,0,0>::ShapeDerivative(double p[dimension],double m[][dimension])
    {
      m[0][0]= 1.;  m[0][1]= 0.;
      m[1][0]= 0.;  m[1][1]= 1.;
      m[2][0]=-1.;  m[2][1]=-1.;
    }


    template<class array>
    array
    Interpolation<Triangle,1,0,0>::ShapeFunction(const array& p)
    {
      array v(3);
      v[0]=1.-p[0]-p[1];
      v[1]=p[0];
      v[2]=p[1];
      return v;
    }


    template<class array, class matrix>
    matrix
    Interpolation<Triangle,1,0,0>::ShapeDerivative(const array& p)
    {
      matrix m(3,2);
      m(0,0)= 1.;  m(0,1)= 0.;
      m(1,0)= 0.;  m(1,1)= 1.;
      m(2,0)=-1.;  m(2,1)=-1.;
      return m;
    }


    /*
       Five nodes triangle element
                     p3
                   /  |
                 /    |
              p6     p5
             /        |
           p1---p4---p2
   */
   struct Interpolation<Triangle,2,0,0>
   {

       static const unsigned int dimension = 2;
       void ShapeFunction(double[dimension],double[]);
       void ShapeDerivative(double[dimension],double[][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);

   };


   void
   Interpolation<Triangle,2,0,0>::ShapeFunction(double p[dimension],double v[])
   {
     double p0=1.-p[0]-p[1];
     v[0]=(2.*p0-1.)*p0;
     v[1]=(2.*p[0]-1.)*p[0];
     v[2]=(2.*p[1]-1.)*p[1];
     v[3]=4.*p0*p[0];
     v[4]=4.*p[0]*p[1];
     v[5]=4.*p[1]*p0;
   }


   template<class array>
   array
   Interpolation<Triangle,2,0,0>::ShapeFunction(const array& p)
   {
     array v(6);
     double p0=1.-p[0]-p[1];
     v[0]=(2.*p0-1.)*p0;
     v[1]=(2.*p[0]-1.)*p[0];
     v[2]=(2.*p[1]-1.)*p[1];
     v[3]=4.*p0*p[0];
     v[4]=4.*p[0]*p[1];
     v[5]=4.*p[1]*p0;
     return v;
   }


   template<class array, class matrix>
   matrix
   Interpolation<Triangle,2,0,0>::ShapeDerivative(const array& p)
   {
     matrix m(3,2);
     m(0,0)= 1.;  m(0,1)= 0.;
     m(1,0)= 0.;  m(1,1)= 1.;
     m(2,0)=-1.;  m(2,1)=-1.;
     return m;
   }

}

#endif
