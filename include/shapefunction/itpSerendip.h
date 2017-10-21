/*
   Purpose: Define the Gauss-Lengedre quadrature scheme.
            Global variable GaussLengendreRegistered defined here to guarantee
          the registration of this scheme to Quadratory factory.


    Date                    Programmer              Status
    ===========            ============          =============
    Oct 06,2005             YUAN Xi               Original
*/


#ifndef _itpSerendip_h_
#define _itpSerendip_h_

#ifndef _interpolatebase_h_
  #include "interpolateBase.h"
#endif


namespace GNCLib
{
    /*
        Four nodes rectanglar element
            p4--------p3
            |         |
            |         |
            |         |
            p1--------p2
    */
    struct Interpolation<Serendip,1,1,0>
    {
       static const unsigned int dimension = 2;
       static const unsigned int numpoints = 4;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };


    void
    Interpolation<Serendip,1,1,0>::ShapeFunction(double p[dimension],double v[numpoints])
    {
      v[0]=(1.-p[0])*(1.-p[1])*0.25;
      v[1]=(1.+p[0])*(1.-p[1])*0.25;
      v[2]=(1.+p[0])*(1.+p[1])*0.25;
      v[3]=(1.-p[0])*(1.+p[1])*0.25;
    }

    void
    Interpolation<Serendip,1,1,0>::ShapeDerivative(double p[dimension],double m[numpoints][dimension])
    {
      m[0][0]=(-1.+p[1])*0.25;  m[0][1]=(-1.+p[0])*0.25;
      m[1][0]=( 1.-p[1])*0.25;  m[1][1]=(-1.-p[0])*0.25;
      m[2][0]=( 1.+p[1])*0.25;  m[2][1]=( 1.+p[0])*0.25;
      m[3][0]=(-1.-p[1])*0.25;  m[3][1]=( 1.-p[0])*0.25;
    }



    template<class array>
    array
    Interpolation<Serendip,1,1,0>::ShapeFunction(const array& p)
    {
      array v(4);
      v[0]=(1.-p[0])*(1.-p[1])*0.25;
      v[1]=(1.+p[0])*(1.-p[1])*0.25;
      v[2]=(1.+p[0])*(1.+p[1])*0.25;
      v[3]=(1.-p[0])*(1.+p[1])*0.25;
      return v;
    }

    template<class array, class matrix>
    matrix
    Interpolation<Serendip,1,1,0>::ShapeDerivative(const array& p)
    {
      matrix m(4,2);
      m(0,0)=(-1.+p[1])*0.25;  m(0,1)=(-1.+p[0])*0.25;
      m(1,0)=( 1.-p[1])*0.25;  m(1,1)=(-1.-p[0])*0.25;
      m(2,0)=( 1.+p[1])*0.25;  m(2,1)=( 1.+p[0])*0.25;
      m(3,0)=(-1.-p[1])*0.25;  m(3,1)=( 1.-p[0])*0.25;
      return m;
    }


    /*
        Eight nodes rectanglar element
        Be careful on its seqences!
            p4---p7---p3
            |         |
            p8        p6
            |         |
            p1---p5---p2
    */
    struct Interpolation<Serendip,2,2,0>
    {
       static const unsigned int dimension = 2;
       static const unsigned int numpoints = 8;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };

    void
    Interpolation<Serendip,2,2,0>::ShapeFunction(double p[dimension],double v[numpoints])
    {
       v[0]=(1.-p[0])*(1.-p[1])*(-p[0]-p[1]-1.0)*0.25;
       v[1]=(1.+p[0])*(1.-p[1])*( p[0]-p[1]-1.0)*0.25;
       v[2]=(1.+p[0])*(1.+p[1])*( p[0]+p[1]-1.0)*0.25;
       v[3]=(1.-p[0])*(1.+p[1])*(-p[0]+p[1]-1.0)*0.25;
       v[4]=(1-p[0]*p[0])*(1-p[1]);
       v[5]=(1-p[1]*p[1])*(1+p[0]);
       v[6]=(1-p[0]*p[0])*(1+p[1]);
       v[7]=(1-p[1]*p[1])*(1-p[0]);
    }



    template<class array>
    array
    Interpolation<Serendip,2,2,0>::ShapeFunction(const array& p)
    {
      array v(8);
      v[0]=(1.-p[0])*(1.-p[1])*(-p[0]-p[1]-1.0)*0.25;
      v[1]=(1.+p[0])*(1.-p[1])*( p[0]-p[1]-1.0)*0.25;
      v[2]=(1.+p[0])*(1.+p[1])*( p[0]+p[1]-1.0)*0.25;
      v[3]=(1.-p[0])*(1.+p[1])*(-p[0]+p[1]-1.0)*0.25;
      v[4]=(1-p[0]*p[0])*(1-p[1]);
      v[5]=(1-p[1]*p[1])*(1+p[0]);
      v[6]=(1-p[0]*p[0])*(1+p[1]);
      v[7]=(1-p[1]*p[1])*(1-p[0]);

      return v;
    }

    /*
        Twelve nodes rectanglar element
        Be careful on its seqences!
            p4---p10--p9--p3
            |             |
            p11           p8
            |             |
            p12           p7
            |             |
            p1---p5---p6--p2
    */
    struct Interpolation<Serendip,3,3,0>
    {
       static const unsigned int dimension = 2;
       static const unsigned int numpoints = 12;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };


    void
    Interpolation<Serendip,3,3,0>::ShapeFunction(double p[dimension],double v[numpoints])
    {
      v[0]=(1.-p[0])*(1.-p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[1]=(1.+p[0])*(1.-p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[2]=(1.+p[0])*(1.+p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[3]=(1.-p[0])*(1.+p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;

      v[4] =(1-p[1])*(1-p[0]*p[0])*(1.-3.*p[1])*0.28125;
      v[5] =(1-p[1])*(1-p[0]*p[0])*(1.+3.*p[1])*0.28125;
      v[6] =(1+p[0])*(1-p[1]*p[1])*(1.-3.*p[0])*0.28125;
      v[7] =(1+p[0])*(1-p[1]*p[1])*(1.+3.*p[0])*0.28125;
      v[8] =(1+p[1])*(1-p[0]*p[0])*(1.+3.*p[1])*0.28125;
      v[9] =(1+p[1])*(1-p[0]*p[0])*(1.-3.*p[1])*0.28125;
      v[10]=(1-p[0])*(1-p[1]*p[1])*(1.+3.*p[1])*0.28125;
      v[11]=(1-p[0])*(1-p[1]*p[1])*(1.-3.*p[1])*0.28125;
    }


    template<class array>
    array
    Interpolation<Serendip,3,3,0>::ShapeFunction(const array& p)
    {
      array v(12);
      v[0]=(1.-p[0])*(1.-p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[1]=(1.+p[0])*(1.-p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[2]=(1.+p[0])*(1.+p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;
      v[3]=(1.-p[0])*(1.+p[1])*(-10.+9.*(p[0]*p[0]+p[1]*p[1]))*0.03125;

      v[4] =(1-p[1])*(1-p[0]*p[0])*(1.-3.*p[1])*0.28125;
      v[5] =(1-p[1])*(1-p[0]*p[0])*(1.+3.*p[1])*0.28125;
      v[6] =(1+p[0])*(1-p[1]*p[1])*(1.-3.*p[0])*0.28125;
      v[7] =(1+p[0])*(1-p[1]*p[1])*(1.+3.*p[0])*0.28125;
      v[8] =(1+p[1])*(1-p[0]*p[0])*(1.+3.*p[1])*0.28125;
      v[9] =(1+p[1])*(1-p[0]*p[0])*(1.-3.*p[1])*0.28125;
      v[10]=(1-p[0])*(1-p[1]*p[1])*(1.+3.*p[1])*0.28125;
      v[11]=(1-p[0])*(1-p[1]*p[1])*(1.-3.*p[1])*0.28125;

      return v;
    }

    /*
       Eight nodes hexahedral element
              p8-------p7
            / |       / |
           p4--------p3 |            ^x3
           |  |      |  |            |
           | p5------|-p6            |   / x2
           | /       |/              | /
           p1--------p2              |-----> x1
    */
    struct Interpolation<Serendip,1,1,1>
    {
       static const unsigned int dimension = 3;
       static const unsigned int numpoints = 8;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };

   void
   Interpolation<Serendip,1,1,1>::ShapeFunction(double p[dimension],double v[numpoints])
   {
      v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*0.125;
      v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*0.125;
      v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*0.125;
      v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*0.125;
      v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*0.125;
      v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*0.125;
      v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*0.125;
      v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*0.125;
   }



    template<class array>
    array
    Interpolation<Serendip,1,1,1>::ShapeFunction(const array& p)
    {
      array v(8);
      v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*0.125;
      v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*0.125;
      v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*0.125;
      v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*0.125;
      v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*0.125;
      v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*0.125;
      v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*0.125;
      v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*0.125;

      return v;
    }

    /*
      Twenty nodes hexahedral element
               p8----p19---p7
              /|          /|
           p16 |       p15 |
           /  p20      /  p18
          p4---p11----p3   |
          |    |      |    |
          |    p5----p17---p6        ^ x3
         p12   /     p10  /          |
          | p13       | p14          |
          |/          |/             |
          p1----p9 ---p2             -----> x1
   */
   struct Interpolation<Serendip,2,2,2>
   {
       static const unsigned int dimension = 3;
       static const unsigned int numpoints = 20;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);

   };


   void
   Interpolation<Serendip,2,2,2>::ShapeFunction(double p[dimension],double v[numpoints])
   {
     v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*(-p[0]-p[1]-p[2]-2.)*0.125;
     v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*( p[0]-p[1]-p[2]-2.)*0.125;
     v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*( p[0]-p[1]+p[2]-2.)*0.125;
     v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*(-p[0]-p[1]+p[2]-2.)*0.125;
     v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*(-p[0]+p[1]-p[2]-2.)*0.125;
     v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*( p[0]+p[1]-p[2]-2.)*0.125;
     v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*( p[0]+p[1]+p[2]-2.)*0.125;
     v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*(-p[0]+p[1]+p[2]-2.)*0.125;

     v[8] =(1.-p[0]*p[0])*(1.-p[1])*(1.-p[2])*0.25;
     v[9] =(1.-p[2]*p[2])*(1.-p[1])*(1.+p[0])*0.25;
     v[10]=(1.-p[0]*p[0])*(1.-p[1])*(1.+p[2])*0.25;
     v[11]=(1.-p[2]*p[2])*(1.-p[1])*(1.-p[0])*0.25;

     v[12]=(1.-p[1]*p[1])*(1.-p[0])*(1.-p[2])*0.25;
     v[13]=(1.-p[1]*p[1])*(1.+p[0])*(1.-p[2])*0.25;
     v[14]=(1.-p[1]*p[1])*(1.+p[0])*(1.+p[2])*0.25;
     v[15]=(1.-p[1]*p[1])*(1.-p[0])*(1.+p[2])*0.25;

     v[16]=(1.-p[0]*p[0])*(1.+p[1])*(1.-p[2])*0.25;
     v[17]=(1.-p[2]*p[2])*(1.+p[1])*(1.+p[0])*0.25;
     v[18]=(1.-p[0]*p[0])*(1.+p[1])*(1.+p[2])*0.25;
     v[19]=(1.-p[2]*p[2])*(1.+p[1])*(1.-p[0])*0.25;
   }



   template<class array>
   array
   Interpolation<Serendip,2,2,2>::ShapeFunction(const array& p)
   {
     array v(20);
     v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*(-p[0]-p[1]-p[2]-2.)*0.125;
     v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*( p[0]-p[1]-p[2]-2.)*0.125;
     v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*( p[0]-p[1]+p[2]-2.)*0.125;
     v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*(-p[0]-p[1]+p[2]-2.)*0.125;
     v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*(-p[0]+p[1]-p[2]-2.)*0.125;
     v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*( p[0]+p[1]-p[2]-2.)*0.125;
     v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*( p[0]+p[1]+p[2]-2.)*0.125;
     v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*(-p[0]+p[1]+p[2]-2.)*0.125;

     v[8] =(1.-p[0]*p[0])*(1.-p[1])*(1.-p[2])*0.25;
     v[9] =(1.-p[2]*p[2])*(1.-p[1])*(1.+p[0])*0.25;
     v[10]=(1.-p[0]*p[0])*(1.-p[1])*(1.+p[2])*0.25;
     v[11]=(1.-p[2]*p[2])*(1.-p[1])*(1.-p[0])*0.25;

     v[12]=(1.-p[1]*p[1])*(1.-p[0])*(1.-p[2])*0.25;
     v[13]=(1.-p[1]*p[1])*(1.+p[0])*(1.-p[2])*0.25;
     v[14]=(1.-p[1]*p[1])*(1.+p[0])*(1.+p[2])*0.25;
     v[15]=(1.-p[1]*p[1])*(1.-p[0])*(1.+p[2])*0.25;

     v[16]=(1.-p[0]*p[0])*(1.+p[1])*(1.-p[2])*0.25;
     v[17]=(1.-p[2]*p[2])*(1.+p[1])*(1.+p[0])*0.25;
     v[18]=(1.-p[0]*p[0])*(1.+p[1])*(1.+p[2])*0.25;
     v[19]=(1.-p[2]*p[2])*(1.+p[1])*(1.-p[0])*0.25;


     return v;
   }


   /*
      36 nodes hexahedral element

                 p8--p30--p29--p7
                /|           / |
            p24 p31        p23 p28
           p20   |        p19  |
           /    p32      /     |
          p4--p14--p13--p3     p27
          |      |       |     |
          |     p5--p25--p25--p6        ^ x3
        p15    /       p12    /         |
          |  p21         |  p22         |
        p16 p17        p11 p18          |
          |/             |/             |
          p1--p9 --p10--p2               -----> x1
   */
   struct Interpolation<Serendip,3,3,3>
   {
       static const unsigned int dimension = 3;
       static const unsigned int numpoints = 32;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
   };


   void
   Interpolation<Serendip,3,3,3>::ShapeFunction(double p[dimension],double v[numpoints])
   {
     v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;

     v[8] =(1.-p[0]*p[0])*(1-3.*p[0])*(1.-p[1])*(1.-p[2])*0.140625;
     v[9] =(1.-p[0]*p[0])*(1+3.*p[0])*(1.-p[1])*(1.-p[2])*0.140625;
     v[10]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.-p[1])*(1.+p[0])*0.140625;
     v[11]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.-p[1])*(1.+p[0])*0.140625;
     v[12]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.-p[1])*(1.+p[2])*0.140625;
     v[13]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.-p[1])*(1.+p[2])*0.140625;
     v[14]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.-p[1])*(1.-p[0])*0.140625;
     v[15]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.-p[1])*(1.-p[0])*0.140625;

     v[16]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.-p[0])*(1.-p[2])*0.140625;
     v[17]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.+p[0])*(1.-p[2])*0.140625;
     v[18]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.+p[0])*(1.+p[2])*0.140625;
     v[19]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.-p[0])*(1.+p[2])*0.140625;

     v[20]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.-p[0])*(1.-p[2])*0.140625;
     v[21]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.+p[0])*(1.-p[2])*0.140625;
     v[22]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.+p[0])*(1.+p[2])*0.140625;
     v[23]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.-p[0])*(1.+p[2])*0.140625;

     v[24]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.+p[1])*(1.-p[2])*0.140625;
     v[25]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.+p[1])*(1.-p[2])*0.140625;
     v[26]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.+p[1])*(1.+p[0])*0.140625;
     v[27]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.+p[1])*(1.+p[0])*0.140625;
     v[28]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.+p[1])*(1.+p[2])*0.140625;
     v[29]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.+p[1])*(1.+p[2])*0.140625;
     v[30]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.+p[1])*(1.-p[0])*0.140625;
     v[31]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.+p[1])*(1.-p[0])*0.140625;
   }



   template<class array>
   array
   Interpolation<Serendip,3,3,3>::ShapeFunction(const array& p)
   {
     array v(32);
     v[0]=(1.-p[0])*(1.-p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[1]=(1.+p[0])*(1.-p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[2]=(1.+p[0])*(1.-p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[3]=(1.-p[0])*(1.-p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[4]=(1.-p[0])*(1.+p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[5]=(1.+p[0])*(1.+p[1])*(1.-p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[6]=(1.+p[0])*(1.+p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;
     v[7]=(1.-p[0])*(1.+p[1])*(1.+p[2])*(9.*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])-19.)*0.015625;

     v[8] =(1.-p[0]*p[0])*(1-3.*p[0])*(1.-p[1])*(1.-p[2])*0.140625;
     v[9] =(1.-p[0]*p[0])*(1+3.*p[0])*(1.-p[1])*(1.-p[2])*0.140625;
     v[10]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.-p[1])*(1.+p[0])*0.140625;
     v[11]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.-p[1])*(1.+p[0])*0.140625;
     v[12]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.-p[1])*(1.+p[2])*0.140625;
     v[13]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.-p[1])*(1.+p[2])*0.140625;
     v[14]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.-p[1])*(1.-p[0])*0.140625;
     v[15]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.-p[1])*(1.-p[0])*0.140625;

     v[16]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.-p[0])*(1.-p[2])*0.140625;
     v[17]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.+p[0])*(1.-p[2])*0.140625;
     v[18]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.+p[0])*(1.+p[2])*0.140625;
     v[19]=(1.-p[1]*p[1])*(1-3.*p[1])*(1.-p[0])*(1.+p[2])*0.140625;

     v[20]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.-p[0])*(1.-p[2])*0.140625;
     v[21]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.+p[0])*(1.-p[2])*0.140625;
     v[22]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.+p[0])*(1.+p[2])*0.140625;
     v[23]=(1.-p[1]*p[1])*(1+3.*p[1])*(1.-p[0])*(1.+p[2])*0.140625;

     v[24]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.+p[1])*(1.-p[2])*0.140625;
     v[25]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.+p[1])*(1.-p[2])*0.140625;
     v[26]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.+p[1])*(1.+p[0])*0.140625;
     v[27]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.+p[1])*(1.+p[0])*0.140625;
     v[28]=(1.-p[0]*p[0])*(1+3.*p[0])*(1.+p[1])*(1.+p[2])*0.140625;
     v[29]=(1.-p[0]*p[0])*(1-3.*p[0])*(1.+p[1])*(1.+p[2])*0.140625;
     v[30]=(1.-p[2]*p[2])*(1+3.*p[2])*(1.+p[1])*(1.-p[0])*0.140625;
     v[31]=(1.-p[2]*p[2])*(1-3.*p[2])*(1.+p[1])*(1.-p[0])*0.140625;

     return v;
   }





}  //End namespcae GNCLib

#endif
