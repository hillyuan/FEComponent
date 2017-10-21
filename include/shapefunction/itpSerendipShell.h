/*
   Purpose: Interpolate functions for shell elements.

   P.A. : last function argument is supposed thickness of shell


    Date                    Programmer              Status
    ===========            ============          =============
    Oct 06,2005             YUAN Xi               Original
*/


#ifndef _itpSerendipShell_h_
#define _itpSerendipShell_h_

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
    template<unsigned int h>
    struct Interpolation<SerendipShell,1,1,h>
    {
       static const unsigned int dimension = 3;
       static const unsigned int numpoints = 4;

       void ShapeFunction(double[dimension],double[numpoints]);
       void ShapeDerivative(double[dimension],double[numpoints][dimension]);

       template<class array>
       array ShapeFunction(const array&);
       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };


    void
    Interpolation<SerendipShell,1,1,0>::ShapeFunction(double p[dimension],double v[numpoints])
    {
      v[0]=(1.-p[0])*(1.-p[1])*0.25;
      v[1]=(1.+p[0])*(1.-p[1])*0.25;
      v[2]=(1.+p[0])*(1.+p[1])*0.25;
      v[3]=(1.-p[0])*(1.+p[1])*0.25;
    }

    void
    Interpolation<SerendipShell,1,1,0>::ShapeDerivative(double p[dimension],double m[numpoints][dimension])
    {
      m[0][0]=(-1.+p[1])*0.25;  m[0][1]=(-1.+p[0])*0.25;
      m[1][0]=( 1.-p[1])*0.25;  m[1][1]=(-1.-p[0])*0.25;
      m[2][0]=( 1.+p[1])*0.25;  m[2][1]=( 1.+p[0])*0.25;
      m[3][0]=(-1.-p[1])*0.25;  m[3][1]=( 1.-p[0])*0.25;
    }



    template<unsigned int h>
    template<class array>
    array
    Interpolation<SerendipShell,1,1,h>::ShapeFunction(const array& p)
    {
      array v(4);
      v[0]=(1.-p[0])*(1.-p[1])*0.25;
      v[1]=(1.+p[0])*(1.-p[1])*0.25;
      v[2]=(1.+p[0])*(1.+p[1])*0.25;
      v[3]=(1.-p[0])*(1.+p[1])*0.25;
      return v;
    }


    template<unsigned int h>
    template<class array, class matrix>
    matrix
    Interpolation<SerendipShell,1,1,h>::ShapeDerivative(const array& p)
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
    template<unsigned int h>
    struct Interpolation<SerendipShell,2,2,h>
    {
       template<class array>
       array ShapeFunction(const array&);

       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };


    template<unsigned int h>
    template<class array>
    array
    Interpolation<SerendipShell,2,2,h>::ShapeFunction(const array& p)
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
    template<unsigned int h>
    struct Interpolation<SerendipShell,3,3,h>
    {
       template<class array>
       array ShapeFunction(const array&);

       template<class array, class matrix>
       matrix ShapeDerivative(const array&);
    };


    template<unsigned int h>
    template<class array>
    array
    Interpolation<SerendipShell,3,3,h>::ShapeFunction(const array& p)
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



}  //End namespcae GNCLib

#endif
