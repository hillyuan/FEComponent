/*
    Purpose: The class Quadrature defined here provide a general tool for
          numerical integration of a general function, including member
          function of a general class, up to five dimensions. Improper
          integrals are not considered currently. More general class as
          that consider function of different kind of viariables and to
          any dimensions maybe complished with Loki::Functor or Boost::
          Function but not considered here.
            Class Quadrature is templated class with template parameters
          denoting the rank of integration of different dimensions. The
          reason of adopting template is that in general application only
          one instance of this class is needed, or it is a Singleton.
          Templates provide a strong type definition that supposing Singlton.
          The backgroud integration scheme, such as Gauss-Legendre, or
          Gauss-Hermite, are transparent to users. It just known by Quadrature
          Factory.
            Argument type of functions to be integrated can be specified by user
          provide that computation of
              double * user defined function(user defined funcion argument)
          is defined. In this case, user can deal with problems with number
          of arguments of the function greater than 5. Default argument types is
          list of double value. User defined coordinates can also be acqired in
          the case that user defined coordinate class has constructor of
              coordinate(double,double,double,double,double)


       Date                    Programmer              Status
	===========            ============          =============
	July 11,2005             YUAN Xi               Original
*/


#ifndef _quadraturebase_h_
 #include "quadratureBase.h"
#endif
#include <iostream>
#include <stdlib.h>
#include <assert.h>


namespace GNCLib
{

  template
  <
     unsigned int dim,
     unsigned int n1,    QuadratureScheme q1=GaussLegendre,
     unsigned int n2=0,  QuadratureScheme q2=GaussLegendre,
     unsigned int n3=0,  QuadratureScheme q3=GaussLegendre,
     unsigned int n4=0,  QuadratureScheme q4=GaussLegendre,
     unsigned int n5=0,  QuadratureScheme q5=GaussLegendre
  >
  class Quadrature
  {
    private:
      unsigned int               ngauss;
      double**                   abscissas;
      double*                    weights;

    public:
      Quadrature();
      ~Quadrature();

      //Integration of user provided coordinate type function
      template<typename T, class coordinate>
      T Integrate( T (*)(coordinate) ) const;
      template<typename T, class C, class coordinate>
      T Integrate( C, T (C::*)(coordinate) ) const;


      //Integration specified for general function with upto 5 arguments
      template<typename T>
      T Integrate( T (*)(double) ) const;
      template<typename T>
      T Integrate( T (*)(double,double) ) const;
      template<typename T>
      T Integrate( T (*)(double,double,double) ) const;
      template<typename T>
      T Integrate( T (*)(double,double,double,double) ) const;
      template<typename T>
      T Integrate( T (*)(double,double,double,double,double) ) const;

      //Integration specified for general member function with upto 5 arguments
      template<typename T, class C>
      T Integrate( C, T (C::*)(double) ) const;
      template<typename T, class C>
      T Integrate( C, T (C::*)(double,double) ) const;
      template<typename T, class C>
      T Integrate( C, T (C::*)(double,double,double) ) const;
      template<typename T, class C>
      T Integrate( C, T (C::*)(double,double,double,double) ) const;
      template<typename T, class C>
      T Integrate( C, T (C::*)(double,double,double,double,double) ) const;

      //Method to fetch internal data
      int numQuadraturePoints() const {return ngauss;}
      int numQuadratureDimension() const {return dim;}

                            // Return the weights' value of a point denoted.
      double quadratureWeight(unsigned int i) const {return weights[i];}
                            /* This function give a way to return a user defined
                               coordinate. The MUST BE condition is that this
                               coordinate class has construct like
                                       coordinate(double,double,...)
                               It should be noted that, however, it is not a
                              efficient one. Because we should construct and
                              destruct a local coordinate class and construct
                              a coordinate class to return. Following function
                              with the same function is much faster.
                            */
      template<class coordinate>
      coordinate quadraturePoint(const unsigned int) const;
                           // Return the coordinate value of a point denoted.
      void quadraturePoint(unsigned int,double [dim]) const;
  };

  template
  <
     unsigned int dim,
     unsigned int n1,  QuadratureScheme q1,
     unsigned int n2,  QuadratureScheme q2,
     unsigned int n3,  QuadratureScheme q3,
     unsigned int n4,  QuadratureScheme q4,
     unsigned int n5,  QuadratureScheme q5
  >
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>::~Quadrature()
  {
    for(unsigned int i = 0; i < ngauss; i++)
                free(abscissas[i]);
    free(abscissas);

    free(weights);
  }



// ----------Implementation-Quadrature-----------------
                /*
                   Only dim upto 5 is considered. Maybe lack of scalbility but
                   no effective method could be found yet.
                   Algorithm of tensor product is adopted.
                */

  template
  <
     unsigned int dim,
     unsigned int n1,  QuadratureScheme q1,
     unsigned int n2,  QuadratureScheme q2,
     unsigned int n3,  QuadratureScheme q3,
     unsigned int n4,  QuadratureScheme q4,
     unsigned int n5,  QuadratureScheme q5
  >
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5> :: Quadrature()
  {
    unsigned int i,j,k,l,dimension,coorddim;
    unsigned int orders[5];
    QuadratureScheme qs[5];

    if(n1>0) {
        dimension=1;   orders[0]=n1;  qs[0]=q1;
    }
    if(n2>0) {
        dimension=2;   orders[1]=n2;  qs[1]=q2;
    }
    if(n3>0) {
        dimension=3;   orders[2]=n3;  qs[2]=q3;
    }
    if(n4>0) {
        dimension=4;   orders[3]=n4;  qs[3]=q4;
    }
    if(n5>0) {
        dimension=5;   orders[4]=n5;  qs[4]=q5;
    }
 //   assert(dimension==dim);

    ngauss=1;
    for(i=0; i<dimension; i++) {
        assert(orders[i]>0);
        ngauss = ngauss * orders[i];
    }

    weights = (double *)malloc(ngauss * sizeof(double));
    if(weights == NULL) {
        fprintf(stderr, "out of memory in quadrature definition\n");
        exit(1);
    }
    abscissas = (double **)malloc(ngauss * sizeof(double *));
    if(abscissas == NULL) {
         fprintf(stderr, "out of memory in quadrature definition\n");
         exit(1);
    }
    for(i = 0; i < ngauss; i++) {
         abscissas[i] = (double *)malloc(dim * sizeof(double));
         if(abscissas[i] == NULL) {
             fprintf(stderr, "out of memory in quadrature definition\n");
             exit(1);
         }
    }

    double curweight;             //weights in current dimension
    double weighttemp[ngauss];    //intermedia result of weights
    double curcoord[dim];         //coord of points in curr dim
    double coordtemp[ngauss][dim];
    unsigned int curorder, curpoint, curdimpos, pointcount;

    QuadratureBase* qbase;
    //prepare front of weight and coordinate
    curorder=orders[0];
    qbase=QuadratureFactory::Instance().CreateQuadrature(qs[0],orders[0]);
    coorddim=qbase->coordDim;
    for(j=0; j<curorder; j++) {
        weighttemp[j]=qbase->weights[j];
        weights[j]=weighttemp[j];
        for(l=0; l<coorddim; l++) {
            coordtemp[j][l]=qbase->abscissas[j*coorddim+l];
            abscissas[j][l]=qbase->abscissas[j*coorddim+l];
        }
    }
    curpoint = curorder;
    curdimpos = coorddim;

    for(i=1; i<dimension; i++) {
         curorder=orders[i];
         qbase=QuadratureFactory::Instance().CreateQuadrature(qs[i],orders[i]);
         coorddim=qbase->coordDim;

         pointcount=0;
         for(j=0; j<curorder; j++) {
             curweight=qbase->weights[j];
             for(l=0; l<coorddim; l++) {
                 curcoord[l]=qbase->abscissas[j*coorddim+l];
             }
             //tensor product of weighttemp and curweight, coordtemp and curcoord
             for(k=0; k<curpoint; k++) {
                 for(l=0; l<curdimpos; l++) {
                     abscissas[pointcount][l]=coordtemp[k][l];
                 }
                 for(l=0; l<coorddim; l++) {
                     abscissas[pointcount][curdimpos+l]=curcoord[l];
                 }
                 weights[pointcount]=curweight*weighttemp[k];
                 pointcount++;
             }
         }

         // restore above tensor products into coordtemp and weighttemp
         curpoint = pointcount;
         curdimpos = curdimpos + coorddim;
         for(j=0; j<curpoint; j++) {
             weighttemp[j]=weights[j];
             for(k=0; k<curdimpos; k++) {
                 coordtemp[j][k] = abscissas[j][k];
             }
         }
     }

     for(i=0; i<ngauss; i++) {
     //    weights[i]   = arrtemp1[i];
         std::cout << i<<"   "<<weights[i] <<"   "<<abscissas[i][0]<<"   "<<abscissas[i][1]<<std::endl;
     }

  };


  template
  <
     unsigned int dim,
     unsigned int n1,  QuadratureScheme q1,
     unsigned int n2,  QuadratureScheme q2,
     unsigned int n3,  QuadratureScheme q3,
     unsigned int n4,  QuadratureScheme q4,
     unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class coordinate>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(coordinate) ) const
  {
      coordinate a;

      T integral=0.0;
      for(int i=0;i<ngauss;i++) {
          a=abscissas[i];
          integral=integral+weights[i]*func(a);
      }

      return integral;
  }

  template
  <
     unsigned int dim,
     unsigned int n1,  QuadratureScheme q1,
     unsigned int n2,  QuadratureScheme q2,
     unsigned int n3,  QuadratureScheme q3,
     unsigned int n4,  QuadratureScheme q4,
     unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C, class coordinate>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate(C aObject,  T (C::*func)(coordinate) ) const
  {
      coordinate a;

      T integral=0.0;
      for(int i=0;i<ngauss;i++) {
          a=abscissas[i];
          integral=integral+weights[i]*(aObject.*func)(a);
      }

      return integral;
  }




  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(double) ) const
  {
     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*func(abscissas[i][0]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(double,double) ) const
  {
     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*func(abscissas[i][0],abscissas[i][1]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(double,double,double) ) const
  {
     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*func(abscissas[i][0],abscissas[i][1],abscissas[i][2]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(double,double,double,double) ) const
  {
     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*func(abscissas[i][0],abscissas[i][1],
                                           abscissas[i][2],abscissas[i][3]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( T (*func)(double,double,double,double,double) ) const
  {
     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*func(abscissas[i][0],abscissas[i][1],abscissas[i][2],
                                           abscissas[i][3],abscissas[i][4]);
     }

     return integral;
  }




  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( C aObject,  T (C::*func)(double) ) const
  {
     assert(dim>=1);

     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*(aObject.*func)(abscissas[i][0]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( C aObject,  T (C::*func)(double,double) ) const
  {
     assert(dim>=2);

     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*(aObject.*func)(abscissas[i][0],abscissas[i][1]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate(C aObject,  T (C::*func)(double,double,double) ) const
  {
     assert(dim>=3);

     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*(aObject.*func)(abscissas[i][0],abscissas[i][1],
                                                                  abscissas[i][2]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( C aObject,  T (C::*func)(double,double,double,double)) const
  {
     assert(dim>=4);

     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
         integral=integral+weights[i]*(aObject.*func)(abscissas[i][0],abscissas[i][1],
                                           abscissas[i][2],abscissas[i][3]);
     }

     return integral;
  }

  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<typename T, class C>
  T
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::Integrate( C aObject,  T (C::*func)(double,double,double,double,double) ) const
  {
     assert(dim>=5);

     T integral=0.0;
     for(int i=0;i<ngauss;i++) {
      //   a=abscissas[i];
         integral=integral+weights[i]*(aObject.*func)(abscissas[i][0],abscissas[i][1],
                                      abscissas[i][2],abscissas[i][3],abscissas[i][4]);
     }

     return integral;
  }


  template
  <
    unsigned int dim,
    unsigned int n1,  QuadratureScheme q1,
    unsigned int n2,  QuadratureScheme q2,
    unsigned int n3,  QuadratureScheme q3,
    unsigned int n4,  QuadratureScheme q4,
    unsigned int n5,  QuadratureScheme q5
  >
  template<class coordinate>
  coordinate
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::quadraturePoint(const unsigned int i) const
  {
      switch (dim)
      {
      case 1:
      default:
          return coordinate(abscissas[i][0]);
          break;
      case 2:
          return coordinate(abscissas[i][0],abscissas[i][1]);
          break;
      case 3:
          return coordinate(abscissas[i][0],abscissas[i][1],abscissas[i][2]);
          break;
      case 4:
          return coordinate(abscissas[i][0],abscissas[i][1],
                            abscissas[i][2],abscissas[i][3]);
          break;
      case 5:
          return coordinate(abscissas[i][0],abscissas[i][1],
                            abscissas[i][2],abscissas[i][3],abscissas[i][4]);
          break;
      }
  }


  template
  <
   unsigned int dim,
   unsigned int n1,  QuadratureScheme q1,
   unsigned int n2,  QuadratureScheme q2,
   unsigned int n3,  QuadratureScheme q3,
   unsigned int n4,  QuadratureScheme q4,
   unsigned int n5,  QuadratureScheme q5
  >
  void
  Quadrature<dim,n1,q1,n2,q2,n3,q3,n4,q4,n5,q5>
  ::quadraturePoint(unsigned int np, double coord[dim]) const
  {
      coord=abscissas[np];
  }




struct QuadratureTest
{
   double afunction(double x,double y) {return x*y;}
};

double bfunction(double x,double y) {return x*y;}

double cfunction(double x) {return x;}

}




