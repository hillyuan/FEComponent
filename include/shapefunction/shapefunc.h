/*
    Purpose: Classes of shape function


        Date                    Programmer              Status
        ===========            ============          =============
        Oct 04,2005             YUAN Xi               Original
*/


#ifndef _shapefunction_h_
#define _shapefunction_h_


#include <vector>
#include <functional>
#include <stdexcept>

enum InterpolationScheme
{
    Langrange,
    Hermit,
    Serendip,
    Triangle,
    Tetrahedron,
    TrianglarSerendip,
    TrianglarLangerange,
    Hiearchary
};


namespace GNCLib
{
  namespace Private
  {

  /*
     Default coordinate used in Qudrature class.
     Users defined coordinate can also be obtained from Quadrature class
    only if it with the same constructor. -See function quadraturePoint
    below.
  */
  template<unsigned int dim>
  struct Array
  {
                               /**
                               * Provide a way to get the
                               * dimension of an object without
                               * explicit knowledge of it's
                               * data type. Implementation is
                               * this way instead of providing
                               * a function <tt>dimension()</tt>
                               * because now it is possible to
                               * get the dimension at compile
                               * time without the expansion and
                               * preevaluation of an inlined
                               * function; the compiler may
                               * therefore produce more
                               * efficient code and you may use
                               * this value to declare other
                               * data types.
                               */
       static const unsigned int dimension = dim;
                                 /**
                                 * Store the values in a simple
                                 * array.  For <tt>dim==0</tt> store
                                 * one element, because otherways
                                 * the compiler would choke.  We
                                 * catch this case in the
                                 * constructor to disallow the
                                 * creation of such an object.
                                 */
       double values[(dim!=0) ? (dim) : 1];


                                 /**
                                 * Standard constructor. Creates
                                 * an origin.
                                 */
      Array() {};


                                /**
                                *  Constructor for one to five dimensional points.
                                *  This function is only implemented for <tt>dim==1
                                *  -5</tt>
                                */
      explicit Array(unsigned int x)
      {values[0]=x;}

      Array(const double x, const double y)
      {values[0]=x;   values[1]=y;}

      Array(const double x, const double y, const double z)
      {values[0]=x;   values[1]=y;  values[2]=z;}

      Array(const double x, const double y, const double z, const double w)
      {values[0]=x;   values[1]=y;  values[2]=z;   values[3]=w;}

      Array(const double x, const double y, const double z, const double w, const double u)
      {values[0]=x;   values[1]=y;  values[2]=z;   values[3]=w;  values[4]=u;}


                                  /**
                                  * Read access to the <tt>index</tt>th
                                  * Array.
                                  */
      double  operator () (const unsigned int index) const
      {
         Assert(index<dim);
         return values[index];
      }

      double  operator [] (const unsigned int index) const
      {
         Assert(index<dim);
         return values[index];
      }

      Array &operator= (const Array& obj)
      {
          for(register int i=0; i<dim; i++) {
              values[i]=obj.values[i];
          }
      }


   };


   }  // end of namespace Private

   std::unary_function< std::vector<double>, std::vector<double> > afunc;

    /*
      Shape function class are generated after interpolation scheme and order
    of each dimension are denoted. It is a singlton also.
    */
    template
    <
      InterpolationScheme q,
      int n1,
      int n2,
      int n3
    >
    class Interpolation
    {
      private:
        int npoints;

      public:
        Interpolation();

        template<class array>
        array  shapefunc( const array& ) const;
   //     std::unary_function< std::vector<double>, std::vector<double> > shapefunc;

        template<class array, class matrix>
        matrix shapederive(const array&) const;
    };


    template
    <>
    class Interpolation<Serendip,1,1,0>
    {
      private:
        int npoints;

      public:
        Interpolation():npoints(4) {};

        template<class array>
        array  shapefunc( const array& ) const;
   //     std::unary_function< std::vector<double>, std::vector<double> > shapefunc;

        template<class array, class matrix>
        matrix shapederive(const array&) const;
    };



    // --------Implementation------------
    template
    <>
    template<class array>
    array
    Interpolation<Serendip,1,1,0>::shapefunc( const array& p) const
    {
      array v(4);
      v[0]=(1.-p[0])*(1.-p[1])/4.;
      v[1]=(1.+p[0])*(1.-p[1])/4.;
      v[2]=(1.+p[0])*(1.+p[1])/4.;
      v[3]=(1.-p[0])*(1.+p[1])/4.;
      return v;
    }



}

#endif
