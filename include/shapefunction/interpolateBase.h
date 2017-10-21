/*
    Purpose: Build classes of 1D numerical integration and , as an except,
        2D Hammer integration. Those classes are used to generate quadrature
        classes used by the user and as conseqence, transprent to users.
           This is a part of Library of Generic Numerical Components *GNCLib*.
           You need only add your new derived class of QuadratureBase to
        accomplish different quadrature scheme. ID of each quadrature scheme
        is randam number generated automatically. Dangerous in generating same
        ID then?


        Date                    Programmer              Status
        ===========            ============          =============
        Oct. 06,2005             YUAN Xi               Original
*/


#ifndef _interpolatebase_h_
#define _interpolatebase_h_

enum InterpolationScheme
{
    Langrange,
    Hermit,
    Serendip,
    Triangle,
    CrossedTriangle,
    Tetrahedron,
    TrianglarSerendip,
    TrianglarLangerange,
    Hiearchary
};


namespace GNCLib
{
    /*
      Type-safy declearation of interpolation scheme. This templated class is
      virtual just to defining common interface of substantilized class. The
      interface provides shape function, shape derivative, and, if needed
      second shape derivative.
      This class is virtual one because we do not provide details of its
      member function. Materialization of its member function are realized
      in its specified class of this templated class . Although it can also
      be realized by complete specification of those memeber functions . But
      we here introduced specified class instead. Because it is less compiler
      -depedent and we can define those specified  class into a Singleton to
      save the computer resource as well.
      Number of interpolation points and dimensionality of interpolation are
      defined as static constants so as to strictly confine the dimension of
      arrays that its member function return.
    */
    template<InterpolationScheme, unsigned int, unsigned int, unsigned int>
    struct Interpolation
    {
        static const unsigned int dimension = 0;   //virtual value
        static const unsigned int numpoints = 0;

        void ShapeFunction(double[dimension],double[numpoints]);
        void ShapeDerivative(double[dimension],double[numpoints][dimension]);

        template<class array>
        array ShapeFunction(const array&);
        template<class array, class matrix>
        matrix ShapeDerivative(const array&);
    };

    /*
      Values of shape functions and shape derivatives in some specified points,
      which are generally defined by quadrature strategy,are stored in this table.
      Those values can be fetched in numerical caculation with no-calling of
      above shape functions to save computation time. It is much a optimized
      way considering the fact that the memeory used by this table is generally
      disregardably small. But when considering solid physical problem, following
      one would better.
      It is a const singleton. Template classes used to define it is supposed
      Singleton with memeber function Instance() also?
    */
    template<class CInterpolation, class CQuadrature>
    class InterpolationTable
    {
      private:
        unsigned int numgauss, numintps, dimension;
        double** TShapeFunction;
        double*** TShapeDerivative;

      private:
        InterpolationTable();
        ~InterpolationTable();

      public:
        static InterpolationTable& Instance();

        void ShapeFunction(unsigned int i, double[]) const;
        void ShapeDerivative(unsigned int i, double[][CInterpolation::dimension]) const;
    };

} // End of namespace

#endif
