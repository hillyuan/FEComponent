/* ---------------------------------------------------------------------
 *
 *  File name: Interpolation.h
 *
 *  Author   : YUAN X.
 *
 *  Version  :
 *
 * ---------------------------------------------------------------------

   DESCRIPTION:

      Define the virtual basis of Interpolation Strategy

     During applicaions on numerical calculation, only one instance of
   one specific interpolation stategy is needed,ie., it is a Singlton. In
   our case, we need to consider the time of initiation of a strategy,
   because such strategys maybe quite a lot, list up all of them is noisy
   although it is quite small (with no data member). On the other
   hand, we ignore the destruction of those instances bacause of its
   ignorable smallness. We also take consider the life cycle of them
   because I cannot image any dependence between them and any other
   Singltons, as that pointed out by A. Andrei, would appear. In a word,
   we left the destruction of them at the end of application.
     S. Meyers stategy of Singlton design is adopted.Those singlton
   behaviors is implemented by th Loki libaray developded by A. Andrei.

   Revision History:

*/


#ifndef _felib_interpolation_h
#define _felib_interpolation_h

#include <functional>

enum ConvexType
{
    TRI3,
    REC4,
    HEX8
};



namespace GNCLib
{



                    /*
                       Define interpolation function for all type of fields.
                       Templated function objects (Functors) other than
                       templated functions are adopted due to the strong
                       suggestion of Scott Meyers:Effective STL:40, 46.

                       All fucntion are inlined. Consider the fact that all
                       objects here are Singlet (You must declear it!), no
                       much worry about the resource used even if it is quite
                       complex but high computation efficiency can be expected.
                    */

   template<int n1,int n2,int n3>
   struct SerendipShapeFunction
   : public std::unary_function<array, array>
   {
      array operator()(const array&) const;
   };

   template<int n1,int n2, int n3>
   array
   SerendipShapeFunction<n1,n2,n3>::operator()(const array& p) const
   {
    array v(3);

    int n=p.size();
    switch (n)
    {
    case 2:      // local coordinate
      v[0]=1.-p[0]-p[1];
      v[1]=p[0];
      v[2]=p[1];
      break;
    case 3:      // local area coordinate
      v = p;
      break;
    default:
      std::cout << "Error in definition of local ccordinate";
      break;
    }

    return v;
   }




    template<ConvexType,InterpolationScheme,class array>
    struct ShapeFunction
    : public std::unary_function<array, array>
    {
      inline
      array operator()(const array&) const;
    };

    template<CellType,InterpolationScheme,class array, class matrix>
    struct ShapeDerivative
    : public std::unary_function<matrix, array>
    {
      inline
      matrix operator()(const array&) const;
    };



  // --------Implementation-TRI3------------
  template<>
  array
  ShapeFunction<TRI3,Lagrange>::operator()(const array& p) const
  {
    array v(3);

    int n=p.size();
    switch (n)
    {
    case 2:      // local coordinate
      v[0]=1.-p[0]-p[1];
      v[1]=p[0];
      v[2]=p[1];
      break;
    case 3:      // local area coordinate
      v = p;
      break;
    default:
      std::cout << "Error in definition of local ccordinate";
      break;
    }

    return v;
  }


  template<>
  matrix
  ShapeDerivative<TRI3,Lagrange>::operator()(const array& p) const
  {
    matrix m(3,2);
    m(0,0)= 1.;  m(0,1)= 0.;
    m(1,0)= 0.;  m(1,1)= 1.;
    m(2,0)=-1.;  m(2,1)=-1.;
    return m;

  }


  // --------Implementation-REC4------------
  template<>
  array
  ShapeFunction<REC4,Lagrange>::operator()(const array& p) const
  {
    array v(4);
    v[0]=(1.-p[0])*(1.-p[1])/4.;
    v[1]=(1.+p[0])*(1.-p[1])/4.;
    v[2]=(1.+p[0])*(1.+p[1])/4.;
    v[3]=(1.-p[0])*(1.+p[1])/4.;
    return v;
  }

  template<>
  matrix
  ShapeDerivative<REC4,Lagrange>::operator()(const array& p) const
  {
    matrix m(4,2);
    m(0,0)=(-1.+p[1])/4.0;  m(0,1)=(-1.+p[0])/4.0;
    m(1,0)=( 1.-p[1])/4.0;  m(1,1)=(-1.-p[0])/4.0;
    m(2,0)=( 1.+p[1])/4.0;  m(2,1)=( 1.+p[0])/4.0;
    m(3,0)=(-1.-p[1])/4.0;  m(3,1)=( 1.-p[0])/4.0;
    return m;
  }

  // --------Implementation-HEX8------------
  template<>
  array
  ShapeFunction<HEX8,Lagrange>::operator()(const array& p) const
  {
    // natural coordinates of eight nodes
    // Be careful on its seqences
    matrix m(8,3);
    m(0,0)=-1.0;  m(0,1)=-1.0;  m(0,2)=-1.0;
    m(1,0)= 1.0;  m(1,1)=-1.0;  m(1,2)=-1.0;
    m(2,0)= 1.0;  m(2,1)= 1.0;  m(2,2)=-1.0;
    m(3,0)=-1.0;  m(3,2)= 1.0;  m(3,2)=-1.0;
    m(4,0)=-1.0;  m(4,1)=-1.0;  m(4,2)= 1.0;
    m(5,0)= 1.0;  m(5,1)=-1.0;  m(5,2)= 1.0;
    m(6,0)= 1.0;  m(6,1)= 1.0;  m(6,2)= 1.0;
    m(7,0)=-1.0;  m(7,2)= 1.0;  m(7,2)= 1.0;

    array v(8);
    for(int i=0;i<8;i++)
     v[i]=(1.+p[0]*m(i,0))*(1.+p[1]*m(i,1))*(1.+p[2]*m(i,2));

    return v;
  }

  template<>
  matrix
  ShapeDerivative<HEX8,Lagrange>::operator()(const array& p) const
  {
    // natural coordinates of eight nodes
    // Be careful on its seqences
    matrix m(8,3);
    m(0,0)=-1.0;  m(0,1)=-1.0;  m(0,2)=-1.0;
    m(1,0)= 1.0;  m(1,1)=-1.0;  m(1,2)=-1.0;
    m(2,0)= 1.0;  m(2,1)= 1.0;  m(2,2)=-1.0;
    m(3,0)=-1.0;  m(3,2)= 1.0;  m(3,2)=-1.0;
    m(4,0)=-1.0;  m(4,1)=-1.0;  m(4,2)= 1.0;
    m(5,0)= 1.0;  m(5,1)=-1.0;  m(5,2)= 1.0;
    m(6,0)= 1.0;  m(6,1)= 1.0;  m(6,2)= 1.0;
    m(7,0)=-1.0;  m(7,2)= 1.0;  m(7,2)= 1.0;

    matrix m1(8,3);
    for(int i=1;i<8;i++)
    {
      m1(i,0)=m(i,0)*(1.+p[1]*m(i,1))*(1.+p[2]*m(i,2));
      m1(i,1)=m(i,1)*(1.+p[0]*m(i,0))*(1.+p[2]*m(i,2));
      m1(i,2)=m(i,2)*(1.+p[0]*m(i,0))*(1.+p[1]*m(i,1));
    }
    return m1;
  }


} // namespace Interpolation

#endif
