#include <iostream>
#include <stdlib.h>

//#ifndef  _felib_interpolation_h
//#include "Interpolation.h"
//#endif
#include "./quadrature/quadrature.h"
#include "./quadrature/qHammer.h"
#include "./quadrature/qGaussLegendre.h"
//#include "./shapefunction/interpolateBase.h"
//#include "./shapefunction/itpTriangle.h"
//#include "./shapefunction/itpSerendipShell.h"
//#include "./shapefunction/itpSerendip.h"
#include "./physical/physicalBase.h"
#include "./physical/strainFunction.h"
//#include "constitutive.h"
//#include "node.h"
//#include "element.h"
#ifndef _felib_init_h
// #include "init.h"
#endif



using namespace std;
using namespace GNCLib;

int main(int argc, char *argv[])
{

  QuadratureTest a;
  Quadrature< 3, 2,GaussLegendre,2,GaussLegendre,2,GaussLegendre > q;
  typedef double (QuadratureTest::*PFN)(double,double);
  PFN pfn=&QuadratureTest::afunction;
 // Loki::Functor<double,TYPELIST_2(double,double)> cmd(&a,&QuadratureTest::afunction);
//  boost::function<double (double,double) > pfn1;
//  pfn1=&bfunction;
  double aa=q.Integrate<double>( &bfunction );
//  aa=q.Integrate<double>( pfn1 );
//  double aa=q.Integrate<double>( bind(&bfunction, _1, 2) );
  QuadratureTest atest;
  double b=q.Integrate( atest, pfn );
//  b=q.Integrate<double>( mem_fun< double,QuadratureTest,Private::Coordinate<2> >(&QuadratureTest::afunction) );
  cout << "Sucess" << aa<< " "<<b << endl;
  return 0;
  /*array p(2,0.5),v(3);
  v=felib::InterpolationScheme::ShapeFunction<TRI3,Lagrange>()(p);

  std::vector<felib::Element*> eles(2);
  Element* ele=IsoparametricElement<TRI3, Lagrange, Deformation, Deformation::PLANESTRESS,
                                    Deformation::ELASTICITY,GaussFull>();

  std::cout << v[0]<< std::endl;
  std::cout<< v[1] << std::endl << v[2];
  std::cout << b << std::endl;

  system("PAUSE");*/

}

/*
class Thing {
public:
  void suspend(){std::cout<<"test";}
  int iii(int i){return i;}
};

typedef void (Thing::*PFN)();
typedef int (Thing::*PFN1)();

void integrate(PFN p)
{
  Thing a;
  (a.*p)();
}

int integ(Thing* a)
{
  PFN1 p=&Thing::iii;
  return (a->*p)(1);
}

main()
{
    PFN pfn = &Thing::suspend;
    integrate(pfn);
};
*/
