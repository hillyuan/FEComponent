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
        July 11,2005             YUAN Xi               Original
*/


#ifndef _quadraturebase_h_
#define _quadraturebase_h_


#include <string>
#include <map>
#include <stdexcept>

enum QuadratureScheme
{
    GaussLegendre,
    Hammer2,                  //2D Hammer--triangle
    Hammer3,                  //3D Hammer--tetgen
    Trapez,
    Simpson,
    Weddle,
    Milne,
    AssumedStrain_Shear
};


namespace GNCLib
{
    /*
      Constant holder of Abscissas and Weights for one-dimensional intergration.
         This class is used to build the Quadrature class, which complete the
      muilti dimensional integration and shouldn't be seen by the end user.
    */
    struct QuadratureBase
    {
        unsigned int ngauss;
        unsigned int coordDim;
        double*  abscissas;
        double*  weights;

        QuadratureBase() {};
        QuadratureBase(unsigned int n , unsigned int c)
        : ngauss(n), coordDim(c)
        {
           weights    = new double[ngauss];
           abscissas  = new double[ngauss*coordDim];
        }

        ~QuadratureBase()
        {
            if (abscissas != NULL) delete [] (abscissas);
            if (weights   != NULL) delete [] (weights);
        }
    };

    /*
        This factory is a Meyer's Singleton aims to generate various
        quadratureBase classes. New quadrature scheme can be regisited
        with its ID in string. Theretically, it maybe wrong in regisitation
        if same ID was used. Doesn't it should happen? Anyway, use long
        ID is suggested.
    */
    class QuadratureFactory
    {
    public:
        typedef QuadratureBase* (*CreateQuadratureCallback)(unsigned int);
    private:
        typedef std::map<QuadratureScheme, CreateQuadratureCallback> CallbackMap;

    public:
        static QuadratureFactory& Instance();

        bool Register(QuadratureScheme, CreateQuadratureCallback);
        bool Unregister(QuadratureScheme);
        QuadratureBase* CreateQuadrature(QuadratureScheme,unsigned int);
    private:
        QuadratureFactory() {};
        QuadratureFactory(const QuadratureFactory&) {};
        ~QuadratureFactory() {};

        CallbackMap callbacks_;
    };

    QuadratureFactory& QuadratureFactory::Instance()
    {
        static QuadratureFactory factory;
        return factory;
    }

    QuadratureBase* QuadratureFactory
    ::CreateQuadrature(QuadratureScheme QuadratureId, unsigned int n)
    {
        CallbackMap::const_iterator i=callbacks_.find(QuadratureId);
        if(i==callbacks_.end())
        {
            throw std::runtime_error("Unknown ID of Quadrature scheme");
        }
        return (i->second)(n);
    }

    bool QuadratureFactory::Register(QuadratureScheme quadratureId,
        CreateQuadratureCallback createFn)
    {
        return callbacks_.insert(
            CallbackMap::value_type(quadratureId,createFn)).second;
    }

    bool QuadratureFactory::Unregister(QuadratureScheme quadratureId)
    {
        return callbacks_.erase(quadratureId) == 1;
    }

} // End of namespace

#endif





