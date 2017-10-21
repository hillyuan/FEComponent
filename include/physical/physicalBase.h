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
        Oct. 08,2005             YUAN Xi               Original
*/


#ifndef _physicalbase_h_
#define _physicalbase_h_

#include "../shapefunction/interpolateBase.h"

enum PhysicalScheme
{
   //Mechanical problems
    Strain_Stress,
    Strain_Stress_AxisSymmetric,    //axis-symmetrical problem
    Strain_Stress_Shell,            //shell, including 1d beam
    Strain_Stress_Membrane,
    VelocityGradient_Stress,
    VelocityGradient_Stress_AxisSymmetric,    //axis-symmetrical problem
    VelocityGradient_Stress_Shell,            //shell, including 1d beam
    VelocityGradient_Stress_Membrane,

  //Thermo-conduct problem

   //electro-magnetic problems
    MagneticFluxDensity_MagneticField,
    ElectricFluxDemSity_ElectricField,
    MagneticFluxDensity_MagneticField_AxisSymmetric,
    ElectricFluxDemSity_ElectricField_AxisSymmetric
};


namespace GNCLib
{
  /*
     Values of shape functions and shape derivatives in some specified points,
     which are generally defined by quadrature strategy, are stored in this table.
     Those values can be fetched in numerical caculation with no-calling of
     above shape functions to save computation time. It is much a optimized
     way considering the fact that the memeory used by this table is generally
     disregardably small.
     It is a const singleton. Template classes used to define it is supposed
     Singleton with memeber function Instance() also?
   */
   template<PhysicalScheme, class CInterpolation, class CQuadrature>
   class PhysicalTable
   {
     private:
       unsigned int numgauss;
       double*** TPhysicalFunction;

     private:
       PhysicalTable();
       ~PhysicalTable();   //only for shell element

     public:
       static PhysicalTable& Instance();

       void PhysicalFunction(double***) const;
   };


} // End of namespace

#endif
