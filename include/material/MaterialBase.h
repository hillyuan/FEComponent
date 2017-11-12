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
        Oct. 12,2017             YUAN Xi               Original
*/


#ifndef _materialbase_h_
#define _materialbase_h_


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
   struct MaterialBase
   {
     public:
       virtual void ConstitutiveMatrix()=0;
   };


} // End of namespace

#endif
