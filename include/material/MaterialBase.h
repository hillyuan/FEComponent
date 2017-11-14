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

#include "factory.hpp"

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
   class MaterialBase
   {
	   public:
          static unsigned n_parameters = 0;      // number material constants
          static unsigned n_status = 0;          // number status descript needed elastic/plastic e.g.
          static unsigned n_status_parameters=0; // intermediate variables needed. e.g 6 plastic strain components
	   
	   protected:
	      std::unordered_map<string, vector<double> > dictionary;
	   
	   public:
          virtual void ConstitutiveMatrix()=0;
   };
   
   // A preprocessor define used by derived classes
	#define REGISTER_CLASS(NAME, TYPE) static Registrar registrar(NAME, [](void) -> MaterialBase * { return new TYPE();});
   
 
   class IsoElastic: public MaterialBase 
   {
	   public:
          static unsigned n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned n_status = 0;          // elastic only
          static unsigned n_status_parameters=0; // no state changs ?
	   
	   public:
		  IsoElastic() {
			  auto instance = BaseFactory<IsoElastic>:: Instance()->Create("ELASTIC");
			  REGISTER_CLASS("ELASTIC", IsoElastic);
		  };
          virtual void ConstitutiveMatrix();
   };
   
   class MisesPlastic: public MaterialBase 
   {
	   public:
          static unsigned n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned n_status = 1;          // 0: elastic/ 1: plastic
          static unsigned n_status_parameters=7; // equivalent strain 1 + plastic strain 6
		  
	   public:
          virtual void ConstitutiveMatrix();
   };


} // End of namespace

#endif
