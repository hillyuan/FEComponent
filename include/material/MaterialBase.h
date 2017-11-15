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


#include <array>
#include "factory.hpp"

namespace GNCLib
{
	enum StainType {Inifinite=0, GreenLagrangian, Log};
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
	      std::string :: name;
	      std::unordered_map<string, vector<double> > dictionary;   // temp, when reading?
	   
	   public:
          virtual void ConstitutiveMatrix()=0;
   };
   
   // A preprocessor define used by derived classes
	#define REGISTER_CLASS(NAME, TYPE) static Registrar registrar(NAME, [](void) -> MaterialBase * { return new TYPE();});
	
   class CompositeMaterial: public MaterialBase
   {
	   std::vector<MaterialBase*> :: materials;
   }
	
   class Elastic: public MaterialBase
   {};
   
 
   class IsotropicElastic: public Elastic 
   {
	   typedef std::array<double, 2>  constants_type;
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   // how to define the fucntion of material_constant, table, function?
	   public:
          static unsigned n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned n_status = 0;          // elastic only
          static unsigned n_status_parameters=0; // no state changs ?
	   
	   public:
		  IsotropicElastic() {
			  auto instance = BaseFactory<IsotropicElastic>:: Instance()->Create("ELASTIC");
			  REGISTER_CLASS("ELASTIC", IsotropicElastic);
		  };
          virtual void ConstitutiveMatrix();
		  
	   private:
	      material_constants mc;
   };
   
   class OrthotropicElastic: public Elastic
   {
	   typedef std::array<double, 9>  constants_type;
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   
	   public:
          static unsigned n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned n_status = 0;          // elastic only
          static unsigned n_status_parameters=0; // no state changs ?
	   
	   public:
		  OrthotropicElastic() {
			  auto instance = BaseFactory<OrthotropicElastic>:: Instance()->Create("ELASTIC");
			  REGISTER_CLASS("ELASTIC", OrthotropicElastic);
		  };
          virtual void ConstitutiveMatrix();
		  
       private:
	      material_constants mc;
   };
   
   class MonneyRivlin: public Elastic
   {
	   typedef std::array<double, 3>  constants_type;               // C10, C01, D1
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   
	   public:
          static unsigned n_parameters = 3;      
          static unsigned n_status = 0;          
          static unsigned n_status_parameters=0; 
	   
	   public:
		  OrthotropicElastic() {
			  auto instance = BaseFactory<OrthotropicElastic>:: Instance()->Create("ELASTIC");
			  REGISTER_CLASS("ELASTIC", OrthotropicElastic);
		  };
          virtual void ConstitutiveMatrix();
		  
       private:
	      material_constants mc;
   };
   
   struct YieldFunction
   {
	   typedef eigen::TensorFixedSize<double, Size<3,3>> tensor_type;
       typedef std::array<double,6> array_type;
	   
	   virtual array_type NormalDirection(array_type& stress)=0;      // strain for strain type function
	   virtual tensor_type NormalDirection(tensor_type& stress)=0;
	   
	   virtual double FunctionVal(array_type& stress)=0;      // strain for strain type function
	   virtual double FunctionVal(tensor_type& stress)=0;
   };
   
   class Mises: public YieldFunction
   {
	   typedef double  constants_type;                   // initial yielding stress
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	      material_constants mc;
   };
   
   class DruckerPrager: public YieldFunction
   {
	   typedef std::array<double, 2>  constants_type;      // c, fai
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	      material_constants mc;
   };
   
   class MohrCoulomb: public YieldFunction
   {
	   typedef std::array<double, 2>  constants_type;      // c, fai
       typedef std::function< constants_type(std::vetcor<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	      material_constants mc;
   };
   
   template<typename CYield>                 // kinematic type?
   class Plastic: public MaterialBase 
   {
	   typedef std::function< double(std::vetcor<double>)> yield;       //hardening rule here
	   typedef std::function< double(std::vetcor<double>)> kinematic;   //kinematic hardening rule here
	   
	   public:
          static unsigned n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned n_status = 1;          // 0: elastic/ 1: plastic
          static unsigned n_status_parameters=7; // equivalent strain 1 + plastic strain 6
		  
	   public:
          virtual void ConstitutiveMatrix();
		  
	   private:
		  Elastic* Elastic_pt;
		  yield y;
		  kinematic k;
   };


} // End of namespace

#endif
