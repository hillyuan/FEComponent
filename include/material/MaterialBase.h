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

#include <vector>
#include <array>
#include <memory>
#include "misc/registry.h"
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include "misc/VectTree.h"

namespace GNCLib
{
	enum StainType {Inifinite=0, GreenLagrangian, Log};
	typedef Eigen::Matrix<double,6,6> Constitutive_type;
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
          static unsigned const n_parameters = 0;      // number material constants
          static unsigned const n_status = 0;          // number status descript needed elastic/plastic e.g.
          static unsigned const n_status_parameters=0; // intermediate variables needed. e.g 6 plastic strain components
	   
	   protected:
	      std::string  name;
	   
	   public:
          virtual Constitutive_type ConstitutiveMatrix()=0;
   };
	
   /* Material can be constructed by elastic, plastic, visco etc*/
   class Material: public MaterialBase
   {
		std::vector<MaterialBase*>  materials;
		void SetMaterialConstants();
		void SetMaterialConstants(std::vector<double>&);
   };
	
   class Elastic: public MaterialBase
   {
 	    protected:
			XYLIB::CVectorTree material_constants;
			
		public:
			Elastic() {};
   };
  // REGISTER_SUBCLASS(MaterialBase, Elastic);
   
 
   class IsotropicElastic: public Elastic 
   {
	   typedef std::array<double, n_parameters>  constants_type;
    //   typedef std::function< constants_type (std::vector<double>)> material_constants;
	   // how to define the fucntion of material_constant, table, function?
	   public:
          static unsigned const n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned const n_status = 0;          // elastic only
          static unsigned const n_status_parameters=0; // no state changs ?
	   
	   public:
	      IsotropicElastic() {};
		  
          Constitutive_type ConstitutiveMatrix();
		  void SetMaterialConstants();
		  void SetMaterialConstants(std::vector<double>&);
		  
	   private:
		  double E,nu;      // temp var, get from mc in depdends given
   };
   REGISTER_SUBCLASS(Elastic, IsotropicElastic);
   
   class OrthotropicElastic: public Elastic
   {
	   typedef std::array<double, 9>  constants_type;
     //  typedef std::function< constants_type (std::vector<double>)> material_constants;
	   
	   public:
          static unsigned const n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned const n_status = 0;          // elastic only
          static unsigned const n_status_parameters=0; // no state changs ?
	   
	   public:
          Constitutive_type ConstitutiveMatrix();
		  
       private:
	  //    material_constants mc;
   };
   REGISTER_SUBCLASS(Elastic, OrthotropicElastic);
   
   class MooneyRivlin: public Elastic
   {
	   typedef std::array<double, 3>  constants_type;               // C10, C01, D1
     //  typedef std::function< constants_type (std::vector<double>)> material_constants;
	   
	   public:
          static unsigned const n_parameters = 3;      
          static unsigned const n_status = 0;          
          static unsigned const n_status_parameters=0; 
	   
	   public:
          Constitutive_type ConstitutiveMatrix();
		  
       private:
	 //     material_constants mc;
   };
   REGISTER_SUBCLASS(Elastic, MooneyRivlin);
   
   struct YieldFunction
   {
	   typedef Eigen::Tensor<double, 2> tensor_type;
       typedef std::array<double,6> array_type;
	   
	   virtual array_type NormalDirection(array_type& stress)=0;      // strain for strain type function
	   virtual tensor_type NormalDirection(tensor_type& stress)=0;
	   
	   virtual double FunctionVal(array_type& stress)=0;      // strain for strain type function
	   virtual double FunctionVal(tensor_type& stress)=0;
   };
   
   class Mises: public YieldFunction
   {
	   typedef double  constants_type;                   // initial yielding stress
    //   typedef std::function< constants_type (std::vector<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	 //     material_constants mc;
   };
   REGISTER_SUBCLASS(YieldFunction, Mises);
   
   class DruckerPrager: public YieldFunction
   {
	   typedef std::array<double, 2>  constants_type;      // c, fai
    //   typedef std::function< constants_type (std::vector<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	 //     material_constants mc;
   };
   REGISTER_SUBCLASS(YieldFunction, DruckerPrager);
   
   class MohrCoulomb: public YieldFunction
   {
	   typedef std::array<double, 2>  constants_type;      // c, fai
  //     typedef std::function< constants_type (std::vector<double>)> material_constants;
	   
	   array_type NormalDirection(array_type& stress);      // strain for strain type function
	   tensor_type NormalDirection(tensor_type& stress);
	   
	   double FunctionVal(array_type& stress);      // strain for strain type function
	   double FunctionVal(tensor_type& stress);
	   
	   private:
	 //     material_constants mc;
   };
   REGISTER_SUBCLASS(YieldFunction, MohrCoulomb);
   
   template<typename CYield>                 // kinematic type?
   class Plastic: public MaterialBase 
   {
	   typedef std::function< double (std::vector<double>)> yield;       //hardening rule here
	   typedef std::function< double (std::vector<double>)> kinematic;   //kinematic hardening rule here
	   
	   public:
          static unsigned const n_parameters = 2;      // Youngs Modulus, Poisson's ratio
          static unsigned const n_status = 1;          // 0: elastic/ 1: plastic
          static unsigned const n_status_parameters=7; // equivalent strain 1 + plastic strain 6
		  
	   public:
          Constitutive_type ConstitutiveMatrix();
		  
	   private:
		  std::shared_ptr<Elastic> Elastic_pt;
		  yield y;
		  kinematic k;
   };


} // End of namespace

#endif
