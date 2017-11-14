/*==============================================================================

                                    FEComponent

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 2013 Xi YUAN

   This file is part of FECompnent.

   FEComponent is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FEComponent is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with FEComponent. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                        Data structure

  ==============================================================================*/

#ifndef __TL_FACTORY_H
#define __TL_FACTORY_H

#include <string>
#include <unordered_map>
#include <function>


namespace XYLIB
{

    // register a factory function for BaseClass
	template< class BaseClass, class BaseFactory >
	struct Registrar {
		Registrar(std::string name, std::function<BaseClass*(void)> FactoryFunction) {
			BaseFactory::Instance()->RegisterFactoryFunction(name, FactoryFunction);
		}
	};

	// The singleton factory to generate BaseClass
	template< class BaseClass >
	struct BaseFactory
	{
		public:
			// Get the singleton instance of the factory
			static CFactory * Instance() {
				static CFactory factory;
				return &factory
			}

			// register a factory function to create an instance of className
			void RegisterFactoryFunction(std::string name, std::function<BaseClass*(void)> FactoryFunction){
				MapRegistry[name] = FactoryFunction;
			}

			// create an instance of a registered class
			shared_ptr<BaseClass> Create(std::string const& name) {
				BaseClass * instance = nullptr;

				// find name in the registry and call factory method.
				auto it = MapRegistry.find(name);
				if(it != MapRegistry.end())
				instance = it->second();
    
				// wrap instance in a shared ptr and return
				if(instance != nullptr)
					return std::shared_ptr<BaseClass>(instance);
				else
					return nullptr;
			}

		private:
			// a private ctor
			CFactory(){};

			// the registry of factory functions
			std::unordered_map<string, function<BaseClass*(void)>> MapRegistry;

	};
}

#endif
