#! /bin/bash

LIFEV_INSTALL_DIR=/usr/local/lifev

TRILINOS_HOME=/usr/local/trilinos


rm -f CMakeCache.txt
rm -rf CMakeFiles

cmake .. -G "MSYS Makefiles" \
  -DCMAKE_INSTALL_PREFIX:PATH=${LIFEV_INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS="-std=gnu++11 -fpermissive" \
  -DTrilinos_INCLUDE_DIRS:PATH=${TRILINOS_HOME}/include \
  -DTrilinos_LIBRARY_DIRS:PATH=${TRILINOS_HOME}/lib \
  -DLifeV_ENABLE_ALL_PACKAGES:BOOL=ON \
  -DTPL_METIS_INCLUDE_DIRS=/usr/local/parmetis/include \
  -DTPL_ParMETIS_INCLUDE_DIRS=/usr/local/parmetis/include \
  -DLifeV_ENABLE_TESTS:BOOL=ON \

