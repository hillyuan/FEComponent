#! /bin/bash

PARMETIS_INSTALL_DIR=/usr/local/parmetis

rm -f CMakeCache.txt
rm -rf CMakeFiles

cmake .. -G "MSYS Makefiles" \
  -DCMAKE_INSTALL_PREFIX:PATH=${PARMETIS_INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS="-D_MINGW_" \
  -DCMAKE_C_FLAGS="-D_MINGW_" \
  -DMPI_C_LIBRARIES=/mingw64/lib/libmsmpi.a \
  -DMPI_CXX_LIBRARIES=/mingw64/lib/libmsmpi.a \
  -DGKLIB_PATH=/d/programs/parmetis-4.0.3/metis/GKlib \


