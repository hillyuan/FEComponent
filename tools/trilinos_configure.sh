#!/bin/bash

# Set this to the root of your Trilinos source directory.
#TRILINOS_PATH= "../Trilinos"

#
# You can invoke this shell script with additional command-line
# arguments.  They will be passed directly to CMake.
#
EXTRA_ARGS=$@

#
# Each invocation of CMake caches the values of build options in a
# CMakeCache.txt file.  If you run CMake again without deleting the
# CMakeCache.txt file, CMake won't notice any build options that have
# changed, because it found their original values in the cache file.
# Deleting the CMakeCache.txt file before invoking CMake will insure
# that CMake learns about any build options you may have changed.
# Experience will teach you when you may omit this step.
#
rm -f CMakeCache.txt

#
# Enable all primary stable Trilinos packages.
#
cmake .. -G "MSYS Makefiles" \
  -DCMAKE_INSTALL_PREFIX:FILEPATH="/usr/local/Trilinos" \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_FLAGS="-fpermissive" \
  -DCMAKE_FORTRAN_COMPILER=gfortran \
  -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python2 \
  -DGIT_EXEC:FILEPATH=/usr/bin/git \
  -DBLAS_LIBRARY_DIRS="/mingw64/lib" \
  -DTPL_BLAS_LIBRARIES="/mingw64/libopenblas.a" \
  -DTPL_LAPACK_LIBRARIES="/mingw64/liblapack.a" \
  -DBoost_INCLUDE_DIRS="C:\boost_1_59_0" \
  -DCMAKE_BUILD_TYPE:STRING=RELEASE \
  -DTPL_ENABLE_MPI:BOOL=ON \
  -DTrilinos_ENABLE_Fortran=ON \
  -DTrilinos_ENABLE_SEACAS=OFF \
  -DTrilinos_ENABLE_Kokkos=OFF \
  -DTrilinos_ENABLE_Anasazi=ON \
  -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
$EXTRA_ARGS \
#$TRILINOS_PATH