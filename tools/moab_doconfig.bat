set SuiteSparse_DIR=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install
:: set MOAB_DIR= "C:/MyProgram/tools/moab/vs2019")
set TRILINOS_DIR=C:/MyProgram/solver/myTrilinos/vsbuild/install

del CMakeCache.txt
rd CMakeFiles

cmake .. -G "Visual Studio 16 2019" -Wno-dev   ^
	-DEIGEN3_DIR=C:/MyProgram/solver/eigen  ^
    -DCURL_INCLUDE_DIR=C:/MyProgram/tools/libcurl/include   ^
	-DCURL_LIBRARY_RELEASE=C:/MyProgram/tools/libcurl/lib    ^
	-DENABLE_BLASLAPACK=OFF  ^
	-DENABLE_FORTRAN:BOOL=OFF
::    -DBLAS_LIBRARIES=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install/lib64/lapack_blas_windows/libblas.lib   ^
::    -DLAPACK_LIBRARIES=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install/lib64/lapack_blas_windows/liblapack.lib ^
::  -DBLAS_LIBRARIES=C:/MyProgram/solver/OpenBLAS/mingw/lib/libopenblas.lib   
::  -DLAPACK_LIBRARIES=C:/MyProgram/solver/OpenBLAS/mingw/lib/libopenblas.lib
  

