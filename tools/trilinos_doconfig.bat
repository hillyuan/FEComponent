set SuiteSparse_DIR=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install
set HDF5_ROOT=C:/MyProgram/tools/hdf5-1.10.5/build/_CPack_Packages/win64/ZIP/HDF5-1.10.5-win64

del CMakeCache.txt
rmdir /Q /s CMakeFiles

cmake .. ^
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON   ^
  -DTrilinos_ENABLE_EXAMPLES:BOOL=ON  ^
  -DTrilinos_ENABLE_TESTS:BOOL=ON  ^
  -DTrilinos_ENABLE_Epetra=ON  ^
  -DTrilinos_ENABLE_Kokkos=OFF ^
  -DTrilinos_ENABLE_Anasazi=ON ^
  -DTrilinos_ENABLE_Sundance=ON ^
  -DTPL_BLAS_LIBRARIES=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install/lib64/lapack_blas_windows/libblas.lib   ^
  -DTPL_LAPACK_LIBRARIES=C:/MyProgram/solver/suitesparse-metis-for-windows/build/install/lib64/lapack_blas_windows/liblapack.lib ^
  -DTPL_ENABLE_HDF5=ON  ^
  -DHDF5_LIBRARY_DIRS='C:/MyProgram/tools/HDF_Group/HDF5/1.10.5/lib'      ^
  -DHDF5_INCLUDE_DIRS='C:/MyProgram/tools/HDF_Group/HDF5/1.10.5/include'
  

