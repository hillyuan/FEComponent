#!/bin/bash

export PETSC_DIR=/usr/local/petsc
export LIBMESH_DIR=/usr/local/libmesh

./configure  \
        --with-methods="opt"  \
		--prefix="$LIBMESH_DIR" \
		--enable-unique-ptr \
		--disable-warnings \
		--disable-netcdf \
		--disable-exodus  \
		--disable-metis \
		--disable-parmetis \
		--disable-fparser \
		--enable-petsc-required \
		
# add "enablepetsc=yes" at line 33062 of file 'configure'. Mingw builded petsc library not discerned!