#!/usr/bin/env python2

configure_options = [
  # Blas autodetec with cygwin blas at /usr/lib/liblapack,a,libblas.a
  '--with-mpi-include=/mingw64/include',
  '--with-mpi-lib=/mingw64/lib/libmsmpi.a',
  '--with-ar=/usr/bin/ar' ,
  '--with-shared-libraries=0',
  '--with-debugging=0',
  '--with-visibility=0',
  '--prefix=/usr/local/petsc',
  'FOPTFLAGS=-O3 -fno-range-check',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
