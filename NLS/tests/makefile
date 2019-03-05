# Application constants

NLS     = .
NLSOBJ  = ../obj
NLSAPP  = ../bin
NLSLIB  = ..

LINSYS = ../external_libs/linsys

# Compilation parameters
CC      = g++
POLICE  = -Wall
#DEBUG  = -g -D_MM_
DEBUG   = -O3
INCLUDE = -I. -I$(NLSLIB)/include -I$(LINSYS)/include
CFLAGS  = $(INCLUDE) -D_UNIX_ $(POLICE) $(DEBUG)
STDLFLAGS  := r

SYSLIBS =       \
	../lib/libnls.a	\
	-lm  


APPMAINMOD         =       \
	main	\
	modfunc_1d         \
	modfunc_2d         \
	modstp           \
	modst           \
	modbeam2        \
	modbeam3	\
	utl

APPOBJS    = $(APPMAINMOD:%=$(NLSOBJ)/%.o)

$(NLSAPP)/tnls      : $(APPOBJS)
			  $(CC) -o $@ $(APPOBJS) $(SYSLIBS)

$(NLSOBJ)/%.o   : $(NLS)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

# make clean
.PHONY: clean
clean:
	rm -f $(LIBOBJS) $(APPOBJS)
