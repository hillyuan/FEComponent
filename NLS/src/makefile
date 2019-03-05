# Aplication constants
NLS     = .
NLSOBJ  = ../obj
NLSLIB  = ../lib

CTRL    = $(NLS)/ctrl
UNISCH  = $(NLS)/ctrl/unisch
ORCTRL  = $(NLS)/ctrl/orctrl
MODEL   = $(NLS)/model
MATHPCK = $(NLS)/mathpack
MEMPCK  = $(NLS)/mempack
LINSYS  = $(NLS)/linsys

# Compilation parameters
CC      = g++
GC      = gcc
AR      = ar
RANLIB  = ranlib
POLICE  = -Wall
#DEBUG  = -g -D_MM_
DEBUG   = -O3
INCLUDE = -I../include
CFLAGS  = $(INCLUDE) -D_UNIX_ $(POLICE) $(DEBUG)
STDLFLAGS  := r

# Library modules

CTRLMOD         =       \
	ctrl

UNISCHMOD       =       \
	unisch          \
	dispctrl        \
	nrctrl          \
	workctrl        \
	mgdctrl        \
	alctrl          \
	srctrl          \
	strnctrl        \
	orctrl          \
	gdctrl

ORCTRLMOD       =       \
	oorctrl

MODELMOD        =       \
	model

MATHPCKMOD      =       \
	mathpack

MEMPCKMOD       =       \
	mempack

LINSYSMOD       =       \
	linearsystem	\
	crout	\
	pcgprofile	\
	pcgcsr \
	linsys  \
	sparsetools

LIBOBJS    = 	\
	  $(CTRLMOD:%=$(NLSOBJ)/%.o)    \
	  $(UNISCHMOD:%=$(NLSOBJ)/%.o)  \
	  $(ORCTRLMOD:%=$(NLSOBJ)/%.o)  \
	  $(MODELMOD:%=$(NLSOBJ)/%.o)   \
	  $(MEMPCKMOD:%=$(NLSOBJ)/%.o)  \
	  $(MATHPCKMOD:%=$(NLSOBJ)/%.o)	\
	  $(LINSYSMOD:%=$(NLSOBJ)/%.o)

$(NLSLIB)/libnls.a : $(LIBOBJS)
	$(AR) $(STDLFLAGS) $@ $(LIBOBJS)
	-$(RANLIB) $@

$(NLSOBJ)/%.o   : $(CTRL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(UNISCH)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(ORCTRL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(MODEL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(MATHPCK)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(MEMPCK)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(NLSOBJ)/%.o   : $(LINSYS)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<
		  
$(NLSOBJ)/%.o   : $(LINSYS)/%.c
		  $(GC) -c $(CFLAGS) -o $@ $<

# make clean
.PHONY: clean
clean:
	rm -f $(LIBOBJS) $(APPOBJS)
