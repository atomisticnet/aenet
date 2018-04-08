#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#               Main Makefile for the AENET package
#-----------------------------------------------------------------------
#+ This file is part of the AENET package.
#+
#+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
#+
#+ This Source Code Form is subject to the terms of the Mozilla Public
#+ License, v. 2.0. If a copy of the MPL was not distributed with this
#+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
#+
#+ This program is distributed in the hope that it will be useful, but
#+ WITHOUT ANY WARRANTY; without even the implied warranty of
#+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#+ Mozilla Public License, v. 2.0, for more details.
#-----------------------------------------------------------------------
# 2015-05-26 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

#----------------------------- libraries ------------------------------#

LIB = ../lib/liblbfgsb.a

#------------------------------ objects -------------------------------#

OBJECTS  = aenet.o aeio.o chebyshev.o constants.o feedforward.o geometry.o \
           input.o io.o lclist.o optimize.o parallel.o potential.o random.o \
           sfbasis.o sfsetup.o sortlib.o symmfunc.o timing.o trainset.o xsflib.o

PICOBJECTS  = aenet_pic.o aeio_pic.o chebyshev_pic.o constants_pic.o \
              feedforward_pic.o geometry_pic.o input_pic.o io_pic.o \
              lclist_pic.o optimize_pic.o parallel_pic.o potential_pic.o \
              random_pic.o sfbasis_pic.o sfsetup_pic.o sortlib_pic.o \
              symmfunc_pic.o timing_pic.o trainset_pic.o xsflib_pic.o

#------------------------------ targets -------------------------------#

BIN      = ../bin/
EXT      = ./ext/

TARGET1  = $(BIN)generate.x
TARGET2  = $(BIN)train.x
TARGET3  = $(BIN)predict.x
TARGETS  = $(TARGET1) $(TARGET2) $(TARGET3)

#----------------------------- versioning -----------------------------#

version := $(shell cat VERSION)

#------------------------------- rules --------------------------------#

.SUFFIXES: .f90 .o .mod $(SUFFIXES)
.PHONY : veryclean clean tags all generate train predict lib

all      : $(TARGETS)
generate : $(TARGET1)
train    : $(TARGET2)
predict  : $(TARGET3)
lib      : libaenet.a libaenet.so

%.o : %.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o : %.F90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o : $(EXT)%.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%_pic.o : %.f90
	$(FC) $(FCFLAGS) -fPIC -o $@ $<

%_pic.o : %.F90
	$(FC) $(FCFLAGS) -fPIC -o $@ $<

%_pic.o : $(EXT)%.f90
	$(FC) $(FCFLAGS) -fPIC -o $@ $<

$(BIN)%.x : %.f90 $(OBJECTS) Makefile.inc
	$(LD) $(LDFLAGS) -o $@-$(version)-$(SUFFIX) $< $(OBJECTS) $(LIB) $(NUMLIB)

$(BIN)%.x : %.F90 $(OBJECTS) Makefile.inc
	$(LD) $(LDFLAGS) -o $@-$(version)-$(SUFFIX) $< $(OBJECTS) $(LIB) $(NUMLIB)

libaenet.a : $(OBJECTS) Makefile.inc
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

libaenet.so : $(PICOBJECTS) Makefile.inc
	$(LIBLD) $(LIBFLAGS) $(PICOBJECTS) $(LIBLIB) -o $@

veryclean: clean
	for f in $(TARGETS); do rm -f $$f-*; done
	if [ -e libaenet.a ]; then rm -f libaenet.a; fi
	if [ -e libaenet.so ]; then rm -f libaenet.so; fi
	if [ -e TAGS ]; then rm -f TAGS; fi
	for f in $(TESTS); do rm -f $$f-*; done

clean :
	for f in $(OBJECTS); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in $(PICOBJECTS); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in *.mod; do rm -f $$f; done
	for f in *~; do rm -f $$f; done

# requires exuberant-ctags
tags :
	rm -f TAGS
	ctags -e -R ../src ../lib

# unit tests (build with 'make tests')
include ./tests/Makefile.inc

# aenet tools (build with 'make tools')
include ./tools/Makefile.inc

#----------------------- explicit dependencies ------------------------#

aeio.o        : io.o
aenet.o       : aeio.o chebyshev.o geometry.o io.o lclist.o potential.o sfbasis.o sfsetup.o timing.o
geometry.o    : aeio.o constants.o io.o xsflib.o
input.o       : aeio.o geometry.o io.o
lclist.o      : sortlib.o
optimize.o    : feedforward.o io.o parallel.o random.o sortlib.o trainset.o
parallel.o    : aeio.o io.o input.o feedforward.o potential.o sfsetup.o trainset.o
potential.o   : aeio.o geometry.o io.o feedforward.o sfsetup.o trainset.o
sfbasis.o     : io.o chebyshev.o
sfsetup.o     : aeio.o io.o sfbasis.o symmfunc.o
trainset.o    : aeio.o io.o sfsetup.o
xsflib.o      : io.o
