#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#         Makefile for tools included in the AENET package
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

TLS      = ./tools/
TLSBIN   = ../tools/

TOOL1    = $(TLSBIN)trnset2ASCII.x
TOOL2    = $(TLSBIN)neighbors.x
TOOL3    = $(TLSBIN)fingerprint.x
TOOLS    = $(TOOL1) $(TOOL2) $(TOOL3)

#------------------------------- rules --------------------------------#

.PHONY: tools

tools : $(TOOLS)

%.o : $(TLS)%.f90
	$(FC) $(FCFLAGS) $< -o $*.o

$(TLSBIN)%.x : %.o $(OBJECTS) Makefile
	$(LD) $(LDFLAGS) -o $@ $< $(OBJECTS) $(LIB) $(NUMLIB)
