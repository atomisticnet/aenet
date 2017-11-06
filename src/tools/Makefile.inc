#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#         Makefile for tools included in the AENET package
#-----------------------------------------------------------------------
#+ This file is part of the AENET package.
#+
#+ Copyright (C) 2012-2017 Nongnuch Artrith and Alexander Urban
#+
#+ This program is free software: you can redistribute it and/or modify
#+ it under the terms of the GNU General Public License as published by
#+ the Free Software Foundation, either version 3 of the License, or
#+ (at your option) any later version.
#+
#+ This program is distributed in the hope that it will be useful, but
#+ WITHOUT ANY WARRANTY; without even the implied warranty of
#+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#+ General Public License for more details.
#+
#+ You should have received a copy of the GNU General Public License
#+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
