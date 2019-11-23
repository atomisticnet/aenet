#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#      Makefile for the unit tests included in the AENET package
#-----------------------------------------------------------------------
#+ This file is part of the AENET package.
#+
#+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
#+
#+ This Source Code Form is subject to the terms of the Mozilla Public
#+ License, v. 2.0. If a copy of the MPL was not distributed with this
#+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
#+
#+ This program is distributed in the hope that it will be useful, but
#+ WITHOUT ANY WARRANTY; without even the implied warranty of
#+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#+ Mozilla Public License, v. 2.0, for more details.
#+ ---------------------------------------------------------------------
#+ If you make use of AENET for your publication, please cite:
#+ [1] N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
#+ [2] J. Behler and M. Parrinello, Phys. Rev. Lett. 98 (2007) 146401.
#+
#+ If you used the Chebyshev descriptor, please cite:
#+ [3] N. Artrith, A. Urban, and G. Ceder, PRB 96 (2017) 014112.
#-----------------------------------------------------------------------
# 2015-05-26 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

TST      = ./tests/
TSTBIN   = ../unittests/

TEST1    = $(TSTBIN)test_xsflib.x
TEST2    = $(TSTBIN)test_lclist.x
TEST3    = $(TSTBIN)test_geometry.x
TEST4    = $(TSTBIN)test_feedforward.x
TEST5    = $(TSTBIN)test_symmfunc.x
TEST6    = $(TSTBIN)test_sfsetup.x
TEST7    = $(TSTBIN)test_sfbasis.x
TESTS    = $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5) $(TEST6) $(TEST7)

#------------------------------- rules --------------------------------#

.PHONY: test tests

test  : tests
tests : $(TESTS)
	@for t in $(TESTS); do if [ -x $$t ]; then $$t; fi; done

%.o : $(TST)%.f90 unittest.o
	$(FC) $(FCFLAGS) $< -o $*.o

$(TSTBIN)%.x : %.o $(OBJECTS) unittest.o Makefile
	$(LD) $(LDFLAGS) -o $@ $< $(OBJECTS) unittest.o $(LIB) $(NUMLIB)
