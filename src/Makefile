#-----------------------------------------------------------------------
#     This Makefile just prints a list of supported environments
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
# 2015-10-16 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

.PHONY: all
all: makefiles/README
	@cat makefiles/README

clean :
	rm -f *.o
	for f in *.mod; do rm -f $$f; done
	for f in *~; do rm -f $$f; done
