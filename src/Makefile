#-----------------------------------------------------------------------
#     This Makefile just prints a list of supported environments
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
# 2015-10-16 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

.PHONY: all
all: makefiles/README
	@cat makefiles/README

clean :
	rm -f *.o
	for f in *.mod; do rm -f $$f; done
	for f in *~; do rm -f $$f; done
