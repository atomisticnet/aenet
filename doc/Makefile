#-----------------------------------------------------------------------
#  This Makefile generates the AENET documentation in several formats
#-----------------------------------------------------------------------
#  Generating the PDF and TXT documentation only requires Emacs.
#
#  $ make
#
#  Generating the Markdown version of the doc also requires Pandoc
#  (https://pandoc.org), as the vanilla org-mode-to-markdown export is
#  not giving satisfactory results. Along with the doc, also the
#  README.org file in the main directory is converted to GitHub-flavored
#  Markdown to be useful for the GitHub repo, as GitHub's org-mode parser
#  is quite buggy.
#
#  $ make md
#-----------------------------------------------------------------------
#+ This file is part of the AENET package.
#+
#+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
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
# 2018-04-07 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

EMACS = emacs
PANDOC = pandoc --standalone
INFILEDIR = input-files

INFILES = $(INFILEDIR)/aenet-predict.py-input.txt \
          $(INFILEDIR)/fingerprint-setup.txt \
          $(INFILEDIR)/generate.x-input.txt \
          $(INFILEDIR)/predict.x-input.txt \
          $(INFILEDIR)/train.x-input.txt

.PHONY: pdf txt clean

all : pdf txt
pdf : aenet-doc.pdf
txt : aenet-doc.txt $(INFILES)
md : aenet-doc.md ../README.md

%.pdf : %.org $(INFILES:.txt=.org)
	$(EMACS) $< --batch -f org-latex-export-to-pdf --kill

%.txt : %.org
	$(EMACS) $< --batch -f org-ascii-export-to-ascii --kill

%.md : %.org
	$(EMACS) $< --batch -f org-org-export-to-org --kill
	$(PANDOC) --from org+smart --to gfm -o $@ $<.org
	rm -f $<.org

$(INFILEDIR)/%.txt : $(INFILEDIR)/%.org
	$(EMACS) $< --batch -f org-ascii-export-to-ascii --kill

clean:
	@rm -f aenet-doc.tex *~ $(INFILEDIR)/*~

veryclean: clean
	@rm -f aenet-doc.pdf aenet-doc.txt $(INFILES)
