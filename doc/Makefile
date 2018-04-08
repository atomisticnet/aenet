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
#+ This Source Code Form is subject to the terms of the Mozilla Public
#+ License, v. 2.0. If a copy of the MPL was not distributed with this
#+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
#+
#+ This program is distributed in the hope that it will be useful, but
#+ WITHOUT ANY WARRANTY; without even the implied warranty of
#+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#+ Mozilla Public License, v. 2.0, for more details.
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
