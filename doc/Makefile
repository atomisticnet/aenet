EMACS = emacs
PANDOC = pandoc --smart --standalone
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
md : aenet-doc.md

%.pdf : %.org $(INFILES:.txt=.org)
	$(EMACS) $< --batch -f org-latex-export-to-pdf --kill

%.txt : %.org
	$(EMACS) $< --batch -f org-ascii-export-to-ascii --kill

%.md : %.org
	$(EMACS) $< --batch -f org-org-export-to-org --kill
	$(PANDOC) --from org --to markdown -o $@ $<.org
	rm -f $<.org

$(INFILEDIR)/%.txt : $(INFILEDIR)/%.org
	$(EMACS) $< --batch -f org-ascii-export-to-ascii --kill

clean:
	@rm -f aenet-doc.tex *~ $(INFILEDIR)/*~

veryclean: clean
	@rm -f aenet-doc.pdf aenet-doc.txt $(INFILES)
