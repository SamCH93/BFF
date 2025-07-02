all: pdf code

# names of .Rnw .bib. and .bst files
FILE=bff
BIB=bibliography
BIBSTYLE=apalikedoiurl
.PHONY: all clean pdf pdf2 texdiff code

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib $(BIBSTYLE).bst
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## use texdiff to show difference to previous version
texdiff: tex old/$(FILE)-V1.tex old/ld.cfg
	-latexdiff -c old/ld.cfg --graphics-markup=none old/$(FILE)-V1.tex $(FILE).tex > $(FILE)-changes.tex
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-bibtex $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes

## extract R code
code: $(FILE).Rnw
	Rscript -e "knitr::purl('$(FILE).Rnw')" --vanilla

## clean
clean:  
	-rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).tex  $(FILE).bbl \
	$(FILE).out $(FILE).brf $(FILE).tex $(FILE).synctex.gz \
	$(FILE)-changes.aux  $(FILE)-changes.blg  $(FILE)-changes.log \
	$(FILE)-changes.tex  $(FILE)-changes.bbl $(FILE)-changes.out \
	$(FILE)-changes.brf $(FILE)-changes.synctex.gz
	-rm -r ./figure/
