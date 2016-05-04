# Makefile to save the result from the SRP009615 example
# Modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ggbio/vignettes/Makefile

HTMLS= SRP009615-results.html recount-quickstart.html

all: $(HTMLS)

clean:
	cp fakeVignettes/*.Rmd .

%.html: %.Rmd
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(rmarkdown); render('$*.Rmd')"
