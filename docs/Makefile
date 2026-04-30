# Makefile for PSRCHIVE Sphinx + Breathe documentation
#
# Builds Doxygen XML (via override configs), then runs Sphinx.
# Usage:  make html

SPHINXBUILD   = sphinx-build
SPHINXOPTS    =
BUILDDIR      = _build
DOXDIR        = _doxygen

.PHONY: help clean epsic util breathe sphinx rtd

help:
	@echo "  make rtd     – build Doxygen XML then Sphinx HTML"
	@echo "  make xml     – build Doxygen XML only"
	@echo "  make html    – build Doxygen HTML only"
	@echo "  make clean   – remove all generated output"

Doxyfile:
	doxygen -g

epsic util: Doxyfile
	@mkdir -p $(DOXDIR)/xml/$@
	doxygen $@.dox > $(DOXDIR)/$@.dox.out 2> $(DOXDIR)/$@.dox.err
	@grep '^Error' $(DOXDIR)/$@.dox.err || echo 2> /dev/null

breathe: epsic util
	breathe-apidoc -q -o api/epsic -p epsic -m -f -g class $(DOXDIR)/xml/epsic
	breathe-apidoc -q -o api/util -p util -m -f -g class,struct $(DOXDIR)/xml/util
	cp overrides/api_epsic_classlist api/epsic/classlist.rst
	cp overrides/api_util_classlist api/util/classlist.rst
	cp overrides/api_util_structlist api/util/structlist.rst

sphinx: breathe
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . $(BUILDDIR)/html

rtd: sphinx
	@echo "Build finished. Open $(BUILDDIR)/html/index.html"

clean:
	rm -rf $(BUILDDIR) $(DOXDIR) api Doxyfile Doxyfile.bak
