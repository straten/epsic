# Makefile for PSRCHIVE Sphinx + Breathe documentation
#
# Builds Doxygen XML (via override configs), then runs Sphinx.
# Usage:  make html

SPHINXBUILD   = sphinx-build
SPHINXOPTS    =
BUILDDIR      = _build
DOXDIR        = _doxygen

.PHONY: help clean epsic util apidocs docs

help:
	@echo "  make docs    – build Doxygen XML, Breathe RST, then Sphinx HTML"
	@echo "  make apidocs – build Doxygen XML, Breathe RST only"
	@echo "  make epsic   – build Doxygen XML for epsic core only"
	@echo "  make util    - build Doxygen XML for epsic util only"
	@echo "  make clean   – remove all generated output"

Doxyfile:
	doxygen -g

epsic util: Doxyfile
	@mkdir -p $(DOXDIR)/xml/$@
	doxygen $@.dox > $(DOXDIR)/$@.dox.out 2> $(DOXDIR)/$@.dox.err
	@grep '^Error' $(DOXDIR)/$@.dox.err || echo 2> /dev/null

apidocs: epsic util
	breathe-apidoc -q -o api/epsic -p epsic -m -f -g class $(DOXDIR)/xml/epsic
	breathe-apidoc -q -o api/util -p util -m -f -g class,struct $(DOXDIR)/xml/util
	cp overrides/api_epsic_classlist api/epsic/classlist.rst
	cp overrides/api_util_classlist api/util/classlist.rst
	cp overrides/api_util_structlist api/util/structlist.rst

docs: apidocs
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo "Build finished. Open $(BUILDDIR)/html/index.html"

clean:
	rm -rf $(BUILDDIR) $(DOXDIR) api Doxyfile Doxyfile.bak
