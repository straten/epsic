# Makefile for PSRCHIVE Sphinx + Breathe documentation
#
# Builds Doxygen XML (via override configs), then runs Sphinx.
# Usage:  make html

SPHINXBUILD   = sphinx-build
SPHINXOPTS    =
BUILDDIR      = _build
DOXDIR        = _doxygen

.PHONY: help clean html xml sphinx rtd

help:
	@echo "  make rtd     – build Doxygen XML then Sphinx HTML"
	@echo "  make xml     – build Doxygen XML only"
	@echo "  make html    – build Doxygen HTML only"
	@echo "  make clean   – remove all generated output"

html xml:
	@mkdir -p $(DOXDIR)/html
	@cp commands.sty $(DOXDIR)/html
	doxygen $@.dox > $(DOXDIR)/$@.dox.out 2> $(DOXDIR)/$@.dox.err
	@grep '^Error' $(DOXDIR)/$@.dox.err || echo 2> /dev/null

sphinx:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . $(BUILDDIR)/html

rtd: xml sphinx
	@echo "Build finished. Open $(BUILDDIR)/html/index.html"

clean:
	rm -rf $(BUILDDIR) $(DOXDIR) api
