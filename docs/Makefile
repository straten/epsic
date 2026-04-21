# Makefile for PSRCHIVE Sphinx + Breathe documentation
#
# Builds Doxygen XML (via override configs), then runs Sphinx.
# Usage:  make html

SPHINXBUILD   = sphinx-build
SPHINXOPTS    =
BUILDDIR      = _build

# Doxygen override configs (build order matches Documentation/Makefile)
DOX_PROJECTS  = epsic

.PHONY: help clean doxygen sphinx html

help:
	@echo "  make html    – build Doxygen XML then Sphinx HTML"
	@echo "  make doxygen – build Doxygen XML only"
	@echo "  make clean   – remove all generated output"

doxygen:
	@for proj in $(DOX_PROJECTS); do \
		echo "==> Doxygen: $$proj"; \
		mkdir -p $(BUILDDIR)/$$proj; \
		doxygen $(BUILDDIR)/$$proj.dox.override \
			> $(BUILDDIR)/$$proj.out 2> $(BUILDDIR)/$$proj.err; \
		grep '^Error' $(BUILDDIR)/$$proj.err || true; \
	done

sphinx:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . $(BUILDDIR)/html

html: doxygen sphinx
	@echo "Build finished. Open $(BUILDDIR)/html/index.html"

clean:
	rm -rf $(BUILDDIR) $(BUILDDIR)/*/
