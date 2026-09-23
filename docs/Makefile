SPHINXBUILD ?= python -m sphinx
SPHINXOPTS ?= -n -W --keep-going
SOURCEDIR = .
BUILDDIR = _build

.PHONY: help html doctest linkcheck clean
help html doctest linkcheck clean:
	$(SPHINXBUILD) -M $@ $(SOURCEDIR) $(BUILDDIR) $(SPHINXOPTS)
