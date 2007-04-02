PERVER = 1.01
PERLABELED = Percolator v $(PERVER)
PERUNLABELED = Percolator unofficial version
BUILDDIR = Release

all: $(BUILDDIR)/percolator


$(BUILDDIR)/percolator: *.cpp *.h
	$(MAKE) -C $(BUILDDIR)

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h
