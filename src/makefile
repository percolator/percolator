PERVER = 1.01
PERLABELED = Percolator v $(PERVER)
PERUNLABELED = Percolator unofficial version
BUILDDIR = Release
LIBFILE = $(BUILDDIR)/libpercolator.a
LIBOBJS = $(BUILDDIR)/[^m]*.o

all: compile library

compile: $(BUILDDIR)/percolator

library: LIBFILE

$(BUILDDIR)/percolator: *.cpp *.h
	$(MAKE) -C $(BUILDDIR)

LIBFILE: compile
	${AR} r ${LIBFILE} ${LIBOBJS}

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h

