PERVER = 1.01
PERLABELED = Percolator v $(PERVER)
PERUNLABELED = Percolator unofficial version
BUILDDIR = Release
LIBFILE = $(LIBDIR)/libpercolator.a
LIBOBJS = $(BUILDDIR)/[^m]*.o
LIBDIR = $(BUILDDIR)/lib

all: compile library

compile: $(BUILDDIR)/percolator

library: LIBFILE

$(BUILDDIR)/percolator: *.cpp *.h
	$(MAKE) -C $(BUILDDIR)

LIBFILE: compile
	@mkdir -p ${LIBDIR}
	${AR} r ${LIBFILE} ${LIBOBJS}

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h

