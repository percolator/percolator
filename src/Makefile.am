distdir = $(PACKAGE)_$(VERSION)
PERLABELED = Percolator v $(VERSION)
PERUNLABELED = Percolator unofficial version

#if DEBUG
#  AM_CPPFLAGS = -DDEBUG -g -O0
#  AM_CFLAGS = -DDEBUG -g -O0
#else !DEBUG
#  AM_CPPFLAGS = -O3
#  AM_CFLAGS = -O3
#endif !DEBUG

noinst_LTLIBRARIES = libpercolator.la
libpercolator_la_SOURCES = \
  Caller.cpp \
  Caller.h \
  DataSet.cpp \
  DataSet.h \
  Globals.cpp \
  Globals.h \
  IntraSetRelation.cpp \
  IntraSetRelation.h \
  Normalizer.cpp \
  Normalizer.h \
  Option.cpp \
  Option.h \
  PercolatorCInterface.cpp \
  PercolatorCInterface.h \
  ResultHolder.cpp \
  ResultHolder.h \
  SanityCheck.h \
  SanityCheck.cpp \
  Scores.cpp \
  Scores.h \
  SetHandler.cpp \
  SetHandler.h \
  SqtSanityCheck.h \
  SqtSanityCheck.cpp \
  StdvNormalizer.cpp \
  StdvNormalizer.h \
  UniNormalizer.cpp \
  UniNormalizer.h \
  f2c.h \
  gcvspl.c \
  gcvspl.h \
  ssl.cpp \
  ssl.h

bin_PROGRAMS = percolator
percolator_SOURCES = main.cpp
percolator_LDADD = libpercolator.la
#percolator_CFLAGS = $(AM_CFLAGS)

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h
