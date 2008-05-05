if PERCOLATOR
distdir = $(PACKAGE)_$(VERSION)
else
distdir = qvality_$(VERSION)
endif
PERLABELED = Percolator v $(VERSION)
PERUNLABELED = Percolator unofficial version

CFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary
CXXFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary
CPPFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary

libqvality_la_SOURCES = \
  BaseSpline.cpp \
  BaseSpline.h \
  LogisticRegression.cpp \
  LogisticRegression.h \
  Option.cpp \
  Option.h \
  PosteriorEstimator.cpp \
  PosteriorEstimator.h \
  Transform.h 

libArrayLibrary_la_SOURCES = \
  ArrayLibrary/ArrayLibrary.h \
  ArrayLibrary/MatrixFunctions.cpp \
  ArrayLibrary/MatrixFunctions.h \
  ArrayLibrary/Numerical.cpp \
  ArrayLibrary/Numerical.h \
  ArrayLibrary/PackedMatrix.cpp \
  ArrayLibrary/PackedMatrix.h \
  ArrayLibrary/PackedSetFunctions.cpp \
  ArrayLibrary/PackedSetFunctions.h \
  ArrayLibrary/PackedVec.cpp \
  ArrayLibrary/PackedVec.h \
  ArrayLibrary/Random.cpp \
  ArrayLibrary/Random.h \
  ArrayLibrary/VectorFunctions.cpp \
  ArrayLibrary/VectorFunctions.h


if PERCOLATOR
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
  ssl.cpp \
  ssl.h

percolator_SOURCES = main.cpp
percolator_LDADD = libpercolator.la libqvality.la libArrayLibrary.la

PERCOLATOR_PROGRAM = percolator
PERCOLATOR_LIBRARY = libpercolator.la

endif

noinst_LTLIBRARIES = libqvality.la libArrayLibrary.la $(PERCOLATOR_LIBRARY)
bin_PROGRAMS = qvality $(PERCOLATOR_PROGRAM)
EXTRA_PROGRAMS = percolator

qvality_SOURCES = qvality_main.cpp
qvality_LDADD = libqvality.la libArrayLibrary.la

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h
