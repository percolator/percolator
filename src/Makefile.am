distdir = $(PACKAGE)_$(VERSION)_src

PERLABELED = Percolator v $(VERSION)
PERUNLABELED = Percolator unofficial version

AM_CFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP
AM_CXXFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP
AM_CPPFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP

perc_files = \
  Caller.cpp \
  Caller.h \
  DataSet.cpp \
  DataSet.h \
  DescriptionOfCorrect.cpp \
  DescriptionOfCorrect.h \
  FeatureNames.cpp \
  FeatureNames.h \
  MassHandler.cpp \
  MassHandler.h \
  Normalizer.cpp \
  Normalizer.h \
  Option.cpp \
  Option.h \
  PercolatorCInterface.cpp \
  PercolatorCInterface.h \
  PSMDescription.cpp \
  PSMDescription.h \
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
  svm.h \
  svm.cpp \
  ssl.cpp \
  ssl.h
 

toolkit_files = \
  MSToolkit/MSConvertFile.cpp \
  MSToolkit/MSObject.cpp \
  MSToolkit/MSObject.h \
  MSToolkit/MSReader.cpp \
  MSToolkit/MSReader.h \
  MSToolkit/MSSingleScan.cpp \
  MSToolkit/MSTDemo.cpp \
  MSToolkit/MSToolkitTypes.h \
  MSToolkit/RAMP/ramp_base64.cpp \
  MSToolkit/RAMP/ramp_base64.h \
  MSToolkit/RAMP/ramp.cpp \
  MSToolkit/RAMP/ramp.h \
  MSToolkit/Spectrum.cpp \
  MSToolkit/Spectrum.h \
  MSToolkit/zLib/adler32.c \
  MSToolkit/zLib/compress.c \
  MSToolkit/zLib/crc32.c \
  MSToolkit/zLib/crc32.h \
  MSToolkit/zLib/deflate.c \
  MSToolkit/zLib/deflate.h \
  MSToolkit/zLib/gzio.c \
  MSToolkit/zLib/infback.c \
  MSToolkit/zLib/inffast.c \
  MSToolkit/zLib/inffast.h \
  MSToolkit/zLib/inffixed.h \
  MSToolkit/zLib/inflate.c \
  MSToolkit/zLib/inflate.h \
  MSToolkit/zLib/inftrees.c \
  MSToolkit/zLib/inftrees.h \
  MSToolkit/zLib/trees.c \
  MSToolkit/zLib/trees.h \
  MSToolkit/zLib/uncompr.c \
  MSToolkit/zLib/zconf.h \
  MSToolkit/zLib/zlib.h \
  MSToolkit/zLib/zutil.c \
  MSToolkit/zLib/zutil.h \
  MSToolkit/SQLite/sqlite3.c \
  MSToolkit/SQLite/sqlite3ext.h \
  MSToolkit/SQLite/sqlite3.h

EXTRA_DIST = \
  ArrayLibrary/Array.cpp \
  ArrayLibrary/Arithmetic.cpp \
  ArrayLibrary/OrderedArray.cpp


if PERCOLATOR
nodist_libpercolator_la_SOURCES = $(perc_files)

nodist_libmstoolkit_la_SOURCES = $(toolkit_files)

libmstoolkit_la_CPPFLAGS = ${CPPFLAGS} -I./MSToolkit -I./MSToolkit/zLib -I./MSToolkit/RAMP -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -w -pthread

nodist_percolator_SOURCES = main.cpp
percolator_LDADD = libpercolator.la libqvality.la libArrayLibrary.la libmstoolkit.la -ldl -lpthread

PERCOLATOR_PROGRAM = percolator
PERCOLATOR_LIBRARY = libpercolator.la libmstoolkit.la

EXTRA_DIST +=  $(perc_files) $(toolkit_files) main.cpp

endif

libqvality_la_SOURCES = \
  BaseSpline.cpp \
  BaseSpline.h \
  Globals.cpp \
  Globals.h \
  LogisticRegression.cpp \
  LogisticRegression.h \
  Option.cpp \
  Option.h \
  PosteriorEstimator.cpp \
  PosteriorEstimator.h \
  Transform.h 



libArrayLibrary_la_SOURCES = \
  ArrayLibrary/ArrayTypes.h \
  ArrayLibrary/Array.h \
  ArrayLibrary/Arithmetic.h \
  ArrayLibrary/ArrayLibrary.h \
  ArrayLibrary/MatrixFunctions.cpp \
  ArrayLibrary/MatrixFunctions.h \
  ArrayLibrary/Numerical.cpp \
  ArrayLibrary/Numerical.h \
  ArrayLibrary/OrderedArray.h \
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


noinst_LTLIBRARIES = libqvality.la libArrayLibrary.la $(PERCOLATOR_LIBRARY)
bin_PROGRAMS = qvality $(PERCOLATOR_PROGRAM)

qvality_SOURCES = qvality_main.cpp
qvality_LDADD = libqvality.la libArrayLibrary.la

mark:
	echo $(PERUNLABELED) $(PERLABELED)
	perl -pi -e 's/$(PERUNLABELED)/$(PERLABELED)/g' *.cpp *.h

unmark:
	echo $(PERLABELED) $(PERUNLABELED)
	perl -pi -e 's/$(PERLABELED)/$(PERUNLABELED)/g' *.cpp *.h
