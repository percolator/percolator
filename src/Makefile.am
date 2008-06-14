distdir = $(PACKAGE)_$(VERSION)_src

PERLABELED = Percolator v $(VERSION)
PERUNLABELED = Percolator unofficial version

AM_CFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP
AM_CXXFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP
AM_CPPFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary -I$(srcdir)/MSToolkit -I$(srcdir)/MSToolkit/zLib -I$(srcdir)/MSToolkit/RAMP

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

EXTRA_DIST = \
  ArrayLibrary/Array.cpp \
  ArrayLibrary/Arithmetic.cpp \
  ArrayLibrary/OrderedArray.cpp


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

if PERCOLATOR
libpercolator_la_SOURCES = \
  Caller.cpp \
  Caller.h \
  DataSet.cpp \
  DataSet.h \
  DescriptionOfCorrect.cpp \
  DescriptionOfCorrect.h \
  IntraSetRelation.cpp \
  IntraSetRelation.h \
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
  ssl.cpp \
  ssl.h

libmstoolkit_la_SOURCES = \
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
  MSToolkit/zLib/zutil.h

libmstoolkit_la_CPPFLAGS = ${CPPFLAGS} -I./MSToolkit -I./MSToolkit/zLib -I./MSToolkit/RAMP -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

percolator_SOURCES = main.cpp
percolator_LDADD = libpercolator.la libqvality.la libArrayLibrary.la libmstoolkit.la

PERCOLATOR_PROGRAM = percolator
PERCOLATOR_LIBRARY = libpercolator.la libmstoolkit.la

endif

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
