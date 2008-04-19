distdir = $(PACKAGE)_$(VERSION)

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
percolator_SOURCES = ../main.cpp
percolator_LDADD = libpercolator.la
percolator_CFLAGS = $(AM_CFLAGS)
