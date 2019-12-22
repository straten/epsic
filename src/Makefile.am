
lib_LTLIBRARIES = libepsic.la

libepsic_la_SOURCES = mode.cpp sample.cpp \
	superposed.cpp composite.cpp disjoint.cpp coherent.cpp

pkginclude_HEADERS = mode.h modulated.h sample.h smoothed.h

bin_PROGRAMS = epsic
epsic_SOURCES = epsic.cpp
LDADD = $(top_builddir)/libepsic.la @HEALPIX_LIBS@

libepsic_la_LIBADD = util/libutil.la 
AM_CPPFLAGS = -I$(top_srcdir)/util -I$(top_builddir)/util @HEALPIX_CFLAGS@
SUBDIRS = util


