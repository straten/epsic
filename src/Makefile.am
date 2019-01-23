
lib_LTLIBRARIES = libepsic.la
libepsic_la_SOURCES = mode.C
pkginclude_HEADERS = mode.h modulated.h sample.h smoothed.h

bin_PROGRAMS = epsic
epsic_SOURCES = epsic.C
LDADD = $(top_builddir)/libepsic.la

libepsic_la_LIBADD = util/libutil.la 
AM_CPPFLAGS = -I$(top_srcdir)/util
SUBDIRS = util


