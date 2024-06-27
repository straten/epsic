
lib_LTLIBRARIES = libepsic.la

libepsic_la_SOURCES = mode.cpp sample.cpp \
	superposed.cpp composite.cpp disjoint.cpp coherent.cpp covariant.cpp \
	square_modulated_mode.cpp

pkginclude_HEADERS = mode.h modulated.h sample.h smoothed.h covariant.h

bin_PROGRAMS = epsic
epsic_SOURCES = epsic.cpp
LDADD = $(top_builddir)/libepsic.la @HEALPIX_LIBS@

libepsic_la_LIBADD = util/libutil.la 

# version-info should match AC_INIT version in configure.ac
libepsic_la_LDFLAGS = -version-info 1:0:1

AM_CPPFLAGS = -I$(top_srcdir)/util -I$(top_builddir)/util @HEALPIX_CFLAGS@
SUBDIRS = util

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool

ACLOCAL_AMFLAGS = -I config

clean-local:
	perl -pi -e 's|mode\.C|mode\.cpp|g' .deps/mode.Plo
	perl -pi -e 's|epsic\.C|epsic\.cpp|g' .deps/epsic.Po

include $(top_srcdir)/config/Makefile.local

