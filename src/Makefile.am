
SUBDIRS                 = units third genutil fft

lib_LTLIBRARIES         = libpsrutil.la
libpsrutil_la_SOURCES  =
libpsrutil_la_LIBADD   = third/libthird.la units/libunits.la fft/libpsrfft.la \
			genutil/libgenutil.la tempo/libtempo.la @FLIBS@

if HAVE_CFITSIO
  SUBDIRS += fitsutil
  libpsrutil_la_LIBADD += fitsutil/libfitsutil.la
endif

SUBDIRS += tempo

if HAVE_TEMPO2
  SUBDIRS += tempo2
  libpsrutil_la_LIBADD += tempo2/libtempo_2.la
endif

if HAVE_QT
  SUBDIRS += qtutil
  libpsrutil_la_LIBADD += qtutil/libqtutil.la
endif

SUBDIRS += resources
libpsrutil_la_LIBADD += resources/libresources.la

include $(top_srcdir)/config/Makefile.local

