
noinst_LTLIBRARIES = libtrue_math.la

libtrue_math_la_SOURCES = true_math.c

pkginclude_HEADERS = true_math.h

TESTS = test_true_math

check_PROGRAMS = $(TESTS)

test_true_math_SOURCES = test_true_math.C

LDADD = libtrue_math.la

# disable all compiler optimizations, in particular -ffast-math
CFLAGS=

