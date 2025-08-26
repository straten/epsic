dnl @synopsis EPSIC_COMPLEX_TEMPLATES
dnl
dnl Check whether std::complex<> templates work with user-defined types.
dnl define HAVE_COMPLEX_TEMPLATES
dnl
dnl @category Cxx
dnl @author Paul Demorest (demorest@gmail.com)
dnl @version 2023-08-18
dnl @license AllPermissive

AC_DEFUN([EPSIC_COMPLEX_TEMPLATES],
[AC_CACHE_CHECK(for complex templates accepting user-defined types,
ac_cv_cxx_complex_templates,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([

// template class X
template<typename T> class X { 
  public:
    T val; 
    X(T _val=0) { val=_val; }
    bool operator == (const X &b) const { return val==b.val; }
    // The "canonical" operator approach looks like this:
    //const X& operator*= (const X& b) { val*=b.val; return *this; }
    //friend const X operator * (X a, const X& b) { return a*=b; }
    // This one is a bit more compact:
    friend const X operator * (X a, const X& b) { return X(a.val*b.val); }
    friend const X operator / (X a, const X& b) { return X(a.val/b.val); }
    friend const X operator + (X a, const X& b) { return X(a.val+b.val); }
    friend const X operator - (X a, const X& b) { return X(a.val-b.val); }
};

#include <complex>

],
[
// Test if complex multiply on this type works
std::complex< X<double> > a, b, c;
c = a * b;
// Test if complex norm works
c = std::norm(c);
],
 ac_cv_cxx_complex_templates=yes,
 ac_cv_cxx_complex_templates=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_complex_templates" = yes; then
  AC_DEFINE(HAVE_COMPLEX_TEMPLATES,,
            [define if std::complex works with user-defined types])
fi
])

