dnl @synopsis EPSIC_DEFAULT_PARTIAL_SPECIALIZATION
dnl
dnl If the compiler can match a partial specialization of a template 
dnl with template argument when the template argument has a default argument, 
dnl define HAVE_DEFAULT_PARTIAL_SPECIALIZATION.
dnl
dnl @category Cxx
dnl @author Willem van Straten <vanstraten.willem@gmail.com>
dnl @version 2005-02-12
dnl @license AllPermissive

AC_DEFUN([EPSIC_DEFAULT_PARTIAL_SPECIALIZATION],
[AC_CACHE_CHECK(for default template argument in partial specialization,
ac_cv_cxx_default_partial_specialization,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([

// template class X
template<typename A> class X { };

// template specializiation for template argument
template<typename A, template<class> class C>
class X< C<A> > 
{
  public: class Y {};
};

// a function that returns X<T>::Y, given an argument of type T
template<typename T>
typename X<T>::Y func (T& t)
{ 
  return typename X<T>::Y();
}

// a template class with a default argument
template<class T, class U=T> class Test2
{
};

],
[

// an instance of the Test2<> template class
Test2<int> mi; 

// ... can the compiler do this?
X< Test2<int> >::Y md = func (mi);

],
 ac_cv_cxx_default_partial_specialization=yes,
 ac_cv_cxx_default_partial_specialization=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_default_partial_specialization" = yes; then
  AC_DEFINE(HAVE_DEFAULT_PARTIAL_SPECIALIZATION,,
            [define if partial specialization accepts default template arg])
fi
])

