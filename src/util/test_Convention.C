/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Convention.h"

#include <sstream>
#include <string>

using namespace std;

template<typename T>
void test_io (T value, string match)
{
  ostringstream ost;
  ost << value;

  if (ost.fail())
    throw string ("fail after operator <<");

  if (ost.str() != match)
    throw string ("operator << match failed");

  T output;

  stringstream iost;
  iost << match;

  if (iost.fail())
    throw string ("fail after operator << string!");

  iost >> output;

  if (iost.fail())
    throw string ("fail after operator >>");

  if (output != value)
    throw string ("operator >> match failed");
}

int main () try {
  
  cerr << "testing Basis Circular" << endl;
  test_io (Convention::Circular, "cir");

  cerr << "testing Basis Linear" << endl;
  test_io (Convention::Linear, "lin");
  
  cerr << "testing Basis Elliptical" << endl;
  test_io (Convention::Elliptical, "ell");

  cerr << "testing Hand Left" << endl;
  test_io (Convention::Left, "-1");

  cerr << "testing Hand Right" << endl;
  test_io (Convention::Right, "+1");

  cerr << "testing Argument Conventional" << endl;
  test_io (Convention::Conventional, "+1");

  cerr << "testing Argument Conjugate" << endl;
  test_io (Convention::Conjugate, "-1");

  cerr << "all tests passed" << endl;

  return 0;
}
catch (string& error) {
  cerr << error << endl;
  return -1;
}
