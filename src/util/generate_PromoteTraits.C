/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
/*

  This program generates the file, PromoteTraits.h

  This header file implements a type promotion template similar to
  that of Todd Veldhuizen, http://www.oonumerics.org/blitz/traits.html.
  However, by using some template template arguments, the template can
  be applied to other templates, like std::complex

  Willem van Straten 21 December 2004

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include <stdio.h>

using namespace std;

int main (int argc, char** argv)
{
  bool verbose = false;

  if (argc > 1 && string(argv[1]) == "-v")
    verbose = true;

  // it is assumed that header_PromoteTraits.h is copied to PromoteTraits.h
  // before this program is run
  ofstream out ("PromoteTraits.h", ios::app);
  
  if (!out) {
    cerr << "generate_PromoteTraits: could not open PromoteTraits.h" << endl;
    perror ("");
    return -1;
  }

  vector<string> types;

  // types are added in order of increasing precedence
  types.push_back( "float" );
  types.push_back( "double" );
  types.push_back( "long double" );
  types.push_back( "MEAL::ScalarMath" );

  for (unsigned i=0; i<types.size(); i++)
    for (unsigned j=0; j<types.size(); j++)

      out << "template<>\n"
        "class PromoteTraits< " << types[i] << ", " << types[j] << " > {\n"
	"public:\n"
	"  typedef " << types[max(i,j)] << " promote_type;\n"
	"};\n" << endl;

#if !defined(HAVE_BEST_PARTIAL_SPECIALIZATION) | !defined(HAVE_DEFAULT_PARTIAL_SPECIALIZATION)
  out << "#define PROMOTE_TRAITS_SPECIALIZE 1\n" << endl;
#endif

#if !defined(HAVE_COMPLEX_TEMPLATES)
  out << "#define COMPLEX_SPECIALIZE 1\n" << endl;
#endif

  out << "#endif // !__epsic_PromoteTraits_h\n" << endl;

  out.close ();

  return 0;

}
