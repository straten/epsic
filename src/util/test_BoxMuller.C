/***************************************************************************
 *
 *   Copyright (C) 2024 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "BoxMuller.h"

#include <algorithm>
#include <vector>

#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;

/*
 * This test passes when the code compiles, which verifies that
 * BoxMuller can be used as a generator function object
 */

int main ()
{
  unsigned npts = 10;
  vector<float> data (npts, 0.0);

  std::generate(data.begin(), data.end(), BoxMuller());

  std::ostringstream os;
  for (const auto& datum : data)
    os << datum << ' ';

  std::string expect = "-0.243308 -0.734373 -0.901952 -0.0282791 1.51791 -2.58562 -0.504347 0.626208 0.724619 0.114231 ";

  std::string got = os.str();
 
  if (got != expect)
  {
    std::cerr << "got='" << got << "'" << std::endl;
    std::cerr << "exp='" << expect << "'" << std::endl;

    std::cerr << "BoxMuller test FAIL" << std::endl;
    return -1;
  }

  std::cerr << "BoxMuller test PASS" << std::endl;
  return 0;
}

