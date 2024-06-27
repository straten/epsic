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

  std::generate(data.begin(), data.end(), BoxMuller(13));

  std::ostringstream os;
  for (const auto& datum : data)
    os << datum << ' ';

  std::string expect = "-0.0190504 0.892552 0.225053 0.484419 -1.50899 1.84724 -1.08279 -0.145933 -0.567948 0.78989 ";
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

