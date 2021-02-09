/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include <iostream>
#include <string>

#include "random.h"

template <class Type1, class Type2, class Element>
class MatrixTest {
  
 public:
  
  MatrixTest ()
    { verbose = false; scale=1e3; tolerance = 1e-12; random_init(); }

  virtual ~MatrixTest () {}
  
  virtual void runtest (unsigned nloops);
  
  bool verbose;

  double scale;

  double tolerance;

};

template <class T, class U>
T max (const T& a, const U& b)
{
  if (a>b)
    return a;
  else
    return b;
}

template <class Type1, class Type2, class Element>
void test_matrix (const Type1& t1_1, const Type1& t1_2, const Type2& t2_null,
		  const Element& p, double scale, double tolerance, bool verbose)
{
  Type1 temp1;
  Type1 temp2;

  Element one (1.0);

  //
  // test operator =
  //
  if (verbose)
  std::cerr << "\n*************** Testing operator =" << std::endl;

  if (verbose)
  std::cerr << "Type1 = Type1" << std::endl;
  temp1 = t1_1;

  if (verbose) std::cerr << "Type2 = Type1" << std::endl;
  Type2 t2_1 = temp1;

  if (verbose) std::cerr << "Type2 = Type2" << std::endl;
  Type2 t2_2 = t2_1;
 
  if (verbose) std::cerr << "Type1 = Type2" << std::endl;
  temp1 = t2_2;

  if (verbose) std::cerr << "Type1 == Type1" << std::endl;
  if (t1_1 != temp1)
    throw std::string ("test_matrix "
		 "Error casting from Type1 to Type2");

  if (verbose) std::cerr << "Type1 == Type2" << std::endl;
  if (t1_1 != t2_1)
    throw std::string ("test_matrix "
		 "Error casting from Type1 to Type2");

  //
  // test operator + and - (implicitly tests operator+= and -=)
  //
  if (verbose)
    std::cerr << "\n*************** Testing operator + and -" << std::endl;

  t2_2 = t1_2;

  if (verbose) std::cerr << "Type1 + Type2" << std::endl;
  temp1 = t1_1 + t2_2;

  t2_2 = temp1;

  if (verbose) std::cerr << "Type2 - Type1" << std::endl;
  t2_1 = t2_2 - t1_2;

  if (verbose) std::cerr << "Type1 == Type2" << std::endl;

  if (norm(t1_1 - t2_1)/norm(t2_1) > tolerance)
  {
    std::cerr << "t1_1 = " << t1_1 << " != "
	 << "t2_1 = " << t2_1 << std::endl;

    throw std::string ("test_matrix "
		  "Error adding and subtracting Type1 and Type2");
  }


  //
  // test operator * and / Element (implicitly tests operator*= and /=)
  //
  if (verbose)
    std::cerr << "\n*************** Testing operator * and / Element" 
              << std::endl;

  Element dfactor;
  random_value (dfactor, scale);

  t2_1 = t1_1;

  if (verbose) std::cerr << "Element * Type2" << std::endl;
  t2_2 = dfactor * t2_1;

  if (verbose) std::cerr << "Type2 / Element" << std::endl;
  t2_1 = t2_2 / dfactor;

  if (verbose) std::cerr << "Type1 == Type2" << std::endl;
  if (norm(t2_1 - t1_1)/norm(t2_1) > tolerance)
  {
    std::cerr << "Error Type2=" << t2_1 
	 << " * / Element=" << dfactor
	 << " != Type1=" << t1_1 << std::endl;
    throw std::string ("test_matrix Error in Type2 * / Element");
  }

  //
  // test operator Type1 * Type1 (implicitly tests operator*=)
  //
  if (verbose) std::cerr << "\n*************** Testing operator * Type1" << std::endl;

  if (verbose) std::cerr << "Type1::identity=" << Type1::identity() << std::endl;

  if (verbose) std::cerr << "Type1 * Type1::identity" << std::endl;
  t2_1 = t1_1 * Type1::identity();

  if (verbose) std::cerr << "Type1 == Type2" << std::endl;
  if (t1_1 != t2_1)
    throw std::string ("test_matrix Error Type1 * Type1::identity");


  //
  // test complex conjugate
  //
  if (verbose) std::cerr << "\n*************** Testing conjugate" << std::endl;

  if (verbose) std::cerr << "conj (t1_1 * t1_2)" << std::endl;
  t2_1 = conj (t1_1 * t1_2);

  if (verbose) std::cerr << "conj (t1_1) * conj (t1_2)" << std::endl;
  t2_2 = conj (t1_1) * conj (t1_2);

  if (verbose) std::cerr << "Type2 == Type2" << std::endl;
  if (norm(t2_1 - t2_2)/norm(t2_1) > tolerance)
  {
    std::cerr << "1=" << t2_1 << "\n2=" << t2_2 << std::endl;
    throw std::string ("test_matrix Error Conjugate");
  }

  //
  // test Hermitian
  //
  if (verbose) std::cerr << "\n*************** Testing Hermitian" << std::endl;

  t2_1 = t1_1;
  t2_2 = t1_2;

  if (verbose) std::cerr << "herm (t2_1 * t2_2)" << std::endl;
  temp1 = herm (t2_1 * t2_2);

  if (verbose) std::cerr << "herm (t2_2) * herm (t2_1)" << std::endl;
  temp2 = herm (t2_2) * herm (t2_1);

  if (verbose) std::cerr << "Type2 == Type2" << std::endl;

  if ( norm(temp2 - temp1)/norm(temp2) > tolerance) {
    std::cerr << "Error Hermitian"
      "\n   M1=" << t1_1 <<
      "\n   M2=" << t1_2 <<
      "\n  herm(M1*M2)=" << temp1 <<
      "\n  herm(M2)*herm(M1)=" << temp2 << std::endl;

    throw std::string ("test_matrix Error Hermitian");
  }

  //
  // test det(Type1) function
  //
  if (verbose) std::cerr << "\n*************** Testing determinant" << std::endl;

  t2_1 = t1_1;

  if (verbose) std::cerr << "det(Type2)" << std::endl;
  Element det1 (det(t2_1));

  if (verbose) std::cerr << "Type2 /= sqrt(det)" << std::endl;
  t2_1 /= sqrt(det1);

  if (verbose) std::cerr << "det(Type2)" << std::endl;
  Element det2 (det(t2_1));

  if (norm(det2 - one) > tolerance)
    throw std::string ("test_matrix "
		 "Error normalizing by square root of determinant"); 


  //
  // test norm and trace
  //
  if (verbose)
    std::cerr << "\n*************** Testing norm and trace" << std::endl;

  t2_2 = t1_1;

  if (verbose) std::cerr << "norm (Type2)" << std::endl;
  double variance1 = norm(t2_2);

  if (verbose) std::cerr << "herm (Type2)" << std::endl;
  t2_1 = t2_2*herm(t2_2);

  if (verbose) std::cerr << "trace (Type2)" << std::endl;
  Element variance2 (trace (t2_1));

  if ( norm(variance1 - variance2) > tolerance )
  {
    std::cerr << "Q=" << t2_2 << std::endl;
    std::cerr << "herm(Q)=" << herm(t2_2) << std::endl;
    std::cerr << "Q*herm(Q)=" << t2_1 << std::endl;
    std::cerr << "norm(Q) = " << variance1 << std::endl;
    std::cerr << "trace (Q * herm(Q)) = " << variance2 << std::endl;

    throw std::string ("test_matrix "
		 "Error norm(Q) != trace (Q * herm(Q))");
  }

  //
  // test Type1::inv method
  //
  if (verbose) std::cerr << "\n*************** Testing inverse" << std::endl;

  if (verbose) std::cerr << "Type1::inv" << std::endl;
  t2_2 = inv(t1_1);

  if (verbose) std::cerr << "Type2 * Type2" << std::endl;
  t2_1 = Type2(t1_1) * t2_2;

  double threshold = tolerance * max (norm(t1_1), norm(t2_2));

  if (verbose) std::cerr << "norm (Type2)" << std::endl;
  variance1 = norm (Type2::identity() - t2_1);

  if ( variance1 > threshold ) {
    std::cerr << "Error multipying Type1 by its inverse"
      "\n   Q=" << t1_1 <<
      "\n   inv(Q)=" << t2_2 <<
      "\n   Q*inv(Q)=" << t2_1 << std::endl;
    
    throw std::string ("test_matrix "
		  "Error multipying Type1 by its inverse");
  }
}

template <class Type1, class Type2, class Element>
void MatrixTest<Type1, Type2, Element>::runtest (unsigned nloops)
{
  //
  // test constructors
  //
  if (verbose)
    std::cerr << "\n*************** Testing constructors" << std::endl;

  if (verbose)
    std::cerr << "Type1 null construct" << std::endl;
  Type1 temp;

  if (verbose)
    std::cerr << "Element array (static)" << std::endl;

  Element random_elements[4];
  for (int i=0; i<4; i++)
    random_value (random_elements[i], scale);

  if (verbose)
    std::cerr << "Type1 construct from 4 Elements" << std::endl;
  Type1 t1_1 (random_elements[0], random_elements[1],
	      random_elements[2], random_elements[3]);

  if (verbose)
    std::cerr << "Type1 construct from Type1" << std::endl;
  Type1 t1_2 (t1_1);

  if (verbose)
    std::cerr << "TEST: Type1 == Type1" << std::endl;
  if (t1_1 != t1_2)
    throw std::string ("MatrixTest::runtest "
		  "Error in constructing Type1 from Type1");

  if (verbose)
    std::cerr << "Type1 = Type1" << std::endl;
  temp = t1_1;

  if (verbose)
    std::cerr << "TEST: Type1 == Type1" << std::endl;
  if (t1_1 != temp)
    throw std::string ("MatrixTest::runtest "
		 "Error setting Type1 equal to Type1");

  for (unsigned iloop=0; iloop<nloops; iloop++) {

    Type1 t1_1, t1_2;

    random_vector (t1_1, scale);
    random_vector (t1_2, scale);

    Type2 t2_null;
    Element null;

    try {
      test_matrix (t1_1, t1_2, t2_null, null, scale, tolerance, verbose);
    }
    catch (std::string& error) {
      std::cerr << "MatrixTest::runtest error on loop "<< iloop <<"/"<< nloops 
	   << std::endl;
      throw error;
    }

  }
}
