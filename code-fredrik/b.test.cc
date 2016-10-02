#include <iostream>
#include <functional>
#include <armadillo>
#include <string>

// include code from b.cc (except main function)
#define NO_MAIN
#include "b.cc"
#undef NO_MAIN

using namespace std;
using namespace arma;

const std::string CHECK_MARK = "✓";
const std::string BALLOT_X = "✗";

struct test { std::function<bool()> func; std::string name; };

const double TOLERANCE = 1e-13;

/**
 * Test maxoff() on a symmetric matrix
**/
bool test_maxoff_sym() {
  // create known symmetric matrix
  mat M(5, 5, arma::fill::zeros);

  M(0,0) = 3;
  M(1,1) = 9;
  M(2,2) = 17;
  M(3,3) = 0;
  M(4,4) = -9;
  M(2,0) = M(0,2) = 9;
  M(2,1) = M(1,2) = -5;
  M(3,0) = M(0,3) = -11;
  M(3,1) = M(1,3) = 10;
  M(3,2) = M(2,3) = -7;
  M(4,0) = M(0,4) = 5;

  // expected values
  size_t exp_k = 3, exp_l = 0; // maximal element position
  double exp_a = 121;          // square of maximal element

  // find maximal element
  size_t k, l;
  double a = maxoff(M, k, l);

  // compare (up to order of k and l, since M is symmetric)
  return (abs(a - exp_a) < TOLERANCE && ((k == exp_k && l == exp_l) || (k == exp_l && l == exp_k)));
}

/**
 * Check that jacobi_step() conserves orthonormality of S.
**/
bool test_ortho() {
  // if a matrix preserves scalar product (or orthonormality) then it must be an orthogonal matrix,
  // so S * S^T = I

  // set of random, pre-calculated tau, k and l values (where k != l)
  const double tau_values[] = { 0.5 };
  const double kl_values[][2] = { { 1, 2 } };
  const size_t count = 1;

  // initial P matrix (identity) and constant identity matrix I
  mat P(5, 5, arma::fill::eye);
  const mat I(5, 5, arma::fill::eye);

  // for each step, P' = P * S where S is calculated from a random tau value (see above list)
  for(size_t i = 0; i < count; i++) {
    apply_rot_col(kl_values[i][0], kl_values[i][1], tau_values[i], P);

    // check that P * P^T = I, within tolerance
    if(arma::abs(P * P.t() - I).max() >= TOLERANCE)
      return false;
  }

  // for each step, P' = S^T * P where S is calculated from a random tau value (see above list)
  for(size_t i = 0; i < count; i++) {
    apply_rot_row(kl_values[i][0], kl_values[i][1], tau_values[i], P);

    // check that P * P^T = I, within tolerance
    if(arma::abs(P * P.t() - I).max() >= TOLERANCE)
      return false;
  }

  return true;
}

int main() {
  // define tests
  const struct test tests[] = {
    { test_maxoff_sym, "Maximal diagonal test on symmetric matrix" },
    { test_ortho, "Orthonormality test" },
  };
  const size_t test_count = 2;

  // run tests
  std::cout << "Running tests:" << std::endl;
  for(size_t i = 0; i < test_count; i++) {
    auto test = tests[i];

    std::cout << "  " << test.name << std::endl;

    if(test.func()) {
      std::cout << "    " << CHECK_MARK << " test passed" << std::endl;
    } else {
      std::cout << "    " << BALLOT_X << " test failed" << std::endl;
    }
  }
}
