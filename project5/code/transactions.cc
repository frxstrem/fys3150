#include <iostream>
#include <cstdio>
#include <random>
#include <algorithm>
#include <iterator>

#include "transactions.hh"

using namespace std;

// pseudorandom number generator
default_random_engine G;

// seed PRNG with true random number
unsigned int seed_simulation() {
  random_device source;
  unsigned int seed = source();
  G.seed(seed);
  return seed;
}

/**
 * Simulate transactions.
 *
 * Arguments:
 *  size_t N    - Number of agents
 *  size_t K    - Number of transactions
 *  double m0   - Initial amount of money for each agent
 *  double l    - Saving fraction λ (default: 0.0)
 *  double a    - Neighbor preference parameter α (default: 0.0)
 *  double g    - Previous interaction preference parameter γ (default: 0.0)
 *
 * Returns:
 *  Vector with the amount of money each agent has after all transactions
**/

// 5a)
vector<double> simulate_transactions(size_t N, size_t K, double m0) {
  // uniform distributions
  uniform_int_distribution<> idist(0, N-1);  // discrete distribution on [0,N)
  uniform_real_distribution<> rdist(0, 1);   // real distribution on [0,1)

  // list of agents' money
  vector<double> m(N, m0);

  // do switches
  for(size_t k = 0; k < K; k++) {
    // pick two random, different agents
    int i = idist(G), j = idist(G);
    while(i == j) j = idist(G);

    // calculate total money of the two agents
    double mtot = m[i] + m[j];

    // redistribute money between the two agents
    double e = rdist(G);
    m[i] = e * mtot;
    m[j] = (1 - e) * mtot;
  }

  return m;
}

// 5c)
vector<double> simulate_transactions(size_t N, size_t K, double m0, double l) {
  // uniform distributions
  uniform_int_distribution<> idist(0, N-1);  // discrete distribution on [0,N)
  uniform_real_distribution<> rdist(0, 1);   // real distribution on [0,1)

  // list of agents' money
  vector<double> m(N, m0);

  // do switches
  for(size_t k = 0; k < K; k++) {
    // pick two random, different agents
    int i = idist(G), j = idist(G);
    while(i == j) j = idist(G);

    // redistribute money between the two agents
    double e = rdist(G);
    double rm = (1 - l) * (e * m[j] - (1 - e) * m[i]);
    m[i] += rm;
    m[j] -= rm;
  }

  return m;
}

// // 5d)
// vec simulate_transactions(size_t N, size_t K, double m0, double l, double a) {
//   // uniform distributions
//   uniform_int_distribution<> idist(0, N-1);  // discrete distribution on [0,N)
//   uniform_real_distribution<> rdist(0, 1);   // real distribution on [0,1)
//
//   // list of agents' money
//   vec m = m0 * vec(N, fill::ones);
//
//   // do switches
//   for(size_t k = 0; k < K; k++) {
//     // pick two random agents
//     // with probability p_ij ∝ |m_i - m_j|^-α
//
//     // redistribute money between the two agents
//     double e = rdist(G);
//     double rm = (1 - l) * (e * m[j] - (1 - e) * m[i]);
//     m[i] += rm;
//     m[j] -= rm;
//   }
//
//   return m;
// }
//
// // 5e)
// vec simulate_transactions(size_t N, size_t K, double m0, double l, double a, double g) {
//   // uniform distributions
//   uniform_int_distribution<> idist(0, N-1);  // discrete distribution on [0,N)
//   uniform_real_distribution<> rdist(0, 1);   // real distribution on [0,1)
//
//   // list of agents' money
//   vec m = m0 * vec(N, fill::ones);
//
//   // do switches
//   for(size_t k = 0; k < K; k++) {
//     // pick two random, different agents
//     int i = idist(G), j = idist(G);
//     while(i == j) j = idist(G);
//
//     // redistribute money between the two agents
//     double e = rdist(G);
//     double rm = (1 - l) * (e * m[j] - (1 - e) * m[i]);
//     m[i] += rm;
//     m[j] -= rm;
//   }
//
//   return m;
// }
