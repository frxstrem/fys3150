#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

int main(int argc, char **argv) {
  if(argc <= 2) {
    fprintf(stderr, "Usage: %s OUT-FILE NUM-AGENTS\n", argv[0]);
    return 0;
  }

  const char *filename = argv[1];
  int N = atoi(argv[2]);
  size_t K = 10000000;
  double m0 = 1.0;

  // pseudorandom number generator
  default_random_engine G;

  // seed PRNG with true random number
  random_device source;
  unsigned int seed = source();
  G.seed(seed);

  // uniform distributions
  uniform_int_distribution<> idist(0, N-1);  // discrete distribution on [0,N)
  uniform_real_distribution<> rdist(0, 1);   // real distribution on [0,1)

  // list of agents' money
  vector<double> m(N, m0);

  // sum of m^2
  double m2sq = N * m0 * m0;

  // lists of measured variance at different steps
  std::vector<int> step;
  std::vector<double> variance;

  // do switches
  size_t powerOfTen = 1;
  for(size_t k = 0; k < K; k++) {
    // pick two random, different agents
    int i = idist(G), j = idist(G);
    while(i == j) j = idist(G);

    // calculate total money of the two agents
    double mtot = m[i] + m[j];

    // redistribute money between the two agents
    double e = rdist(G);
    m2sq += (e * e + (1 - e) * (1 - e)) * mtot * mtot - m[i] * m[i] - m[j] * m[j];
    m[i] = e * mtot;
    m[j] = (1 - e) * mtot;

    if(k >= 100 * powerOfTen)
      powerOfTen *= 10;
    if(k % powerOfTen == 0) {
      // calculate variance
      double V = m2sq / N - m0 * m0;

      step.push_back(k);
      variance.push_back(V);
    }
  }

  // save to file
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "k\tV\n");
  for(size_t i = 0; i < step.size() && i < variance.size(); i++) {
    fprintf(fp, "%d\t%.3E\n", step[i], variance[i]);
  }
  fclose(fp);
}
