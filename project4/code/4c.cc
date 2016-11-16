#include <iostream>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "misc.hh" // number, integer
#include "IsingState.hh"
#include "metropolis.hh"

template <class Generator>
void runCalc(FILE *fp, size_t K, number t, IsingState &state, Generator &g);

int main(int argc, char **argv ) {
  // number of Monte Carlo cycles
  size_t K = 1e9;

  // random number generator
  std::default_random_engine g;

  // t = 1.0, ordered initial state
  {
    FILE *fp = fopen("4c-1.0-ordered.dat", "w");
    IsingState state = IsingState(20, 20);
    runCalc(fp, K, 1.0, state, g);
    fclose(fp);
  }

  // t = 1.0, randomized initial state
  {
    FILE *fp = fopen("4c-1.0-random.dat", "w");
    IsingState state = IsingState::randomize(20, 20, g);
    runCalc(fp, K, 1.0, state, g);
    fclose(fp);
  }

  // t = 2.4, ordered initial state
  {
    FILE *fp = fopen("4c-2.4-ordered.dat", "w");
    IsingState state = IsingState(20, 20);
    runCalc(fp, K, 2.4, state, g);
    fclose(fp);
  }

  // t = 2.4, randomized initial state
  {
    FILE *fp = fopen("4c-2.4-random.dat", "w");
    IsingState state = IsingState::randomize(20, 20, g);
    runCalc(fp, K, 2.4, state, g);
    fclose(fp);
  }
}

template <class Generator>
void runCalc(FILE *fp, size_t K, number t, IsingState &state, Generator &g) {
  size_t M = state.M;
  size_t N = state.N;

  fprintf(fp, "n\tE\tabsM\taccepted\n");

  // cumulative energy and absolute magnetization of the Monte Carlo cycles
  number Ecuml = 0, absMcuml = 0;
  size_t totalAccepted = 0;
  size_t powerOfTen = 1;
  for(size_t n = 1; n <= K; n++) {
    bool accepted = metropolisStep(state, t, g);

    if(accepted)
      totalAccepted++;

    number E = state.energy();
    number m = state.magnetization();

    Ecuml    += E;
    absMcuml += (m >= 0 ? m : -m);

    // print mean energy and mean absolute magnetization to file, 90 steps for each power of ten
    if(n >= 100 * powerOfTen)
      powerOfTen *= 10;
    if(n % powerOfTen == 0) {
      fprintf(fp, "%lu\t%.5E\t%.5E\t%lu\n", n, (double)Ecuml / n, (double)absMcuml / n, totalAccepted);
    }
  }
}
