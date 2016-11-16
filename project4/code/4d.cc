#include <iostream>
#include <random>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "misc.hh" // number, integer
#include "IsingState.hh"
#include "metropolis.hh"

template <class Generator>
void runCalc(FILE *fp, size_t Ksteady, size_t K, number t, IsingState &state, Generator &g);

int main(int argc, char **argv ) {
  // number of Monte Carlo cycles
  size_t Ksteady = 1e7;
  size_t K = 1e8;

  // random number generator
  std::default_random_engine g;

  // t = 1.0
  {
    FILE *fp = fopen("4d-1.0.dat", "w");
    IsingState state = IsingState::randomize(20, 20, g);
    runCalc(fp, Ksteady, K, 1.0, state, g);
    fclose(fp);
  }

  // t = 2.4
  {
    FILE *fp = fopen("4d-2.4.dat", "w");
    IsingState state = IsingState::randomize(20, 20, g);
    runCalc(fp, Ksteady, K, 2.4, state, g);
    fclose(fp);
  }
}

template <class Generator>
void runCalc(FILE *fp, size_t Ksteady, size_t K, number t, IsingState &state, Generator &g) {
  size_t M = state.M;
  size_t N = state.N;

  // skip until steady state
  for(size_t n = 0; n < Ksteady; n++)
    metropolisStep(state, t, g);

  // store energies in bins
  std::map<integer, integer> Ebins;
  number Eavg = 0, E2avg = 0;
  for(size_t n = 0; n < K; n++) {
    metropolisStep(state, t, g);

    integer E = state.energy();

    Ebins[E] = Ebins[E] + 1;

    Eavg += (number)E / K;
    E2avg += (number)(E * E) / K;
  }

  // calculate variance σE^2 and standard deviation σE
  number varE = E2avg - (Eavg * Eavg);
  number stdE = sqrt(varE);

  printf("t=%8.4g, var(E)=%8.4g, σ_E=%8.4g\n", (double)t, (double)varE, (double)stdE);

  // save bins in file
  fprintf(fp, "E\tN\tP\n");
  for(std::map<integer, integer>::const_iterator it = Ebins.cbegin(); it != Ebins.cend(); it++) {
    integer E = it->first;  // energy
    integer n = it->second; // number of states counted with energy

    fprintf(fp, "%.5E\t%.5E\t%.5E\n", (double)E, (double)n, (double)n / K);
  }
}
