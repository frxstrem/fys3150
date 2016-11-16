#include <iostream>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "misc.hh" // number, likely/unlikely macros
#include "IsingState.hh"
#include "metropolis.hh"

int main(int argc, char **argv) {
  // number of Monte Carlo cycles
  size_t K = 100000;

  // temperature interval and number of data points
  number t0 = 0, t1 = 10;
  size_t tN = 100;

  // random number generator
  std::default_random_engine g;

  FILE *fp = fopen("4b.dat", "w");
  fprintf(fp, "t\tE\tabsM\tCv\tchi\n");

  for(int k = 1; k <= tN; k++) {
    number t = t0 + (t1 - t0) * k / tN;

    IsingState state = IsingState::randomize(2, 2, g);

    // averages of energy, energy^2, magnetization, abs(magnetization) and magnetization^2
    number Eavg=0, E2avg=0, Mavg=0, absMavg=0, M2avg=0;

    for(size_t n = 0; n < K; n++) {
      metropolisStep(state, t, g);

      number E = state.energy();
      number m = state.magnetization();

      // add to cumulative values of energy and magnetization
      Eavg += E / K;
      E2avg += E * E / K;
      Mavg += m / K;
      absMavg += (m >= 0 ? m : -m) / K;
      M2avg += m * m / K;
    }

    // calculate specific heat capacity (in units of k) and susceptibility (in units of 1/J)
    number Cv = (E2avg - (Eavg * Eavg)) / (t * t);
    number chi = (M2avg - (Mavg * Mavg)) / t;

    // print to file
    fprintf(fp, "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\n", (double)t, (double)Eavg, (double)absMavg, (double)Cv, (double)chi);
  }

  fclose(fp);
}
