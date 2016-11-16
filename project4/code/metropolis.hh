#ifndef METROPOLIS_HH
#define METROPOLIS_HH

#include <random>
#include "misc.hh" // number, integer
#include "IsingState.hh"

template <class Generator>
bool metropolisStep(IsingState &state, number t, Generator &g) {
  size_t M = state.M;
  size_t N = state.N;

  // uniform distribution on [0,1]
  std::uniform_real_distribution<> dist(0, 1);

  // uniform integer distributions on [0,M) and [0,N)
  std::uniform_int_distribution<> iDist(0, M - 1), jDist(0, N - 1);

  // get random indices to flip
  size_t i = iDist(g);
  size_t j = jDist(g);

  // flip spin and calculate energy difference (in units J)
  integer dE = state.flip(i, j);

  bool acceptState = false;

  if(dE <= 0) {
    // new state has lower energy, so accept this new state
    acceptState = true;
  } else {
    // new state has higher energy, so only accept with some probability
    number w = exp(- dE / t);
    number r = dist(g);

    if(r <= w)
      acceptState = true;
  }

  if(!acceptState)
    state.flip(i, j);

  return acceptState;
}

#endif // METROPOLIS_HH
