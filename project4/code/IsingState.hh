#ifndef ISING_STATE_HH
#define ISING_STATE_HH

#include <random>
#include "misc.hh" // number, integer

class IsingState {
public:
  // dimensions M,N
  const size_t M, N;

private:
  // total number of spins S = M * N
  const size_t S;

  // list of spins (size: S)
  integer *s;

  // calculated energy and magnetization
  integer mEnergy, mMagnetization;

public:
  IsingState(size_t M, size_t N)
    : M(M), N(N), S(M*N), mEnergy(0), mMagnetization(0) {
    s = new integer[S];
    for(size_t k = 0; k < S; k++)
      s[k] = 1;

    recalculateEnergyAndMagnetization();
  }

  IsingState(const IsingState &state)
    : M(state.M), N(state.N), S(state.M * state.N), mEnergy(state.mEnergy), mMagnetization(state.mMagnetization) {
      s = new integer[S];
      for(size_t k = 0; k < S; k++)
        s[k] = state.s[k];
  }

  ~IsingState() {
    delete[] s;
  }

  template <class Generator>
  static IsingState randomize(size_t M, size_t N, Generator &g) {
    std::uniform_int_distribution<> dist(0, 1);

    IsingState state(M, N);
    for(size_t k = 0; k < state.S; k++)
        state.s[k] = 2 * dist(g) - 1; // uniform probability of { -1, 1 }

    state.recalculateEnergyAndMagnetization();
    return state;
  }

  inline integer get(int i, int j) const {
    size_t k = ((j+N) % N) * M + ((i+M) % M);
    return s[k];
  }

  /**
   * Make flip, return change in energy
  **/
  integer flip(int i, int j) {
    size_t k = ((j+N) % N) * N + ((i+M) % M);

    // get old spin at (i,j)
    integer oldSpin = get(i,j);

    // flip spin
    s[k] *= -1;

    // calculate change in energy
    integer dEnergy = 2 * oldSpin * ( get(i+1,j) + get(i-1,j) + get(i,j+1) + get(i,j-1) );

    // calculate change in magnetization
    integer dMagnetization = - 2 * oldSpin;

    mEnergy += dEnergy;
    mMagnetization += dMagnetization;

    return dEnergy;
  }

  integer energy() const { return mEnergy; }
  integer magnetization() const { return mMagnetization; }

  void recalculateEnergyAndMagnetization() {
    // calculate energy
    mEnergy = 0;
    for(size_t i = 0; i < M; i++) {
      for(size_t j = 0; j < N; j++) {
        mEnergy -= get(i,j) * get(i+1,j);
        mEnergy -= get(i,j) * get(i,j+1);
      }
    }

    // calculate magnetiztaion
    mMagnetization = 0;
    for(size_t i = 0; i < M; i++)
      for(size_t j = 0; j < N; j++)
        mMagnetization += get(i,j);
  }
};

#endif // ISING_STATE_HH
