#include <iostream>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <mpi.h>

#include "misc.hh" // number, likely/unlikely macros
#include "IsingState.hh"
#include "metropolis.hh"

struct result {
  double t;
  double Eavg, absMavg, Cv, chi;
};

template <class Generator>
void runModelAtTemp(size_t M, size_t N, size_t K, number t, result &r, Generator g);

int main(int argc, char **argv) {
  int NProcesses, RankProcess;

  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
  MPI_Comm_rank (MPI_COMM_WORLD, &RankProcess);

  // read arguments
  const char *filename;
  FILE *fp = NULL;
  int NSpins, MonteCarloCycles;
  double InitialTemp, FinalTemp, TempStep;
  if(RankProcess == 0 && argc <= 5) {
    printf("Bad Usage: %s read output file, Number of spins, MC cycles, initial and final temperature and tempurate step\n", argv[0]);
    exit(1);
  }
  if(RankProcess == 0) {
    const char *filename = argv[1];
    NSpins = atoi(argv[2]);
    MonteCarloCycles = atoi(argv[3]);
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);

    // open file in main process
    // printf("%s\n", filename);
    fp = fopen(filename, "w");
    fprintf(fp, "t\tE\tabsM\tCv\tchi\n");
  }

  // broadcast to all nodes
  MPI_Bcast (&MonteCarloCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&NSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // random number generator
  std::default_random_engine g;

  // calculate number of steps
  ssize_t numberOfTempSteps = 1 + (FinalTemp - InitialTemp + 1e-10) / TempStep;

  // arrays of results
  result *combinedResults = new result[numberOfTempSteps];
  result *tmpResults = new result[numberOfTempSteps];

  for(size_t k = RankProcess; k < numberOfTempSteps; k += NProcesses) {
    number t = InitialTemp + TempStep * k;
    runModelAtTemp(NSpins, NSpins, MonteCarloCycles, t, tmpResults[k], g);
  }

  if(RankProcess == 0) {
    // combine results from all other processes
    for(size_t k = 0; k < numberOfTempSteps; k += NProcesses)
      combinedResults[k] = tmpResults[k];

    for(int p = 1; p < NProcesses; p++) {
      MPI_Recv(tmpResults, sizeof(result) * numberOfTempSteps, MPI_BYTE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for(size_t k = p; k < numberOfTempSteps; k += NProcesses)
        combinedResults[k] = tmpResults[k];
    }

    for(size_t k = 0; k < numberOfTempSteps; k++) {
      const result &r = combinedResults[k];
      fprintf(fp, "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\n", (double)r.t, (double)r.Eavg, (double)r.absMavg, (double)r.Cv, (double)r.chi);
    }
  } else {
    // send result to main process
    MPI_Send(tmpResults, sizeof(result) * numberOfTempSteps, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  }

  delete[] combinedResults;
  delete[] tmpResults;

  MPI_Finalize();
  if(fp != NULL)
    fclose(fp);
  return 0;
}

template <class Generator>
void runModelAtTemp(size_t M, size_t N, size_t K, number t, result &r, Generator g) {
  IsingState state = IsingState::randomize(M, N, g);

  // averages of energy, energy^2, magnetization, abs(magnetization) and magnetization^2
  number Eavg=0, E2avg=0, Mavg=0, absMavg=0, M2avg=0;

  // scale K by M*N
  K *= M*N;

  size_t n = 0;

  for(; n < K / 2; n++) {
    if(unlikely(n % 1000000 == 0))
      printf("t = %6.3f: %10lu / %-10lu (%5.1f%%) [discard]\n", (double)t, n, K, 100. * n / K);

    // run first half of steps but don't include in calculations of energy etc.
    metropolisStep(state, t, g);
  }

  for(; n < K; n++) {
    if(unlikely(n % 1000000 == 0))
      printf("t = %6.3f: %10lu / %-10lu (%5.1f%%)\n", (double)t, n, K, 100. * n / K);

    // run second half of steps
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

  r.t = t;
  r.Eavg = Eavg;
  r.absMavg = absMavg;

  // calculate specific heat capacity (in units of k) and susceptibility (in units of 1/J)
  r.Cv = (E2avg - (Eavg * Eavg)) / (t * t);
  r.chi = (M2avg - (Mavg * Mavg)) / t;
}
