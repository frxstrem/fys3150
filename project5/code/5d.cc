#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <armadillo>
#include <random>

#include <mpi.h>

#include "transactions.hh"
#include "histograms.hh"

using namespace std;
using namespace arma;

int main(int argc, char **argv) {
  int mpiSize, mpiRank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if(argc <= 6) {
    if(mpiRank == 0) {
      fprintf(stderr, "Usage: %s OUTPUT-FILE NUM-AGENTS NUM-TRANSACTIONS NUM-RUNS INITIAL-MONEY [LAMBDA] [ALPHA] [S-FACTOR]\n", argv[0]);
    }
    MPI_Finalize();
    exit(0);
  }

  // read parameters
  const char *filename;
  int N;      // number of agents
  int K;      // number of transactions
  int R;      // number of runs
  double m0;  // initial money for each agent
  double l;   // lambda parameter
  double a;   // alpha parameter
  double S;   // "S-factor", affects accuracy and slowness of simulations with α≠0 or γ≠0 (default: 1.0)
  if(mpiRank == 0) {
    // mandatory arguments
    filename = argv[1];
    N = atoi(argv[2]);
    K = atoi(argv[3]);
    R = atoi(argv[4]);
    m0 = atof(argv[5]);

    // optional arguments
    l  = (argc > 6 ? atof(argv[6]) : 0.0);
    a  = (argc > 7 ? atof(argv[7]) : 0.0);
    S  = (argc > 9 ? atof(argv[9]) : 1.0);
  }

  if(mpiRank == 0) {
    fprintf(stderr, "Running with parameters:\n");
    fprintf(stderr, "  N = %d\n", N);
    fprintf(stderr, "  K = %d\n", K);
    fprintf(stderr, "  R = %d\n", R);
    fprintf(stderr, " m0 = %.3f\n", m0);
    fprintf(stderr, "  λ = %.3f\n", l);
    fprintf(stderr, "  α = %.3f\n", a);
    fprintf(stderr, "  S = %.1f\n", S);
  }

  // broadcast parameters to all nodes
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&R, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&l, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&S, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // reseed simulation with true random seed
  // (so that the different processes don't make the exact same simulations)
  seed_simulation();

  // list of all agents after each run
  vector<double> m(N * R, 1);

  // find start and end indices for range of runs to be done
  // by this process
  size_t runsPerProcess = R / mpiSize;
  size_t startRun = mpiRank * runsPerProcess;
  size_t endRun = startRun + runsPerProcess;
  if(endRun > R) endRun = R;
  if(mpiRank == mpiSize - 1 && endRun < R) endRun = R;

  // do all runs (parallelized)
  for(size_t r = startRun; r < endRun; r++) {
    // simulate transactions
    // (store resulting money distribution in a subvector)
    vector<double> m_dist = simulate_transactions(N, K, m0, l, S, a);

    copy(begin(m_dist), end(m_dist), begin(m) + (r * N));
  }

  // send all data to process #0
  double *mptr = m.data();
  if(mpiRank == 0) {
    for(int p = 1; p < mpiSize; p++) {
      size_t pStartRun = p * runsPerProcess;
      size_t pEndRun   = pStartRun + runsPerProcess;
      if(pEndRun > R) pEndRun = R;
      if(p == mpiSize - 1 && pEndRun < R) pEndRun = R;

      MPI_Recv(&mptr[N * pStartRun], N * (pEndRun - pStartRun), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  } else {
    MPI_Send(&mptr[N * startRun], N * (endRun - startRun), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  // save histogram data to file (for plotting)
  if(mpiRank == 0)
    save_histogram(m, 0.01, filename);

  MPI_Finalize();
  return 0;
}
