#include <iostream>
#include <cstdio>
#include <armadillo>
#include <ctime>
#include <limits>
#include <unistd.h>

#include "comp_eig.hh"

using namespace std;
using namespace arma;

const double rmin = 1e-10;
const double rmax = 60;
const double epsilon = 1e-10;

// potential function (non-intercting)
static constexpr double V0(double r) {
  return r * r;
}

// potential function (interacting)
static constexpr double Vw(double r, double w) {
  return w * w * r * r + 1 / r;
}

// calculate maximal off-diagonal element with respect to absolute value
static double maxoff(const mat &A, size_t &i, size_t &j) {
  i = 1; j = 0; // indices of maximal lement
  double a = 0; // maximal element absolute value

  SizeMat s = size(A);

  for(size_t k = 0; k < s.n_rows; k++) {
    for(size_t l = 0; l < s.n_cols; l++) {
      if(k == l) continue;

      const double x = A(k, l);
      const double y = x * x;
      if(y > a)
        i = k, j = l, a = y;
    }
  }

  return a;
}

static void apply_rot_col(size_t k, size_t l, double t, mat &M) {
  if(k == l) return;

  const double c = 1 / sqrt( 1 + t * t );
  const double s = t * c;

  colvec kvec = c * M.col(k) - s * M.col(l);
  colvec lvec = s * M.col(k) + c * M.col(l);

  M.col(k) = kvec;
  M.col(l) = lvec;
}

static void apply_rot_row(size_t k, size_t l, double t, mat &M) {
  if(k == l) return;

  const double c = 1 / sqrt( 1 + t * t );
  const double s = t * c;

  rowvec kvec = c * M.row(k) - s * M.row(l);
  rowvec lvec = s * M.row(k) + c * M.row(l);

  M.row(k) = kvec;
  M.row(l) = lvec;
}

// solve with Jacobi's method
// A is input matrix, P is eigenvector matrix (with column vector eigenvectors), L is eigenvalue vector
static void jacobi_solve(const mat &A, mat &P, vec &L, size_t &steps, double &step_time) {
  // get N from A matrix
  const size_t N = A.n_rows;
  assert(A.n_cols == N);

  std::cout << "N=" << N << std::endl;

  // B is similar matrix to A through matrix transforms
  mat B = A;

  // S is similarity transform matrix
  mat S(size(A));

  P.eye(size(A));

  clock_t time_start, time_end;
  double time_total = 0;
  steps = 0;
  while(true) {
    // timing
    time_start = clock();

    // calculate maximal off-diagonal element with respect to absolute value
    size_t k, l;
    const double a = maxoff(B, k, l);

    // if less than epsilon, stop
    if(a < epsilon)
      break;

    // calculate sin θ and cos θ
    const double tau = (B(l, l) - B(k, k)) / (2 * B(k, l));
    const double t = (tau >= 0 ? - tau - sqrt( 1 + tau * tau ) : - tau + sqrt( 1 + tau * tau ));

    // apply similarity transform
    //   B = S.t() * B * S;
    apply_rot_col(k, l, t, B);
    apply_rot_row(k, l, t, B);

    // apply S^T to P
    //   P = P * S;
    apply_rot_col(k, l, t, P);

    // timing
    time_end = clock();
    time_total += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    steps++;
  }

  // export steps and time per step
  step_time = time_total / steps;

  std::cout << "JACOBI DONE (N = " << N << ", steps = " << steps << ", time per step = " << (1000 * step_time) << "ms)" << std::endl;

  // find eigenvalues
  L.set_size(N);
  for(size_t i = 0; i < N; i++)
    L(i) = B(i, i);
}

int run_program(size_t N, double w, const std::function<double(double)> &V) {
  const double h = (rmax - rmin) / N;

  // array of ρ valuses
  vec r(N);
  for(size_t i = 0; i < N; i++)
    r(i) = rmin + h * i;

  // calculate e_i, which are all the same
  double e = - 1 / (h * h);

  // calculate d_i, which depend on V(ρ_i)
  vec d(N);
  for(size_t i = 0; i < N; i++)
    d(i) = 2 / (h * h) + V(r[i]);

  // calculate matrix A
  mat A(N, N, arma::fill::zeros);
  for(size_t i = 0; i < N; i++) {
    A(i, i) = d(i);

    if(i > 0)
      A(i, i-1) = e;
    if(i < N - 1)
      A(i, i+1) = e;
  }

  // solve with Armadillo's eig_sym
  vec eigenvalues;
  mat eigenvectors;
  eig_sym(eigenvalues, eigenvectors, A);

  // normalize w.r.t. N
  eigenvectors /= sqrt(h);

  // find lowest eigenvalue
  size_t lowest_index = 0;
  double lowest_eigenvalue = eigenvalues(0);
  for(size_t i = 1; i < N; i++) {
    double l = eigenvalues(i);
    if(l < lowest_eigenvalue) {
      lowest_index = i; lowest_eigenvalue = l;
    }
  }

  // save lowest eigenvector to file
  char filename[20];
  snprintf(filename, 20, "d-%.2lf.dat", w);
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "r          u\n");
  for(size_t i = 0; i < N; i++) {
    fprintf(fp, "%-10.5g %-10.5g\n", r(i), eigenvectors(i, lowest_index));
  }
  fclose(fp);
}

int main(int argc, char **argv) {
  const size_t N = 800;
  const double Wvalues[] = { 0.01, 0.5, 1, 5 };
  const size_t Wlen = sizeof(Wvalues) / sizeof(*Wvalues);

  // run for non-interacting case
  std::cout << "running for non-interactive case" << std::endl;
  run_program(N, 0, V0);

  // run for interacting cases
  for(size_t i = 0; i < Wlen; i++) {
    double W = Wvalues[i];
    std::cout << "running with W = " << W << std::endl;
    const auto V = [W] (double r) -> double { return Vw(r, W); };
    run_program(N, W, V);
  }
}
