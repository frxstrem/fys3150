#include <iostream>
#include <armadillo>
#include <ctime>
#include <cassert>
#include <unistd.h>

using namespace std;
using namespace arma;

const double rmin = 0;
const double rmax = 10;
size_t N = 41;
const double epsilon = 1e-10;

const double h = (rmax - rmin) / (N - 1);

// potential function
static constexpr double V(double r) {
  return r * r;
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

// sort and compare eigenvectors and eigenvalues
static double comp_eig(vec Avalues, mat Avectors, vec Bvalues, mat Bvectors) {
  assert(size(Avalues) == size(Bvalues));
  assert(size(Avectors) == size(Bvectors));
  assert(Avalues.n_elem == Avectors.n_rows && Avectors.is_square());

  size_t M = Avalues.n_elem;

  // sort A eigenvalues (quick and dirty algorithm)
  for(size_t i = 0; i < M; i++) {
    size_t k = i;
    for(size_t j = i + 1; j < M; j++) {
      if(Avalues(j) < Avalues(k))
        k = j;
    }
    if(i == k) continue;
    Avalues.swap_rows(i, k);
    Avectors.swap_cols(i, k);
  }

  // sort B eigenvalues (quick and dirty algorithm)
  for(size_t i = 0; i < M; i++) {
    size_t k = i;
    for(size_t j = i + 1; j < M; j++) {
      if(Bvalues(j) < Bvalues(k))
        k = j;
    }
    if(i == k) continue;
    Bvalues.swap_rows(i, k);
    Bvectors.swap_cols(i, k);
  }

  Bvalues.print();

  // fix antiparallell eigenvectors in A and B
  for(size_t i = 0; i < M; i++) {
    if(Avectors(0, i) * Bvectors(0, i) < 0)
      Bvectors.col(i) = - Bvectors.col(i);
  }

  // now that A and B have been sorted (so eigenvectors and eigenvalues should be in the same order)
  // calculate an error metric based on the dimensionality M
  const double err = norm(Avalues - Bvalues) / sqrt(M) + norm(Avectors - Bvectors) / M;

  // return error metric
  return err;
}

// solve with Jacobi's method
// A is input matrix, P is eigenvector matrix (with column vector eigenvectors), L is eigenvalue vector
static void jacobi_solve(const mat &A, mat &P, vec &L) {
  // assume A is N x N

  // B is similar matrix to A through matrix transforms
  mat B = A;

  // S is similarity transform matrix
  mat S(size(A));

  P.eye(size(A));

  clock_t time_start, time_end;
  double time_total = 0;
  size_t step = 0;
  while(true) {
    // timing
    time_start = clock();

    // calculate maximal off-diagonal element with respect to absolute value
    size_t k, l;
    const double a = maxoff(B, k, l);

    // if less than epsilon, stop
    if(a < epsilon) {
      break;
    }

    // calculate sin θ and cos θ
    const double tau = (B(l, l) - B(k, k)) / (2 * B(k, l));
    const double t = (tau >= 0 ? - tau - sqrt( 1 + tau * tau ) : - tau + sqrt( 1 + tau * tau ));
    const double c = 1 / sqrt( 1 + t * t );
    const double s = t * c;

    // generate rotation matrix
    S.eye();
    S(k, k) = S(l, l) = c;
    S(k, l) = s;
    S(l, k) = -s;

    // apply similarity transform
    B = S.t() * B * S;

    // apply S^T to P
    P = P * S;

    // timing
    time_end = clock();
    time_total += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    step++;
  }

  std::cout << "JACOBI DONE (N = " << N << ", steps = " << step << ", time per step = " << (1000 * time_total / step) << "ms)" << std::endl;

  // find eigenvalues
  L.set_size(N);
  for(size_t i = 0; i < N; i++)
    L(i) = B(i, i);
}

int main(int argc, char **argv) {
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
  mat A(N, N);
  for(size_t i = 0; i < N; i++) {
    A(i, i) = d(i);

    if(i > 0)
      A(i, i-1) = e;
    if(i < N - 1)
      A(i, i+1) = e;
  }

  // solve with Jacobi method
  vec eigenvalues;
  mat eigenvectors;
  jacobi_solve(A, eigenvectors, eigenvalues);

  // solve with Armadillo's eigsys
  vec arma_eigenvalues;
  mat arma_eigenvectors;
  eig_sym(arma_eigenvalues, arma_eigenvectors, A);

  double err = comp_eig(eigenvalues, eigenvectors, arma_eigenvalues, arma_eigenvectors);
  std::cout << err << std::endl;
}
