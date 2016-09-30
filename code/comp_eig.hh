#ifndef COMP_EIG_HH
#define COMP_EIG_HH

#include <armadillo>
#include <cassert>

// sort and compare eigenvectors and eigenvalues
static double comp_eig(arma::vec Avalues, arma::mat Avectors, arma::vec Bvalues, arma::mat Bvectors) {
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

  // fix antiparallell eigenvectors in A and B
  for(size_t i = 0; i < M; i++) {
    if(arma::dot(Avectors.col(i), Bvectors.col(i)) < 0)
      Bvectors.col(i) = - Bvectors.col(i);
  }

  // fix maximal difference between eigenvalues or eigenvectors
  const double values_error = arma::abs(Avalues - Bvalues).max();
  const double vectors_error = arma::abs(Avectors - Bvectors).max();

  // return error metric
  return std::max(values_error, vectors_error);
}

#endif
