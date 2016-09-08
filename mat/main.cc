/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#include <iostream>

#include "Matrix.hh"

using namespace mat;

template <size_t M, size_t N, typename T>
void print_matrix(const Matrix<M,N,T> &mat) {
  for(size_t i = 0; i < M; i++) {
    if(i == 0)
      std::cout << '[';
    else
      std::cout << ' ';

    std::cout << '[';
    for(size_t j = 0; j < N; j++) {
      if(j > 0)
        std::cout << ", ";
      std::cout << mat(i, j);
    }
    std::cout << ']';

    if(i == M - 1)
      std::cout << ']';
    else
      std::cout << ',';
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

int main(int argc, char **argv) {
}
