/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#ifndef TRIDIAGONAL_HH_
#define TRIDIAGONAL_HH_

#ifndef __cplusplus
#error Not compiling with C++
#endif

/**
 * Solve a N×N matrix equation with a tridiagonal matrix.
 * a, b, c contain the lower, main and upper diagonals, respectively.
 * f contains the right-side vector of the equation.
 * a, b, c, f are modified by the function.
 */
template <typename number>
void solve_tridiagonal(size_t N, number *a, number *b, number *c, number *f) {
  // eliminate lower diagonal through forward substitution
  for(size_t i = 0; i < N - 1; i++) { // flops: 5 * (N - 1)
    if(b[i] == 0) continue;
    const number s = a[i] / b[i]; // 1 flop

    a[i]   = 0; // = (b[i] / a[i]) * a[i] - b[i]
    b[i+1] = b[i+1] - s * c[i];   // 2 flops
    f[i+1] = f[i+1] - s * f[i];   // 2 flops
  }
  // now a[i] == 0 for all i

  // eliminate upper diagonal through backward substitution
  for(ssize_t i = N - 2; i >= 0; i--) { // 3 * (N - 1)
    if(b[i+1] == 0) continue;
    const number s = c[i] / b[i + 1]; // 1 flop

    // b[i] = b[i] - (c[i] / b[i + 1]) * a[i] = b[i]
    c[i] = 0; // = c[i] - (c[i] / b[i + 1]) * b[i + 1]
    f[i] = f[i] - s * f[i + 1];       // 2 flops
  }
  // now c[i] == 0 for all i

  // divide diagonal elements
  for(size_t i = 0; i < N; i++) { // 2 * N
    const number s = b[i];

    b[i] /= s;
    f[i] /= s;
  }

  // total flops: 10 * N - 8
}

#endif
