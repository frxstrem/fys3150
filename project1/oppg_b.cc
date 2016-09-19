/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#include <cstdio>
#include <cmath>

#include "tridiagonal.hh"

typedef double number;

/** Function [$f(x)$] */
number F(number x) {
  return 100 * exp(-10 * x);
}

/** Closed-form solution */
number G(number x) {
  return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

int main() {
  // list of N values to test
  const size_t Nlist[] = { 10, 100, 1000 };
  const size_t Nlen = sizeof(Nlist) / sizeof(*Nlist);
  size_t N;

  // generate and solve problem for 10x10, 100x100 and 1000x1000
  for(size_t i = 0, N = Nlist[0]; i < Nlen; N = Nlist[++i]) {
    number h = 1.0 / (N + 1);

    // since N <= 1000, we can keep our arrays on the stack
    // (we can get away with this if we're compiling with g++)
    number x[N+1];
    number a[N+1];
    number b[N+1];
    number c[N+1];
    number f[N+1];

    // generate matrix [$A$] and vector [$\V{b}$]
    for(size_t i = 0; i <= N; i++) {
      x[i] = i * h;

      a[i] = -1;
      b[i] = 2;
      c[i] = -1;
      f[i] = h * h * F(x[i]);
    }

    // solve
    solve_tridiagonal(N + 1, a, b, c, f);

    // now f contains the solution, [$\V{v}$]
    // rename for convinience
    number *v = f;

    // write resulting function to file
    char filename[32];
    snprintf(filename, 32, "plots_b_%ld.dat", N);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "x\tv\n");
    for(size_t i = 0; i <= N; i++)
      fprintf(fp, "%.3lg\t%.3lg\n", x[i], v[i]);
    fclose(fp);
  }
}
