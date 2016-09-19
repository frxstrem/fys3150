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
  const size_t Nlist[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
  const size_t Nlen = sizeof(Nlist) / sizeof(*Nlist);
  size_t N;

  // allocate dynamic memory in heap
  // (since 80MB is too much memory for stack)
  number *x = new number[10000001];
  number *a = new number[10000001];
  number *b = new number[10000001];
  number *c = new number[10000001];
  number *f = new number[10000001];

  // open file to write to
  FILE *fp = fopen("d.dat", "w");
  fprintf(fp, "log h\teps\n");

  // generate and solve problem for 10x10, 100x100 and 1000x1000
  for(size_t i = 0, N = Nlist[0]; i < Nlen; N = Nlist[++i]) {
    number h = 1.0 / (N + 1);

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

    // write maximum relative erro
    number maxeps = 0;
    for(size_t i = 1; i <= N; i++) {
      number u = G(x[i]);

      number eps = log10(fabs((u - v[i]) / u));

      if(eps > maxeps)
        maxeps = eps;
    }
    fprintf(fp, "%.0lf\t%.5E\n", log10(h), maxeps);
  }
  fclose(fp);

  // deallocate memory
  // (strictly speaking not necessary but let's be nice to the OS)
  delete[] x;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] f;
}
