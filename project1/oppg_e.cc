/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#include <cstdio>
#include <cmath>
#include <ctime>
#include <functional>

#include "tridiagonal.hh"
#include "lib.h"

typedef double number;

/** Function [$f(x)$] */
number F(number x) {
  return 100 * exp(-10 * x);
}

/** Closed-form solution */
number G(number x) {
  return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

double timing(std::function<void()> func, int n = 1) {
  clock_t start, finish;
  clock_t total_time = 0;

  for(int i = 0; i < n; i++) {
    start = clock();
    func();
    finish = clock();
    total_time += finish - start;
  }

  return (1.0 * total_time / CLOCKS_PER_SEC / n);
}

int main() {
  // list of N values to test
  const size_t Nlist[] = { 10, 100, 1000 };
  const size_t Nlen = sizeof(Nlist) / sizeof(*Nlist);
  size_t N;

  // number of loops for each matrix size, for timing
  int M = 10;

  // open file to write to
  FILE *fp = fopen("e.dat", "w");
  fprintf(fp, "n\tctime\tlutime\n");

  // allocate A matrix dynamically
  number **A = new number*[1002];
  for(size_t i = 0; i < 1002; i++)
    A[i] = new number[1002];

  // generate and solve problem for 10x10, 100x100 and 1000x1000
  for(size_t i = 0, N = Nlist[0]; i < Nlen; N = Nlist[++i]) {
    printf("N=%ld\n", N);

    clock_t start, finish;
    clock_t total_ctime = 0;
    clock_t total_lutime = 0;
    for(int m = 0; m < M; m++) {
      printf("  M=%d\n", m);
      number h = 1.0 / (N + 1);

      // since N <= 1000, we can keep our arrays on the stack
      // (we can get away with this if we're compiling with g++)
      number x[N+1];
      number a[N+1];
      number b[N+1];
      number c[N+1];
      number f[N+1];
      number f2[N+1];
      int idx[N+1];

      // generate matrix [$A$] and vector [$\V{b}$]
      for(size_t i = 0; i <= N; i++) {
        idx[i] = i;
        x[i] = i * h;

        a[i] = A[i+1][i] = -1;
        b[i] = A[i][i]   =  2;
        c[i] = A[i][i+1] = -1;
        f[i] = f2[i]     = h * h * F(x[i]);
      }

      // solve with custom method
      start = clock();
      solve_tridiagonal(N + 1, a, b, c, f);
      finish = clock();
      total_ctime += finish - start;

      // solve with LU decomposition
      start = clock();
      double lutime = timing([&] () -> void {
        double d;
        ludcmp(A, N + 1, idx, &d);
        lubksb(A, N + 1, idx, f2);
      });
      finish = clock();
      total_lutime += finish - start;
    }

    double ctime = 1.0 * total_ctime / CLOCKS_PER_SEC / M;
    double lutime = 1.0 * total_lutime / CLOCKS_PER_SEC / M;

    // print times
    fprintf(fp, "%ld\t%.5E\t%.5E\n", N, 1e3 * ctime, 1e3 * lutime);
  }
  fclose(fp);

  // deallocate A
  for(size_t i = 0; i < 1002; i++)
    delete[] A[i];
  delete[] A;
}
