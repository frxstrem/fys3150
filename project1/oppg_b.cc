/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#include <cstdio>
#include <cmath>

typedef double number;

/**
 *
 *
 */
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

/** Function [$f(x)$] */
number F(number x) {
  return 100 * exp(-10 * x);
}

/** Closed-form solution */
number G(number x) {
  return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

int main() {
  const size_t Nlist[] = { 10, 100, 1000 };
  size_t N;

  // generate and solve problem for 10x10, 100x100 and 1000x1000
  for(size_t i = 0, N = Nlist[0]; i < 3; N = Nlist[++i]) {
    number h = 1.0 / (N + 1);

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
