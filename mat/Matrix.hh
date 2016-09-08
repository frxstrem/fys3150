/**
 * Copyright (c) 2016 Fredrik Østrem
 * This code is license under MIT license (see LICENSE for details)
 */

#ifndef MATRIX_HH_
#define MATRIX_HH_

#ifndef __cplusplus
#error Not C++
#endif

#include <stdexcept>
#include <functional>
#include <type_traits>

namespace mat {
  typedef double number;

  /**
   * Utility class for commonly used types:
   */
  template <typename T>
  struct NumberType {
    static constexpr T zero     = 0;
    static constexpr T identity = 1;
  };

  /**
   * Matrix class
   */
  template <size_t M, size_t N, typename T = number>
  class Matrix {
  private:
    T data[M][N];

  public:
    /** Create new matrix */
    Matrix() {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          data[i][j] = NumberType<T>::zero;
    }

    /** Create matrix from C array */
    Matrix(const T (&array)[M][N]) {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          data[i][j] = array[i][j];
    }

    /** Copy matrix */
    Matrix(const Matrix &a) {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          data[i][j] = a.data[i][j];
    }

    /** Get diagonal identity matrix */
    static Matrix Identity() {
      Matrix a;
      for(size_t i = 0; i < M && i < N; i++)
        a.data[i][i] = NumberType<T>::identity;
      return a;
    }

    /** Get reference to matrix element */
    T &operator ()(size_t i, size_t j) {
      if(i < 0 || i >= M || j < 0 || j >= N)
        throw std::range_error("Matrix index out of range");

      return data[i][j];
    }

    /** Get reference to matrix element */
    const T &operator ()(size_t i, size_t j) const {
      if(i < 0 || i >= M || j < 0 || j >= N)
        throw std::range_error("Matrix index out of range");

      return data[i][j];
    }

    /** Apply a function to each element of a matrix */
    template <typename U>
    static Matrix apply(const Matrix<M,N,U> &a, auto func) {
      Matrix b;
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          b.data[i][j] = func(a.data[i][j]);
      return b;
    }
    static Matrix &applySelf(Matrix<M,N,T> &a, auto func) {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          a.data[i][j] = func(a.data[i][j]);
      return a;
    }

    /** Apply a function pairwise to each element of two matrices */
    template <typename U, typename V>
    static Matrix apply(const Matrix<M,N,U> &a, const Matrix<M,N,V> &b, auto func) {
      Matrix c;
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          c.data[i][j] = func(a.data[i][j], b.data[i][j]);
      return c;
    }
    template <typename U>
    static Matrix &applySelf(Matrix<M,N,T> &a, const Matrix<M,N,U> &b, auto func) {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          a.data[i][j] = func(a.data[i][j], b.data[i][j]);
      return a;
    }

    /** Apply a function with a left argument to each element of a matrix */
    template <typename U, typename V>
    static Matrix apply(const U &a, const Matrix<M,N,V> &b, auto func) {
      Matrix c;
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          c.data[i][j] = func(a, b.data[i][j]);
      return c;
    }

    /** Apply a function with a right argument to each element of a matrix */
    template <typename U, typename V>
    static Matrix apply(const Matrix<M,N,U> &a, const V &b, auto func) {
      Matrix c;
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          c.data[i][j] = func(a.data[i][j], b);
      return c;
    }
    template <typename U>
    static Matrix &applySelf(Matrix<M,N,T> &a, U &b, auto func) {
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          a.data[i][j] = func(a.data[i][j], b);
      return a;
    }
  };

  /**
   * Matrix operators
   */

  /** Return copy of matrix */
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> operator +(const Matrix<M,N,T> &a) {
    return Matrix<M,N,T>(a);
  }

  /** Add two matrices */
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> operator +(const Matrix<M,N,T> &a, const Matrix<M,N,T> &b) {
    return Matrix<M,N,T>::apply(a, b, std::plus<>());
  }
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> &operator +=(Matrix<M,N,T> &a, const Matrix<M,N,T> &b) {
    return Matrix<M,N,T>::applySelf(a, b, std::plus<>());
  }

  /** Negate a matrix */
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> operator -(const Matrix<M,N,T> &a) {
    return Matrix<M,N,T>::apply(a, std::negate<>());
  }

  /** Subtract two matrices */
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> operator -(const Matrix<M,N,T> &a, const Matrix<M,N,T> &b) {
    return Matrix<M,N,T>::apply(a, b, std::minus<>());
  }
  template <size_t M, size_t N, typename T>
  Matrix<M,N,T> &operator -=(Matrix<M,N,T> &a, const Matrix<M,N,T> &b) {
    return Matrix<M,N,T>::applySelf(a, b, std::minus<>());
  }

  /** Multiply matrix by a scalar */
  template <size_t M, size_t N, typename T, typename U>
  Matrix<M,N,T> operator *(const U &a, const Matrix<M,N,T> &b) {
    return Matrix<M,N,T>::apply(a, b, std::multiplies<>());
  }
  template <size_t M, size_t N, typename T, typename U>
  Matrix<M,N,T> operator *(const Matrix<M,N,T> &a, const U &b) {
    return Matrix<M,N,T>::apply(a, b, std::multiplies<>());
  }
  template <size_t M, size_t N, typename T, typename U>
  Matrix<M,N,T> &operator *=(Matrix<M,N,T> &a, const U &b) {
    return Matrix<M,N,T>::applySelf(a, b, std::multiplies<>());
  }

  /** Divide matrix by a scalar */
  template <size_t M, size_t N, typename T, typename U>
  Matrix<M,N,T> operator /(const Matrix<M,N,T> &a, const U &b) {
    return Matrix<M,N,T>::apply(a, b, std::divides<>());
  }
  template <size_t M, size_t N, typename T, typename U>
  Matrix<M,N,T> &operator /=(Matrix<M,N,T> &a, const U &b) {
    return Matrix<M,N,T>::applySelf(a, b, std::divides<>());
  }

  /** Matrix multiplication */
  template <size_t M, size_t K, size_t N, typename T>
  Matrix<M,N,T> operator *(const Matrix<M,K,T> &a, const Matrix<K,N,T> &b) {
    Matrix<M,N,T> c;

    for(size_t i = 0; i < M; i++)
      for(size_t k = 0; k < K; k++)
        for(size_t j = 0; j < N; j++)
          c(i,j) += a(i,k) * b(k,j);

    return c;
  }
};

#endif // MATRIX_HH_
