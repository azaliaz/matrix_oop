#ifndef __S21_MATRIX_OOP_H__
#define __S21_MATRIX_OOP_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#define EPS 1e-7

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix() noexcept;
  void SetNull();

  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;

  void CreateMatrix();
  void CopyMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const;
  void SubMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  double Determinant() const;
  void ClearMatrix();

  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);

  S21Matrix CreateMinor(int rows, int cols) const;

  double& operator()(int i, int j);
  double operator()(int i, int j) const;
  S21Matrix& operator=(S21Matrix&& other);

  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double num) const;
  bool operator==(const S21Matrix& other) const;
  S21Matrix operator+=(const S21Matrix& other);
  S21Matrix operator-=(const S21Matrix& other);
  S21Matrix operator*=(const S21Matrix& other);
  S21Matrix operator*=(const double num);
};

#endif
