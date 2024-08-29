#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}
int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  CreateMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  CreateMatrix();
  CopyMatrix(other);
}
S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.SetNull();
}
void S21Matrix::SetNull() {
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}
void S21Matrix::CreateMatrix() {
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument("Rows and columns must be positive");
  }
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(other.matrix_[i][j] - matrix_[i][j]) > EPS) return false;
    }
  }
  return true;
}
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.GetRows() || cols_ != other.GetCols()) {
    throw std::length_error("Matrix dimensions must be the same for addition.");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (other.matrix_ == nullptr) {
    throw std::length_error("Matrix doesn't exist");
  }
  if ((rows_ != other.rows_) || (cols_ != other.cols_)) {
    throw std::length_error("Different dimensions");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  if (matrix_ == nullptr) {
    throw std::length_error("Matrix doesn't exist");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.GetRows()) {
    throw std::length_error("Matrix dimensions must agree for multiplication.");
  }
  S21Matrix result(rows_, other.GetCols());
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  (*this) = result;
}
void S21Matrix::SetRows(int rows) {
  if (rows < 1) {
    throw std::length_error(
        "Incorrect input, number of rows should be more than zero");
  }
  if (rows != rows_) {
    S21Matrix tempM(rows, cols_);
    int min_rows = std::min(rows_, rows);
    ;
    for (int i = 0; i < min_rows; i++) {
      for (int j = 0; j < cols_; j++) {
        tempM.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = std::move(tempM);
    tempM.ClearMatrix();
  }
}
void S21Matrix::ClearMatrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = nullptr;
}

void S21Matrix::SetCols(int cols) {
  if (cols < 1) {
    throw std::length_error(
        "Incorrect input, number of rows should be more than zero");
  }
  if (cols != cols_) {
    S21Matrix tempM(rows_, cols);
    int min_cols = std::min(cols_, cols);
    ;
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < min_cols; j++) {
        tempM.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = std::move(tempM);
    tempM.ClearMatrix();
  }
}

S21Matrix S21Matrix::Transpose() const {
  if (matrix_ == nullptr) {
    throw std::length_error("Matrix doesn't exist");
  }
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CreateMinor(int row, int column) const {
  if (row < 0 || row >= rows_ || column < 0 || column >= cols_) {
    throw std::out_of_range("Row or column index out of range.");
  }

  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int i = 0, minor_i = 0; i < rows_; ++i) {
    if (i == row) continue;

    for (int j = 0, minor_j = 0; j < cols_; ++j) {
      if (j == column) continue;

      minor.matrix_[minor_i][minor_j++] = matrix_[i][j];
    }
    ++minor_i;
  }

  return minor;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is ​​not square");
  }

  S21Matrix result(rows_, cols_);

  if (rows_ == 1) {
    result.matrix_[0][0] = 1;
  } else {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        S21Matrix minor = CreateMinor(i, j);

        double det = minor.Determinant();
        result.matrix_[i][j] = pow(-1, (i + j)) * det;
      }
    }
  }

  return result;
}
double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is ​​not square");
  }
  if (rows_ == 1) {
    return matrix_[0][0];
  }
  if (rows_ == 2) {
    return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  }
  double det = 0;
  for (int i = 0; i < rows_; ++i) {
    S21Matrix minor = CreateMinor(0, i);

    double minor_det = minor.Determinant();
    det += (i % 2 == 0 ? 1 : -1) * matrix_[0][i] * minor_det;
  }

  return det;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is not square");
  }

  double det = Determinant();
  if (fabs(det) < EPS) {
    throw std::length_error("Matrix determinant is 0");
  }

  S21Matrix complements = CalcComplements();
  S21Matrix adjugate = complements.Transpose();

  S21Matrix inverse(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      inverse.matrix_[i][j] = adjugate.matrix_[i][j] / det;
    }
  }

  return inverse;
}

double S21Matrix::operator()(int i, int j) const {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::length_error("Incorrect index!");
  }
  return matrix_[i][j];
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::length_error("Incorrect index!");
  }
  return matrix_[i][j];
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    ClearMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    CreateMatrix();
    CopyMatrix(other);
  }
  return *this;
}

// S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
//   S21Matrix tempM(other);
//   *this = std::move(tempM);
//   return *this;
// }
S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    ClearMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}
S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}
S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}
bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}
S21Matrix::~S21Matrix() noexcept { ClearMatrix(); }