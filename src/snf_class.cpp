#include "snf_class.h"

#include <iostream>

SNFClass::SNFClass(const Matrix& A) : A_(A) {
  InitializeIdentityMatrices();
}

void SNFClass::CalculateSNF() {
  int size = A_.size();
  eigen_mat A(size, size);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      A(i, j) = A_[i][j];
    }
  }

  eigen_mat P, S, Q;
  my_DF_algo DF_algo;
  std::tie(P, S, Q) = DF_algo.compute(A);

  P_ = eigen2vector(P);
  Q_ = eigen2vector(Q);
  S_ = eigen2vector(S);

  ExtractDiagonalOfS();
}

std::vector<std::vector<int_t>> SNFClass::eigen2vector(eigen_mat matrix) {
  int size = matrix.rows();
  std::vector<std::vector<int_t>> res(size, std::vector<int_t>(size));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      res[i][j] = matrix(i, j).get_si();
    }
  }

  return res;
}

const Matrix& SNFClass::GetA() const { return A_; }
const Matrix& SNFClass::GetP() const { return P_; }
const Matrix& SNFClass::GetS() const { return S_; }
const Matrix& SNFClass::GetQ() const { return Q_; }
const Vector& SNFClass::GetDiagonalOfS() const {
  return s_diagonal_;
}

void SNFClass::InitializeIdentityMatrices() {
  if (!S_.empty()) {
    size_t size = S_.size();
    P_ = Matrix(size, Vector(size, 0));
    Q_ = Matrix(size, Vector(size, 0));

    for (size_t i = 0; i < size; ++i) {
      P_[i][i] = 1;
      Q_[i][i] = 1;
    }
  }
}

void SNFClass::ExtractDiagonalOfS() {
  s_diagonal_.clear();
  for (size_t i = 0; i < S_.size(); ++i) {
    s_diagonal_.push_back(S_[i][i]);
  }
}