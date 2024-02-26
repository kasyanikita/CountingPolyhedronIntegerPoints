#include "auxiliary_dynamic.h"

AuxiliaryDynamic::AuxiliaryDynamic(const std::vector<std::vector<int_t>>& A,
                                   const std::vector<int_t>& b,
                                   SNFClass& snf)
    : A_(A), b_(b), snf_(snf), h_(A_.size(), std::vector<int_t>(A_.size())) {}

const std::vector<GroupIP::GroupElement>& AuxiliaryDynamic::GetG() const {
  return g_;
}

const std::vector<std::vector<int_t>>& AuxiliaryDynamic::GetH() const {
  return h_;
}

void AuxiliaryDynamic::CalcG() {
  auto P = snf_.GetP();
  auto S = snf_.GetDiagonalOfS();

  int n = b_.size();
  std::vector<GroupElement> g_vec(n + 1, GroupElement(S));
  g_vec[n].assign(mat_vec_mult(P, b_));
  for (int i = 0; i < n; ++i) {
    g_vec[i].assign(get_mat_col(P, i));
  }

  g_ = g_vec;
}

void AuxiliaryDynamic::CalcH() {
  // init variables
  fmpz_mat_t A_flint;
  fmpz_t den;
  fmpz_t det;
  fmpz_mat_t Ainv;
  fmpz_mat_init(Ainv, A_.size(), A_[0].size());
  fmpz_mat_init(A_flint, A_.size(), A_[0].size());

  // std::vector to fmpz_mat_t
  for (int i = 0; i < A_.size(); ++i) {
    for (int j = 0; j < A_[i].size(); ++j) {
      auto val = fmpz_mat_entry(A_flint, i, j);
      *val = A_[i][j];
    }
  }

  // get inverse matrix
  fmpz_mat_det(det, A_flint);
  int res = fmpz_mat_inv(Ainv, den, A_flint);
  if (res == 0) {
    throw std::domain_error("Matrix A is singular");
  }

  // calculate adjugate matrix
  if (*det != *den) {
    std::cout << "det != den" << std::endl;
  }
  fmpz_divexact(den, det, den);
  fmpz_mat_scalar_divexact_fmpz(Ainv, Ainv, den);

  // fmpz_mat_t to std::vector
  for (int i = 0; i < A_.size(); ++i) {
    for (int j = 0; j < A_[i].size(); ++j) {
      auto val = fmpz_mat_entry(Ainv, i, j);
      h_[j][i] = fmpz_get_d(val);
    }
  }
}

void AuxiliaryDynamic::Init() {
    CalcG();
    CalcH();
}