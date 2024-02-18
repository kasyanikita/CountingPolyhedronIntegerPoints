#include "snf_class.h"

#include <iostream>

SNFClass::SNFClass(const std::vector<std::vector<int_t>>& A)
    : A_(A), S_(A) {
  InitializeIdentityMatrices();
}

void SNFClass::CalculateSNF() {
  SNF(S_, P_, Q_);
  ExtractDiagonalOfS();
}

const std::vector<std::vector<int_t>>& SNFClass::GetA() const { return A_; }
const std::vector<std::vector<int_t>>& SNFClass::GetP() const { return P_; }
const std::vector<std::vector<int_t>>& SNFClass::GetS() const { return S_; }
const std::vector<std::vector<int_t>>& SNFClass::GetQ() const { return Q_; }
const std::vector<int_t>& SNFClass::GetDiagonalOfS() const {
  return s_diagonal_;
}

void SNFClass::InitializeIdentityMatrices() {
  if (!S_.empty()) {
    size_t size = S_.size();
    P_ = std::vector<std::vector<int_t>>(size, std::vector<int_t>(size, 0));
    Q_ = std::vector<std::vector<int_t>>(size, std::vector<int_t>(size, 0));

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