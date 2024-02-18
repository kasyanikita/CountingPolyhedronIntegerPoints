#ifndef SNF_CLASS_H_
#define SNF_CLASS_H_

#include <vector>

#include "SmithNormalForm.h"

class SNFClass {
 public:
  explicit SNFClass(const std::vector<std::vector<int_t>>&);

  void CalculateSNF();

  const std::vector<std::vector<int_t>>& GetA() const;
  const std::vector<std::vector<int_t>>& GetP() const;
  const std::vector<std::vector<int_t>>& GetS() const;
  const std::vector<std::vector<int_t>>& GetQ() const;
  const std::vector<int_t>& GetDiagonalOfS() const;

 private:
  std::vector<std::vector<int_t>> A_;
  std::vector<std::vector<int_t>> P_;  // Matrix P.
  std::vector<std::vector<int_t>> S_;  // Matrix S.
  std::vector<std::vector<int_t>> Q_;  // Matrix Q.

  std::vector<int_t> s_diagonal_;  // Diagonal elements of S.

  void InitializeIdentityMatrices();
  void ExtractDiagonalOfS();
};

#endif  // SNF_CLASS_H_
