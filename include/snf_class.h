#ifndef SNF_CLASS_H_
#define SNF_CLASS_H_

#include <eigen3/Eigen/Dense>
#include <vector>

#include "SNF/DF_module.h"
#include "SNF/gmp_supp.h"
#include "global_defs.h"

// #include "SmithNormalForm.h"
using namespace GroupIP;

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
  using snf_int_type = mpz_class;
  using eigen_mat = Eigen::Matrix<snf_int_type, Eigen::Dynamic, Eigen::Dynamic>;
  using my_DF_algo = LatLib::my_diagonal_form_algo<snf_int_type>;

  std::vector<std::vector<int_t>> A_;
  std::vector<std::vector<int_t>> P_;  // Matrix P.
  std::vector<std::vector<int_t>> S_;  // Matrix S.
  std::vector<std::vector<int_t>> Q_;  // Matrix Q.

  std::vector<std::vector<int_t>> eigen2vector(eigen_mat);

  std::vector<int_t> s_diagonal_;  // Diagonal elements of S.

  void InitializeIdentityMatrices();
  void ExtractDiagonalOfS();
};

#endif  // SNF_CLASS_H_
