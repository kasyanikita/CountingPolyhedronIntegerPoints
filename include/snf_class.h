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
  explicit SNFClass(const Matrix&);

  void CalculateSNF();

  const Matrix& GetA() const;
  const Matrix& GetP() const;
  const Matrix& GetS() const;
  const Matrix& GetQ() const;
  const Vector& GetDiagonalOfS() const;

 private:
  using snf_int_type = mpz_class;
  using eigen_mat = Eigen::Matrix<snf_int_type, Eigen::Dynamic, Eigen::Dynamic>;
  using my_DF_algo = LatLib::my_diagonal_form_algo<snf_int_type>;

  Matrix A_;
  Matrix P_;  // Matrix P.
  Matrix S_;  // Matrix S.
  Matrix Q_;  // Matrix Q.

  Matrix eigen2vector(eigen_mat);

  Vector s_diagonal_;  // Diagonal elements of S.

  void InitializeIdentityMatrices();
  void ExtractDiagonalOfS();
};

#endif  // SNF_CLASS_H_
