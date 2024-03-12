#include "counter.h"

using namespace GroupIP;

int main() {
  Matrix A = {{-5, 1}, {1, -3}, {3, 5}};
  Vector b = {-4, -2, 36};

  mpf_class res = 0;
  HyperplaneAvoidSolver hyperplane_avoid_vector(A);
  auto c = hyperplane_avoid_vector.get_vector(A.size());

  for (int i = 0; i < A.size(); ++i) {
    auto Asub = get_sub_matrix(A, i);
    auto bsub = get_sub_vector(b, i);
    res += PolyEvaluation(Asub, bsub, c);
  }

  std::cout << "Integer points = " << res << std::endl;
}