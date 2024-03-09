#ifndef TOOLS_H_
#define TOOLS_H_

#include <flint/fmpz_mat.h>
#include <gmpxx.h>

#include <random>
#include <stdexcept>
#include <vector>

#include "group_element.h"
#include "ExpPoly.h"
#include "global_defs.h"

using namespace GroupIP;

namespace GroupIP {
    
template <class T>
void print_vector(std::vector<T> v, std::string name);
void print_matrix(std::vector<std::vector<int_t>> &M, std::string name);
std::vector<int_t> scalar_vector_mul(int_t s, const std::vector<int_t> &v);
int_t dot_product(const std::vector<int_t> &a, const std::vector<int_t> &b);
std::vector<int_t> get_mat_col(const std::vector<std::vector<int_t>> &M,
                               int k);
std::vector<int_t> mat_vec_mult(const std::vector<std::vector<int_t>> &M,
                                const std::vector<int_t> &v);
int_t calc_s(const GroupElement &g, const GroupElement &g0, uint_t r0);
GroupElement get_group_element_by_index(int_t idx, std::vector<int_t> &S);
ExpPoly vertex_normalize(const ExpPoly &num, std::vector<int_t> &den,
                         const Matrix &A, const Vector &b, const Vector &c);
std::vector<int_t> get_denominator(const Vector &c,
                                   std::vector<GroupElement> &g,
                                   const std::vector<std::vector<int_t>> &h);
}  // namespace GroupIP

int_t get_random_number(int_t min, int_t max);
int_t get_determinant(const GroupIP::Matrix &A);
std::vector<std::vector<int_t>> get_sub_matrix(
    std::vector<std::vector<int_t>> &A, int_t row_to_exclude);
std::vector<int_t> get_sub_vector(std::vector<int_t> &b, int_t row_to_exclude);
std::vector<int_t> gen_rand_vector(int_t n, int_t a, int_t b);
std::vector<std::vector<int_t>> get_adjugate(
    const std::vector<std::vector<int_t>> &A);
std::vector<int_t> adj_c_multiply(std::vector<std::vector<int_t>> &Aadj,
                                  std::vector<int_t> &c);
bool is_any_zero(std::vector<int_t> &v);
bool check_sub(std::vector<std::vector<int_t>> &A_sub, std::vector<int_t> &c);
bool check_subs(std::vector<std::vector<std::vector<int_t>>> &A_subs,
                std::vector<int_t> &c);
std::vector<int_t> get_c_vector(std::vector<std::vector<int_t>> &A,
                                int_t alpha);
std::vector<mpf_class> get_sum_exps(
    int_t n, std::vector<ExpPoly::exp_t> exps,  // std::vector<mpf_class>
    std::vector<ExpPoly::coeff_t> coeffs, std::vector<mpf_class> todd);
std::vector<int_t> get_factorial(int_t n);
void fill_powers(std::vector<std::vector<int_t>> &powers,
                 std::vector<ExpPoly::exp_t> exps);
                 
#endif  // TOOLS_H_