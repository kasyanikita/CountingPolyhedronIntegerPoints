#include "ToddPoly_FFT.h"
#include "Dynamic.h"
#include "ExpPoly.h"

using namespace GroupIP;

std::vector<int_t> get_factorial(int_t n) {
    std::vector<int_t> fact(n + 1, 1);
    for (int i = 2; i <= n; ++i) {
        fact[i] = fact[i - 1] * i;
    }
    return fact;
}

void fill_powers(std::vector<std::vector<int_t>>& powers, std::vector<ExpPoly::exp_t> exps) {
    for (int i = 0; i < exps.size(); ++i) {
        for (int j = 1; j < powers[i].size(); ++j) {
            powers[i][j] = powers[i][j - 1] * exps[i];
        }
    }
}

std::vector<mpf_class> get_sum_exps(int_t n, std::vector<ExpPoly::exp_t> exps, // std::vector<mpf_class>
            std::vector<ExpPoly::coeff_t> coeffs, std::vector<mpf_class> todd)
{
    std::vector<std::vector<int_t>> powers(exps.size(), std::vector<int_t>(n + 1, 1));
    fill_powers(powers, exps);
    // print_matrix(powers, "Powers");
    auto fact = get_factorial(n);

    std::vector<mpf_class> summa;
    for (int i = 0; i < exps.size(); ++i) {
        mpf_class sum = 0;
        for (int j = 0; j <= n; ++j) {
            sum += (powers[i][j] * todd[n - j]) / fact[j];
        }
        // std::cout << "exp: " << exps[i] << ", sum: " << sum << std::endl;
        summa.push_back(coeffs[i] * sum);
    }
    return summa;
}

mpf_class count_simple_cone(std::vector<std::vector<int_t>> &A, std::vector<int_t> &b, std::vector<int_t> &c)
{
    // init dynamic
    Dynamic d(c, A, b);
    int_t n = A.size();
    d.init();
    d.new_start();

    // auto table = d.get_table();
    // for (int i = 0; i < table.size(); ++i) {
    //     for (int j = 0; j < table[i].size(); ++j) {
    //         std::cout << "(" << i << "," << j << "): " << table[i][j] << std::endl;
    //     }
    // }
    auto res = d.get_final_poly();
    auto den = d.get_den();
    // print_vector<int_t>(den, "denominator");

    // get coeffs and exps
    auto poly = res.get_poly();
    std::vector<ExpPoly::coeff_t> coeffs;
    std::vector<ExpPoly::exp_t> exps;
    for (const auto &[exp, coeff] : poly)
    {
        exps.push_back(exp);
        coeffs.push_back(coeff);
    }

    // get todd
    ToddPoly_FFT<int64_t, mpf_class> todd_poly(n, den);
    todd_poly.init();
    auto todd = todd_poly.get_todd();
    // print_vector<mpf_class>(todd, "Todd");

    // get sum alpha
    auto sum_alpha = get_sum_exps(n, exps, coeffs, todd);
    mpf_class numer = 0;
    int_t denom = 1;

    for (int i = 0; i < sum_alpha.size(); ++i) {
        numer += sum_alpha[i];
    }

    for (int i = 0; i < den.size(); ++i) {
        denom *= den[i]; 
    }

    // std::cout << numer / denom << std::endl;
    // print_vector<int_t>(den, "denominator");
    // std::cout << res << std::endl;
    return numer / denom;
}

mpf_class count_integer_points(std::vector<std::vector<int_t>>& A, std::vector<int_t>& b) {
    mpf_class res = 0;
    // auto det = get_determinant(A);
    // auto c = get_c_vector(A, A.size());
    std::vector<int_t> c = {1, -1};
    for (int i = 0; i < A.size(); ++i) {
        auto Asub = get_sub_matrix(A, i);
        auto bsub = get_sub_vector(b, i);
        res += count_simple_cone(Asub, bsub, c);
    }
    return res;
}

int main() {
    // std::vector<std::vector<int_t>> A1 = {
    //     {0, -1}, 
    //     {-1, 0}};
    // std::vector<int_t> b1 = {0, 0};
    // std::vector<int_t> c = {1, 4};
    // auto res = one_case(A1, b1, c);

    // std::vector<std::vector<int_t>> A2 = {
    //     {-1, 0},
    //     {1, 1}};
    // std::vector<int_t> b2 = {0, 5};
    // res = res + one_case(A2, b2, c);

    // std::vector<std::vector<int_t>> A3 = {
    //     {0, -1},
    //     {1, 1}};
    // std::vector<int_t> b3 = {0, 5};
    // res += one_case(A3, b3, c);

    // std::cout << res << std::endl;

    // std::vector<std::vector<int_t>> A2 = {
    //     {0, 1},
    //     {-1, 0}};
    // std::vector<int_t> b2 = {1, 0};
    // res = res + one_case(A2, b2, c);

    // std::vector<std::vector<int_t>> A3 = {
    //     {1, 0},
    //     {0, 1}};
    // std::vector<int_t> b3 = {1, 1};
    // res += one_case(A3, b3, c);

    // std::vector<std::vector<int_t>> A4 = {
    //     {0, -1},
    //     {1, 0}};
    // std::vector<int_t> b4 = {0, 1};
    // res += one_case(A4, b4, c);

    // std::cout << res << std::endl;

    // std::cout << "==============================================\n";
    // std::vector<std::vector<int_t>> A1 = {
    //     {-5, 1},
    //     {1, -3}};
    // std::vector<int_t> b1 = {-4, -2};
    // std::vector<int_t> c = {1, 1};
    // auto res = one_case(A1, b1, c);

    // std::cout << "==============================================\n";
    // std::vector<std::vector<int_t>> A2 = {
    //     {-5, 1},
    //     {3, 5}};
    // std::vector<int_t> b2 = {-4, 36};
    // res += one_case(A2, b2, c);

    // std::cout << "==============================================\n";
    // std::vector<std::vector<int_t>> A3 = {
    //     {1, -3},
    //     {3, 5}};
    // std::vector<int_t> b3 = {-2, 36};
    // res += one_case(A3, b3, c);

    std::vector<std::vector<int_t>> A = {
        {-5, 1},
        {1, -3},
        {3, 5}};
    std::vector<int_t> b = {-4, -2, 36};
    std::cout << count_integer_points(A, b) << std::endl;

    // std::vector<std::vector<int_t>> A = {
    //     {0, -1},
    //     {-1, 0},
    //     {1, 1}};
    // std::vector<int_t> b = {0, 0, 1};
    // std::cout << count_integer_points(A, b) << std::endl;

    // std::cout << (count_integer_points(A, b) - 50015001) << std::endl;
}