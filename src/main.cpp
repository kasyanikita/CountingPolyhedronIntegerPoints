#include "Dynamic.h"
#include "ExpPoly.h"
#include "SmithNormalForm.h"
#include "Tools.h"
#include "flint/fmpz_mat.h"

using namespace GroupIP;
// using namespace MSVP;

// void test_my_DF(size_t retries = 100, size_t n = 100, unsigned int interval = 100)
// {
//     default_random_engine eng;
//     int error_number = 0;
//     for (size_t retr_i = 0; retr_i < retries; retr_i++)
//     {
//         int_mat A = generate_random_mat(n, interval, eng);
//         for (int i = 0; i < A.rows(); ++i)
//         {
//             for (int j = 0; j < A.cols(); ++j)
//             {
//                 std::cout << A(i, j);
//                 std::cout << " ";
//             }
//             std::cout << '\n';
//         }
//         std::cout << "\n\n";
//         int_mat P, S, Q;
//         tie(P, S, Q) = my_DF(A);

//         if (!is_diagonal(S))
//         {
//             cout << "S is not diagonal: " << endl;
//             // cout << S << endl;
//         }

//         if (P * A * Q != S)
//         {
//             error_number++;
//             cout << "DF computation contains errors. P A Q is: " << endl;
//             // cout << P * A * Q << endl;

//             // cout << "S is: " << endl;
//             // cout << S << endl;
//         }
//         for (int i = 0; i < S.rows(); ++i) {
//             for (int j = 0; j < S.cols(); ++j) {
//                 std::cout << S(i, j);
//                 std::cout << " ";
//             }
//             std::cout << '\n';
//         }
//     }
// }

// void check_flint_snf() {
//     std::vector<std::vector<int_t>> A = {
//         {3, 5, 10, 6, 19},
//         {-3, -5, 10, 37, -17},
//         {5, 5, 3, -10, 42},
//         {85, -63, 7, -5, 8},
//         {8, 7, -6, 11, -23}};
//     fmpz_mat_t Af;
//     fmpz_mat_t S;
//     fmpz_mat_init(S, A.size(), A[0].size());
//     fmpz_mat_init(Af, A.size(), A[0].size());
//     for (int i = 0; i < A.size(); ++i) {
//         for (int j = 0; j < A[i].size(); ++j) {
//             auto val = fmpz_mat_entry(Af, i, j);
//             *val = A[i][j];
//         }
//     }
//     fmpz_mat_snf(S, Af);
//     fmpz_mat_print_pretty(S);
// }

void one_case(std::vector<std::vector<int_t>>& A, std::vector<int_t>& b, std::vector<int_t>& c) {
    Dynamic d(c, A, b);
    d.init();
    d.new_start();
}

    int main()
{

    // std::vector<std::vector<int_t>> A = {
    //     {-5, 1},
    //     {1, -3},
    //     {3, 5}
    // };
    // auto c = get_c_vector(A, 9);
    // print_vector<int_t>(c, "c");

    std::cout << "==============================================\n";
    std::vector<std::vector<int_t>> A1 = {
        {-5, 1},
        {1, -3}
    };
    std::vector<int_t> b1 = {-4, -2};
    std::vector<int_t> c = {1, 1};
    one_case(A1, b1, c);

    std::cout << "==============================================\n";
    std::vector<std::vector<int_t>> A2 = {
        {-5, 1},
        {3, 5}};
    std::vector<int_t> b2 = {-4, 36};
    one_case(A2, b2, c);

    std::cout << "==============================================\n";
    std::vector<std::vector<int_t>> A3 = {
        {1, -3},
        {3, 5}};
    std::vector<int_t> b3 = {-2, 36};
    one_case(A3, b3, c);

    // std::vector<std::vector<int_t>> A = {
    //     {-5, 1},
    //     {1, -3},
    //     {3, 5}
    // };
    // auto c = get_c_vector(A, 9);
    // print_vector<int_t>(c, "c");
}