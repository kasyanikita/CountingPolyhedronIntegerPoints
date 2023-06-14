#include <chrono>
#include <fstream>
#include <vector>
#include <iostream>
#include <flint/fmpz_mat.h>
#include "Counter.h"

int_t max_rand = 10;

std::vector<std::vector<int_t>> get_simple_A(int_t n){
    std::vector<std::vector<int_t>> A(n + 1, std::vector<int_t>(n, 0));
    for (int i = 0; i < A[0].size(); ++i) {
        A[0][i] = get_random_number(1, max_rand);
    }
    A[0][1] = max_rand;
    for (int i = 1; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            if (i == j + 1) {
                A[i][j] = -1;
            } else {
                A[i][j] = 0;
            }
        }
    }
    return A;
}

std::vector<int_t> get_simple_b(int_t n)
{
    std::vector<int_t> b(n+1, 0);
    b[0] = get_random_number(1, max_rand);
    return b;
}

void load_to_file(const std::vector<mpf_class>& times, const std::string& path) {
    std::ofstream out(path);
    for (auto x : times)
    {
        out << x << std::endl;
    }
    out.close();
}

void fixed_det(int_t n) {
    std::vector<mpf_class> times;
    for (int i = 2; i <= n; ++i) {
        std::cout << i << std::endl;
        mpf_class sum = 0;
        for (int j = 0; j < 10; ++j) {
            auto A = get_simple_A(i);
            auto b = get_simple_b(i);
            // print_matrix(A, "A");
            // print_vector<int_t>(b, "b");
            auto c = get_c_vector(A, A.size());
            auto begin = std::chrono::high_resolution_clock::now();
            auto res = count_integer_points(A, b, c);
            // std::cout << res << std::endl;
            auto end = std::chrono::high_resolution_clock::now();
            sum += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
        }
        times.push_back(sum / 10);
    }
    print_vector<mpf_class>(times, "times");
    load_to_file(times, "data/time_fixed_det.txt");
}

// void fixed_dim(int_t det)
// {
//     std::vector<mpf_class> times;
//     for (int i = 1; i <= det; ++i)
//     {
//         mpf_class sum = 0;
//         std::cout << i << std::endl;
//         for (int j = 0; j < 1; ++j)
//         {
//             auto A = get_simple_A(20);
//             auto b = get_simple_b(20);
//             b[0] = i;
//             for (int k = 0; k < A[0].size(); ++k) {
//                 A[0][k] = i;
//             }
//             // print_matrix(A, "A");
//             // auto c = get_c_vector(A, A.size() * A.size());
//             auto begin = std::chrono::high_resolution_clock::now();
//             auto c = get_c_vector(A, A.size() * A.size());
//             auto res = count_integer_points(A, b, c);
//             // std::cout << res << std::endl;
//             auto end = std::chrono::high_resolution_clock::now();
//             sum += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
//         }
//         times.push_back(sum / 1);
//     }
//     print_vector<mpf_class>(times, "times");
//     load_to_file(times, "data/time_fixed_dim.txt");
// }

int main() {

    fixed_det(30);
    // fixed_dim(300);

    // std::vector<std::vector<int_t>> A = {
    //     {-5, 1},
    //     {1, -3}};
    // std::vector<int_t> b = {-4, -2};
    // std::vector<int_t> c = {1, 1};
    // Dynamic d(c, A, b);
    // d.init();
    // d.start();
}