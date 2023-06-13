#include "Counter.h"
#include "Tools.h"

void stress_simple_2d_simplex(int_t n) {
    std::vector<std::vector<int_t>> A = {
        {0, -1},
        {-1, 0},
        {1, 1}};

    for (int i = 0; i <= n; ++i) {
        std::vector<int_t> b = {0, 0, i};
        mpf_class ans = (i*i + 3*i) / 2 + 1;
        auto res = count_integer_points(A, b);
        if (!(abs(res - ans) <= 1e-6))
        {
            std::cout << "Error! expected: " << ans << ", result: " << res << '\n';
            return;
        }
    }
    std::cout << "All tests passed!\n";
}

bool is_on_one_line(const std::vector<int_t>& x, const std::vector<int_t>& y) {
    return ((x[0] - x[1]) * y[2] + (y[1] - y[0]) * x[2] == x[0] * y[1] - x[1] * y[0]);
}

std::vector<int_t> get_inequlities(const std::vector<int_t> &x, const std::vector<int_t> &y,
                                    int_t i, int_t j, int_t k)
{
    int_t a = y[j] - y[i];
    int_t b = x[i] - x[j];
    int_t c = x[i] * y[j] - x[j] * y[i];
    if ((a * x[k] + b * y[k]) > c)
    {
        a = -a;
        b = -b;
        c = -c;
    }
    std::vector<int_t> res = {a, b, c};
    return res;
}

int_t brute_force_counter(const std::vector<std::vector<int_t>>& A, const std::vector<int_t>& b,
                        int_t min_x, int_t max_x, int_t min_y, int_t max_y) {
    int_t res = 0;
    for (int x = min_x; x <= max_x; ++x) {
        for (int y = min_y; y <= max_y; ++y) {
            bool flag = true;
            for (int i = 0; i < A.size(); ++i) {
                if (A[i][0] * x + A[i][1] * y > b[i]) {
                    flag = false;
                }
            }
            if (flag) ++res;
        }
    }
    return res;
}

void stress_random_simplex(int_t n) {
    int_t d = 3;
    int_t min_x = -10;
    int_t max_x = 20;
    int_t min_y = -10;
    int_t max_y = 20;
    for (int i = 0; i < n; ++i) {
        std::vector<int_t> x = gen_rand_vector(d, min_x, max_x);
        std::vector<int_t> y = gen_rand_vector(d, min_y, max_y);
        std::vector<std::vector<int_t>> A(d, std::vector<int_t>(d-1, 0));
        std::vector<int_t> b(d);

        if (is_on_one_line(x, y))
        {
            std::cout << "Points are on one line" << std::endl;
            return;
        }
        auto ineq = get_inequlities(x, y, 0, 1, 2);
        b[0] = ineq[d-1];
        for (int i = 0; i < d-1; ++i) {
            A[0][i] = ineq[i];
        }

        ineq = get_inequlities(x, y, 1, 2, 0);
        b[1] = ineq[d - 1];
        for (int i = 0; i < d-1; ++i)
        {
            A[1][i] = ineq[i];
        }

        ineq = get_inequlities(x, y, 2, 0, 1);
        b[2] = ineq[d - 1];
        for (int i = 0; i < d - 1; ++i)
        {
            A[2][i] = ineq[i];
        }

        print_vector<int_t>(b, "b");
        print_matrix(A, "A");
        auto ans = brute_force_counter(A, b, min_x, max_x, min_y, max_y);
        auto res = count_integer_points(A, b);
        if (abs(ans - res) > 1e-6)
        {
            std::cout << "Error! expected: " << ans << ", result: " << res << '\n';
            print_vector<int_t>(b, "b");
            print_matrix(A, "A");
        }
    }
    std::cout << "All tests passed!\n";
}

void stress_simple_3d_simplex(int_t n)
{
    std::vector<std::vector<int_t>> A = {
        {1, 1, 1},
        {-1, 0, 0},
        {0, -1, 0},
        {0, 0, -1}};
    mpf_class ans = 1;
    for (int i = 0; i <= n; ++i)
    {
        std::vector<int_t> b = {i, 0, 0, 0};
        auto res = count_integer_points(A, b);
        if (abs(ans - res) > 1e-6)
        {
            std::cout << "Error! expected: " << ans << ", result: " << res << '\n';
            print_vector<int_t>(b, "b");
            print_matrix(A, "A");
        }
        ans = ans + 3 * (i + 1) + (i * (i - 1)) / 2;
    }
    std::cout << "All tests passed!\n";
}

int main() {
    stress_simple_2d_simplex(10000);
    // stress_random_simplex(10);
    stress_simple_3d_simplex(10000);

    // std::vector<std::vector<int_t>> A = {
    //     {-5, 1},
    //     {1, -3},
    //     {3, 5}
    // };
    // std::vector<int_t> b = {-4, -2, 36};
    // auto res = count_integer_points(A, b);
    // std::cout << "Number of integer points: " << res << std::endl;

    // std::vector<std::vector<int_t>> A = {
    //     {1, 1, 1},
    //     {-1, 0, 0},
    //     {0, -1, 0},
    //     {0, 0, -1}};
    // std::vector<int_t> b = {4, 0, 0, 0};
    // auto res = count_integer_points(A, b);
    // std::cout << "Number of integer points: " << res << std::endl;

    // std::vector<std::vector<int_t>> A = {
    //     {3, 7},
    //     {-22, -7}};
    // std::vector<std::vector<int_t>> A = {
    //     {-23, 4},
    //     {18, 0}};
    // std::vector<std::vector<int_t>> P(2, std::vector<int_t>(2, 0));
    // std::vector<std::vector<int_t>> Q(2, std::vector<int_t>(2, 0));
    // SNF(A, P, Q);
    // print_matrix(A, "SNF");
    // print_matrix(P, "P");

    // std::vector<std::vector<int_t>> A = {
    //     {0, -1},
    //     {-1, 0},
    //     {1, 1}};
    // std::vector<int_t> b = {0, 0, 5};
    // std::cout << brute_force_counter(A, b, 0, 10, 0, 10) << std::endl;
}