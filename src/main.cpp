#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "ToddPolynomial.h"

int main () {
    int n = 2;
    std::vector<mpz_class> v = {1, 2, 10};
    Todd<mpz_class, mpf_class> t(n, v);
    t.init();
    auto todd = t.get_todd();
    for (int i = 0; i <= n; ++i) {
        std::cout << todd[i] << std::endl;
    }
    return 0;
}
