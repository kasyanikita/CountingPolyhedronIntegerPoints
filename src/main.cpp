#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "ToddPolynomial.h"

int main ()
{
    uint64_t n = 20;
    std::vector<mpf_class> v = {2102, 123.5, 123.435};
    Todd<mpz_class, mpf_class> t(n, v);
    for (int i = 0; i < n; ++i) {
        std::cout << t.get_todd(i) << std::endl;
    }
    return 0;
}
