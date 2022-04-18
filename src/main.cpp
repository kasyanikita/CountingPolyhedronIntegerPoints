#include <iostream>
#include "ToddPolynomial.h"

int main() {
    int n = 3;
    Todd t = Todd(n, {1, 7, 4});
    for (int i = 0; i <= n; ++i) {
        std::cout << t.get_todd(i) << " ";
    }
    std::cout << std::endl;
}
