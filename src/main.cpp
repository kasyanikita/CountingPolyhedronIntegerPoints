#include <iostream>
#include "ToddPolynomial.h"

int main() {
    int n = 10;
    Todd t = Todd(n);
    for (int i = 0; i < n + 1; ++i) {
        std::cout << t.get_bernulli(i) << std::endl;
    }
    return 0;
}
