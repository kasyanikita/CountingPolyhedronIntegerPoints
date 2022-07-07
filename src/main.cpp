#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include "ToddFFT.h"

template<typename T>
T random(T range_from, T range_to) {
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<T> distr(range_from, range_to);
    return distr(generator);
}

int main () {
    std::vector<mpz_class> v = {1, 2, 3, 4, 5};
    int n = 2;
    int m = 10;
    Todd<mpz_class, mpq_class> t(m, v);
    t.init();
    auto todd = t.get_todd();
    for (int i = 0; i < 11; ++i) {
        std::cout << mpf_class(todd[i]) << " ";
    }
    std::cout << std::endl;

    std::vector<mpz_class> vfft = {1, 2, 3, 4, 5};
    ToddFFT<mpz_class, mpq_class> tfft(m, vfft);
    tfft.init();
    auto toddfft = tfft.get_todd_mpf();
    for (int i = 0; i < 11; ++i) {
        std::cout << toddfft[i] << " ";
    }
    std::cout << std::endl;
}