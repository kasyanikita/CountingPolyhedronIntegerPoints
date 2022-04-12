#ifndef INCLUDE_TODDPOLYNOMIAL_H_
#define INCLUDE_TODDPOLYNOMIAL_H_

#include <vector>
#include <cstdint>

class Todd {
 private:
    std::vector<double> bernulli;
    std::vector<uint64_t> calc_pascal(int);
    void calc_bernulli(int);
 public:
    explicit Todd(uint64_t = 2, std::vector<double> = {1});
    double get_bernulli(int);
};

#endif  // INCLUDE_TODDPOLYNOMIAL_H_
