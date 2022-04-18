#ifndef INCLUDE_TODDPOLYNOMIAL_H_
#define INCLUDE_TODDPOLYNOMIAL_H_

#include <vector>
#include <cstdint>

class Todd {
 private:
    std::vector<double> xi;
    std::vector<double> bernulli;
    std::vector<double> todd;
    std::vector<uint64_t> calc_pascal(int) const;
    void calc_bernulli(int);
    void calc_todd(uint64_t, const std::vector<double>&);
 public:
    explicit Todd(uint64_t = 2, const std::vector<double>& = {1});
    double get_bernulli(int) const;
    double get_todd(uint64_t) const;
};

#endif  // INCLUDE_TODDPOLYNOMIAL_H_
