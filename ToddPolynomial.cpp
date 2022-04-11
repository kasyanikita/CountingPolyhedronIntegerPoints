#include <iostream>
#include <vector>

class Todd {
private:
    std::vector<double> bernulli;

    std::vector<uint64_t> calc_pascal(int n) {
        std::vector<uint64_t> res(n + 1, 1);
        int i;
        for (i = 0; i < n; ++i) {
            res[i + 1] = res[i] * (n - i) / (i + 1);
        }
        return res;
    }

    void calc_bernulli(int n) {
        bernulli.push_back(1);
        bernulli.push_back(-0.5);
        int i;
        int k;
        double sum;
        for (i = 2; i <= n; ++i) {
            if (i % 2 == 1) {
                bernulli.push_back(0);
                continue;
            }
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (k = 0; k < i; ++k) {
                sum += pascal[k + 2] * bernulli[i - k - 1];
            }
            bernulli.push_back((-1.0 / (i + 1)) * sum);
        }
    }

public:
    Todd(uint64_t m) {
        calc_bernulli(m);
    }
    double get_bernulli(int i) {
        return bernulli[i];
    }
};

int main() {
    int n = 10;
    Todd t = Todd(n);
    for (int i = 0; i < n + 1; ++i) {
        std::cout << t.get_bernulli(i) << std::endl;
    }
    return 0;
}